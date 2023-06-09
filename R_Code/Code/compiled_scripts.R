#### Analysis of food webs
#### Comparison of inferred and empirical webs
#### Jack Shaw

library(tidyverse)
library(ggplot2)
library(igraph)

source("Code/pfim_scripts.R")

theme_set(theme_bw())

#### Load data----

# Select which webs you plan to analyse (default = all)
selected_web_names <- tolower(c(
  "Ythanjacob",
  "Stmarks",
  "Kongsfjorden",
  "Weddell",
  "Chengjiang",
  "Burgess"
))

websizeorder <- tolower(c("Ythanjacob", "Stmarks", "Kongsfjorden", "Weddell","Chengjiang","Burgess"))
modern_webs <- tolower(c("Ythanjacob", "Stmarks", "Kongsfjorden", "Weddell","Chengjiang","Burgess"))


#### Import OG webs----

#Remove trailing white space
og_kon <- read_csv("Data/cirtwill_kongsfjorden_edgelist.csv") %>% rename("res_og_node_name" = res_sp, "con_og_node_name" = con_sp) %>% mutate(across(where(is.character), str_trim)) %>% mutate(across(where(is.character), tolower)) 
og_stm <- read_csv("Data/cirtwill_stmarks_edgelist.csv") %>% rename("res_og_node_name" = res_sp, "con_og_node_name" = con_sp) %>% mutate(across(where(is.character), str_trim))%>% mutate(across(where(is.character), tolower)) 
og_yth <- read_csv("Data/cirtwill_ythan_edgelist.csv") %>% rename("res_og_node_name" = res_sp, "con_og_node_name" = con_sp) %>% mutate(across(where(is.character), str_trim))%>% mutate(across(where(is.character), tolower)) 
og_wed <- read_csv("Data/shaw_weddell_edgelist.csv") %>% dplyr::select(res_og_node_name, con_og_node_name) %>% mutate(across(where(is.character), str_trim))%>% mutate(across(where(is.character), tolower)) 

og_bur <- read_csv("Data/burgess_dunne08_exported_221128.csv") %>% dplyr::select(res_name,con_name) %>% rename("res_og_node_name" = res_name, "con_og_node_name" = con_name) %>% mutate(across(where(is.character), str_trim))
og_che <- read_csv("Data/chengjiang_dunne08_exported_221128.csv") %>% dplyr::select(res_name,con_name) %>% rename("res_og_node_name" = res_name, "con_og_node_name" = con_name) %>% mutate(across(where(is.character), str_trim))

web_list_og <- list(
  kongsfjorden = og_kon,
  stmarks = og_stm,
  ythanjacob = og_yth,
  weddell = og_wed,
  burgess = og_bur,
  chengjiang = og_che
)

# Filter for the webs you're interested in
web_list_og <- web_list_og[c(selected_web_names)]

# Make list of all names to match with metadata
all_web_names <- c()
for (i in names(web_list_og)) {
  temp <- web_list_og[[i]]
  temp <- unique(as.character(as.matrix(temp)))
  all_web_names <- rbind(all_web_names, cbind(node = temp, web = i))
}

#Edit for burgess and chengjiang
all_web_names<-as.data.frame(all_web_names)
all_web_names$node<-ifelse(all_web_names$web %in% c("burgess","chengjiang"),tolower(gsub("[[:punct:][:blank:]]+", " ", all_web_names$node)),all_web_names$node)


#### Import latest metadata----

meta <- read.csv("Data/fw_metadata.csv") %>%
  mutate(across(where(is.character), str_trim)) %>%
  mutate(across(where(is.character), tolower)) %>%
  mutate(web = tolower(web)) %>%
  mutate(
    newest_node_name = ifelse(unchanged_node_name %in% "", updated_node_name, unchanged_node_name),
    oldest_node_name = ifelse(unchanged_node_name %in% "", old_node_name, unchanged_node_name)
  ) %>%
  group_by(newest_node_name) %>% 
  rename("size_bodymass" = og_body_weight)

meta <- meta %>%
  # Filter for the webs you're interested in
  filter(web %in% selected_web_names) %>%
  # Record feeding categories to simplify
  mutate(feeding = gsub("f_mini", "f_surf", feeding))
meta$kingdom[is.na(meta$kingdom)] <- c("kingdom_missing")

# Check that no web-taxa are duplicated, which would mess up merging
# Identify duplicates
id <- paste(meta$web, meta$oldest_node_name)
print(paste(id[id %in% id[duplicated(id)]]))

#Edit for burgess and chengjiang
meta$newest_node_name<-ifelse(meta$web %in% c("burgess","chengjiang"),tolower(gsub("[[:punct:][:blank:]]+", " ", meta$newest_node_name)),meta$newest_node_name)
meta$oldest_node_name<-ifelse(meta$web %in% c("burgess","chengjiang"),tolower(gsub("[[:punct:][:blank:]]+", " ", meta$oldest_node_name)),meta$oldest_node_name)

## Note which taxa need to be dropped
meta <- meta %>%
  # Drop if not recorded in OG web
  filter(tolower(oldest_node_name) %in% tolower(all_web_names$node)) %>%
  # Metazoa not resolved to genus or species
  mutate(drop_taxon_rank = ifelse((kingdom == "Animalia" & is.na(genus)), 1, 0)) %>%
  # Metazoa without functional data
  mutate(drop_functional = ifelse((kingdom == "Animalia" &
    (is.na(tiering) | is.na(motility) | is.na(feeding))), 1, 0)) %>%
  # Metazoa without size info
  mutate(size_bodymass = gsub("NULL", NA, size_bodymass)) %>%
  mutate(drop_size1 = ifelse((kingdom == "Animalia" & is.na(size_max) & is.na(size_bodymass)), 1, 0)) %>%
  mutate(drop_size2 = ifelse((kingdom == "Animalia" & ("unknown" %in% c(size_bodymass, size_max))), 1, 0)) %>%
  mutate(size_select = as.numeric(as.character(gsub("NA", "", paste(size_bodymass, size_max, sep = ""))))) %>%
  # Summarize drops
  mutate(drop_summary = drop_taxon_rank + drop_functional + drop_size1 + drop_size2) %>%
  dplyr::select(-drop_taxon_rank, -drop_functional, -drop_size1, -drop_size2)


meta <- meta %>%
  # Replace taxa with size=0 with the minimum value
  mutate(size_select = ifelse(size_select == 0, min(meta$size_select[meta$size_select > 0], na.rm = TRUE), size_select)) %>%
  # Replace size of matter with minimum
  mutate(size_select = ifelse(kingdom == "Matter", min(meta$size_select[meta$size_select > 0], na.rm = TRUE), size_select))



#### Make standardized food webs----

web_list_updated_taxonomy <- list()
web_list_cleaned_taxonomy <- list()
web_list_cleaned_metazoa <- list()
web_list_fc_base_1 <- list()
web_list_fc_base_2 <- list()

web_list_biglist <- list()

for (i in names(web_list_og)) {

  i_new <- word(i, 1, sep = "_")

  select_web <- web_list_og[[i]] %>% distinct()

  select_meta <- meta %>%
    filter(web == i) %>%
    dplyr::select(newest_node_name, oldest_node_name, kingdom) %>%
    distinct()

  updated_web <- select_web %>%
    left_join(select_meta, by = c("res_og_node_name" = "oldest_node_name")) %>%
    rename("res_new_node_name" = newest_node_name, "res_kingdom" = kingdom) %>%
    left_join(select_meta, by = c("con_og_node_name" = "oldest_node_name")) %>%
    rename("con_new_node_name" = newest_node_name, "con_kingdom" = kingdom)

  web_list_updated_taxonomy[[i]] <- as.matrix(updated_web %>% dplyr::select(res_new_node_name, con_new_node_name))

  # Drop uncoded metazoans
  select_meta_drop <- meta %>%
    filter(web == i) %>%
    filter(drop_summary == 0) %>%
    pull(newest_node_name)
  cleaned_web <- updated_web %>%
    filter(res_new_node_name %in% select_meta_drop & con_new_node_name %in% select_meta_drop)

  web_list_cleaned_taxonomy[[i]] <- as.matrix(cleaned_web %>% dplyr::select(res_new_node_name, con_new_node_name))
  web_list_cleaned_metazoa[[i]] <- as.matrix(cleaned_web %>% filter(res_kingdom == "Animalia" & con_kingdom == "Animalia") %>% dplyr::select(res_new_node_name, con_new_node_name))

  # Group non-metazoans by kingdom
  fc_web <- cleaned_web %>%
    mutate(
      res_node_node_name_fc1 = ifelse(res_kingdom == "Animalia", res_new_node_name, res_kingdom),
      con_node_node_name_fc1 = ifelse(con_kingdom == "Animalia", con_new_node_name, con_kingdom)
    )

  web_list_fc_base_1[[i]] <- as.matrix(fc_web %>% dplyr::select(res_node_node_name_fc1, con_node_node_name_fc1) %>% distinct())

  # Group non-metazoans as single node
  fc_web <- fc_web %>%
    mutate(
      res_node_node_name_fc2 = ifelse(res_kingdom == "Animalia", res_new_node_name, "Non-metazoan"),
      con_node_node_name_fc2 = ifelse(con_kingdom == "Animalia", con_new_node_name, "Non-metazoan")
    )

  web_list_fc_base_2[[i]] <- as.matrix(fc_web[, c("res_node_node_name_fc2", "con_node_node_name_fc2")] %>% distinct())

  test <- (as.data.frame((fc_web %>% dplyr::select(res_node_node_name_fc2, con_node_node_name_fc2))) %>% distinct())


  web_list_biglist[[i]] <- fc_web

  print(i)
}


mstr_web_list <- list(
  ut = web_list_updated_taxonomy,
  ct = web_list_cleaned_taxonomy,
  mo = web_list_cleaned_metazoa,
  fc1 = web_list_fc_base_1,
  fc2 = web_list_fc_base_2
)



#### Description of data structure----

# mstr_web_list<-list(ut=web_list_updated_taxonomy, 
#                     og=web_list_og, 
#                     ct=web_list_cleaned_taxonomy, #Cleaned taxonomy (dropping taxa without known functional information)
#                     mo=web_list_cleaned_metazoa, #Metazoan only
#                     fc1=web_list_fc_base_1, #Group non-metazoans by kingdom
#                     fc2=web_list_fc_base_2) #Group non-metazoans as single node


#### Variables to set up----
load("Code/pfim_metrics.rda")
var_names <- pfwim_metrics %>%
  filter(recommended %in% c(1, 2)) %>%
  # filter(recommended %in% c(1)) %>%
  dplyr::select(variable_short_name, variable_long_name)

# Import food web
trait_combos <- read.csv("Data/trait_categories.csv") %>% filter(trait_type_resource != "sizegp" & trait_type_consumer != "sizegp")
trait_defs <- read.csv("Data/trait_definitions.csv") %>%
  arrange(order) %>%
  mutate(full_trait = paste(order, full_trait, sep = ""))

fig_names <- tibble::tribble(
  ~old, ~new,
  "C", "C",
  "Size", "Size",
  "Diameter", "Diameter",
  "SD generality", "SD gen.",
  "SD vulnerability", "SD vul.",
  "Normalized generality", "Generality",
  "Normalized vulnerability", "Vulnerability",
  "Norm. Btw.", "Btw.",
  "Mean norm. btw.", "Mean btw.",
  "Mean norm. deg.", "Mean deg.",
  "TL (std)", "TL",
  "Mean TL (std)", "Mean TL",
  #"Max TL (std)", "Max TL",
  "OI (std)", "OI",
  #"SOI (std)", "SOI",
  "Mean CPL", "Mean CPL",
  "Norm. mot: Linear chains", "Motif: Lin. ch.",
  "Norm. mot: Omnivory", "Motif: Omn.",
  "Norm. mot: App. comp.", "Motif: Ap. comp.",
  "Norm. mot: Dir. comp.", "Motif: Dir. comp."
)
fig_web_names <- tibble::tribble(
  ~old, ~new,
  "ythanjacob", "Ythan",
  "stmarks", "St Marks",
  "kongsfjorden", "Kongsfjorden",
  "loughhyne", "Lough Hyne",
  "weddell", "Weddell",
  "chengjiang","Chengjiang",
  "burgess","Burgess"
)


#### BASIC STATS: Network level statistics----

stats_netlevel_mstr <- c()
for (i in names(mstr_web_list)) {
  select_list <- mstr_web_list[[i]]

  for (j in names(select_list)) {
    select_fw <- as.matrix(select_list[[j]])
    select_fw <- unique(select_fw[, ])
    ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i)
    stats_netlevel_mstr <- bind_rows(stats_netlevel_mstr, ccc)
  }

  print(i)
}

stats_netlevel <- stats_netlevel_mstr %>%
  dplyr::rename("variable" = "ind", "value" = "values") %>%
  left_join(var_names, by = c("variable" = "variable_short_name")) %>%
  mutate(web = factor(web, levels = websizeorder))


# Table to change names of shortened model webs
abb_names <- tibble::tribble(
  ~short_name, ~long_name,
  "ut", "Original published",
  "ct", "Clean published",
  "fc1", "Multiple basal nodes",
  "fc2", "Single basal node",
  "mo", "Metazoan-only"
)
##SUPFIGURE1##
ggplot(stats_netlevel %>%
  mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = jlookup(web, fig_web_names, "old", "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  filter(!is.na(variable_long_name) & variable_long_name != "Mean deg.") %>%
  mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
  mutate(type = factor(type, levels = c("ut", "ct", "fc1", "fc2", "mo"))) %>%
  mutate(type = jlookup(type, abb_names, matching_col = "short_name", new_values = "long_name")) %>%
  mutate(type = factor(type, levels = abb_names$long_name))) +
  geom_point(aes(x = web, y = value, color = type, shape = type)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable_long_name, scales = "free", ncol = 5)


#### BASIC STATS: Niche normalized----
stats_netlevel_niche_SC <- stats_netlevel %>%
  filter(variable_long_name %in% c("Size", "C")) %>% # filter(type %in% c("fc2","mo","ut") & web %in% modern_webs) %>%
  dplyr::select(-variable) %>%
  pivot_wider(names_from = variable_long_name, values_from = value)



library(foreach)
library(doParallel)
library(R.utils)

print(Sys.time())
argen <- c()
closeAllConnections()
myCluster <- makeCluster(4)
registerDoParallel(myCluster)
stats_netlevel_niche_stats <- c()
reps <- 5
argen <- foreach(i = 1:nrow(stats_netlevel_niche_SC), .combine = rbind, .packages = c("igraph", "tidyverse"), .errorhandling = "remove") %dopar% {
  stats_netlevel_niche_stats <- c()
  temp_df <- niche_model_replicates(stats_netlevel_niche_SC[i, c("Size")] %>% pull(),
    stats_netlevel_niche_SC[i, c("C")] %>% pull(),
    reps = reps,
    return_mean = TRUE
  )
  stats_netlevel_niche_stats <- cbind(temp_df,
    web = stats_netlevel_niche_SC[i, c("web")] %>% pull(),
    type = stats_netlevel_niche_SC[i, c("type")] %>% pull()
  )
  as.data.frame(stats_netlevel_niche_stats)
  stats_netlevel_niche_stats
}

stopCluster(myCluster)
closeAllConnections()
registerDoSEQ()
print(Sys.time())

stats_netlevel_niche_comparison <- stats_netlevel %>% left_join(argen, by = c("variable", "type", "web"))
stats_netlevel_niche_comparison <- as.data.frame(as.matrix(stats_netlevel_niche_comparison)) %>% mutate_at(vars(value, min, med, max), as.numeric)

me_vals <- c()
for (i in 1:nrow(stats_netlevel_niche_comparison)) {
  temp_me <- calc_ME(
    value = stats_netlevel_niche_comparison[i, c("value")],
    model_lower = stats_netlevel_niche_comparison[i, c("min")],
    model_median = stats_netlevel_niche_comparison[i, c("med")],
    model_upper = stats_netlevel_niche_comparison[i, c("max")]
  )
  me_vals <- c(me_vals, temp_me)
}

stats_netlevel_niche_comparison_joined <- cbind(stats_netlevel_niche_comparison, model_error = me_vals) %>% mutate(web = factor(web, levels = websizeorder))

##SUPFIGURE2##
ggplot(
  stats_netlevel_niche_comparison_joined %>%
    mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
    mutate(web = jlookup(web, fig_web_names, "old", "new")) %>%
    mutate(web = factor(web, levels = fig_web_names$new)) %>%
    filter(!is.na(variable_long_name) & variable_long_name != "Mean deg.") %>%
    mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
    mutate(type = factor(type, levels = c("ut", "ct", "fc1", "fc2", "mo"))) %>%
    mutate(type = jlookup(type, abb_names, matching_col = "short_name", new_values = "long_name")) %>%
    mutate(type = factor(type, levels = abb_names$long_name)) %>%
    # mutate(web=factor(web, levels=str_to_title(websizeorder))) %>%
    filter(!is.na(variable_long_name) & !is.na(model_error)),
  aes(x = web, y = model_error, color = type, shape = type)
) +
  geom_hline(yintercept = 1, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkgrey", linetype = "dashed") +
  geom_point() +
  xlab("") +
  ylab("Effect size") +
  theme(axis.text.x = element_text(
    angle = 45,
    # face=websizefont,
    hjust = 1
  ), legend.title = element_blank()) +
  facet_wrap(. ~ variable_long_name, scales = "free", nrow = 3)


#### BASIC STATS: Node-level statistics----
stats_nodelevel_mstr <- c()
for (i in names(mstr_web_list)) {
  sel_web <- mstr_web_list[[i]]

  for (j in names(sel_web)) {
    sel_web2 <- graph_from_edgelist(as.matrix(sel_web[[j]]), directed = TRUE)
    temp_node_stats <- calc_node_stats(sel_web2)

    temp_node_list <- cbind(temp_node_stats,
      web = j,
      type = i
    )

    stats_nodelevel_mstr <- rbind(stats_nodelevel_mstr, temp_node_list)
  }

  print(i)
}

stats_nodelevel <- as.data.frame(stats_nodelevel_mstr) %>%
  pivot_longer(cols = !c("taxon", "web", "type"), names_to = "variable", values_to = "value") %>%
  mutate(value = as.numeric(as.character(value))) %>%
  filter(variable != "size") %>%
  mutate(web = factor(web, levels = websizeorder)) %>%
  mutate(type = factor(type, levels = c("ut", "ct", "fc1", "fc2", "mo")))

ggplot(stats_nodelevel, aes(x = web, y = value)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(variable ~ type, scales = "free")

##SUPFIGURE3##
ggplot(stats_nodelevel %>%
  left_join(var_names, by = c("variable" = "variable_short_name")) %>%
  mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(!is.na(variable_long_name)) %>%
  mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
  mutate(type = factor(type, levels = c("ut", "ct", "fc1", "fc2", "mo"))) %>%
  mutate(type = jlookup(type, abb_names, matching_col = "short_name", new_values = "long_name")) %>%
  mutate(type = factor(type, levels = abb_names$long_name)) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  filter(!is.na(variable_long_name)), aes(x = type, y = value)) +
  geom_boxplot(show.legend = FALSE, aes(color = type)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable_long_name ~ web, scales = "free")

ggplot(stats_nodelevel %>% filter(type %in% c("fc2", "ct")) %>% group_by(variable, web, type) %>% summarize(avg = mean(value, na.rm = TRUE)), aes(x = web, y = avg, color = type)) +
  geom_point() +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1
  ), legend.title = element_blank()) +
  facet_wrap(variable ~ ., scales = "free")


ggplot(stats_nodelevel %>% filter(type %in% c("fc2", "ct")), aes(x = type, y = value, color = type)) +
  geom_boxplot() +
  geom_point(data = stats_nodelevel %>% filter(type %in% c("fc2", "ct")) %>% group_by(variable, web, type) %>% summarize(avg = mean(value, na.rm = TRUE)), aes(x = type, y = avg), color = "black") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1
  ), legend.title = element_blank()) +
  facet_wrap(web ~ variable, scales = "free")


# Generate simplified boxplot stats
stats_nodelevel_summarized <- stats_nodelevel %>%
  group_by(variable, web, type) %>%
  summarize(
    avg = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.25, na.rm = TRUE),
    mid = quantile(value, 0.5, na.rm = TRUE),
    upper = quantile(value, 0.75, na.rm = TRUE)
  ) %>%
  left_join(var_names, by = c("variable" = "variable_short_name")) # %>% filter(!is.na(variable_long_name))


pd <- position_dodge(0.2)
ggplot(stats_nodelevel_summarized %>% filter(type %in% c("fc2", "ct")), aes(x = web, color = type)) +
  geom_point(aes(y = mid), position = pd, size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = pd, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  xlab("") +
  ylab("") +
  facet_wrap(. ~ variable_long_name, scales = "free")


#### LIFE HABIT STATS: Metadata----

vars_for_cor <- c("web", "tiering", "feeding", "motility", "size_select")

all_metadata <- meta %>%
  # filter(web!="sanak")
  filter(drop_summary == 0 & kingdom == "Animalia") %>%
  dplyr::select(c(vars_for_cor)) %>%
  # distinct() %>%
  mutate_all(., as.factor) %>%
  mutate_at(c("size_select"), as.character) %>%
  mutate_at(c("size_select"), as.numeric) %>%
  mutate(size_select = ifelse(size_select <= 0.000000000001, 0.000000000001, size_select))

# Check that there aren't duplicate taxa within webs
temp <- all_metadata %>%
  group_by(web, newest_node_name) %>%
  summarise(count = n_distinct(newest_node_name)) %>%
  filter(count > 1)
nrow(temp) # should be 0



#### LIFE HABIT STATS: Functional diversity statistics----

stats_functional_mstr <- c()
for (i in unique(all_metadata$web)) {
  temp <- as.data.frame(all_metadata) %>%
    filter(web == i) %>%
    dplyr::select(-web, -newest_node_name) %>%
    select_if(colSums(!is.na(.)) > 0) %>%
    drop_na() # Drop NA columns and columns
  fd_sel <- stack(calc_functional_stats(temp))

  stats_functional_mstr <- bind_rows(stats_functional_mstr, cbind(fd_sel, web = i))
  print(i)
}

stats_functional <- stats_functional_mstr %>%
  dplyr::rename("variable" = "ind", "value" = "values") %>%
  mutate(web = factor(web, levels = websizeorder)) %>%
  drop_na(value)
stats_functional$variable <- jlookup(stats_functional$variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")

##SUPFIGURE8##
ggplot(stats_functional %>% filter(variable %in% c(
  "Functional dispersion", "Functional divergence", "Functional evenness",
  "Trait combinations"
)) %>% mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new))) +
  geom_point(aes(x = web, y = value)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1
  ), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free")





#### FUNCTIONAL STRUCTURE: Compare life habits and node positions----


stats_nodehabits <- as.data.frame(stats_nodelevel_mstr) %>%
  filter(type == "fc2" & taxon != "Non-metazoan") %>%
  left_join(all_metadata, by = c("web", "taxon" = "newest_node_name")) %>%
  distinct() %>%
  pivot_longer(c("tiering", "feeding", "motility"), names_to = c("trait"), values_to = c("trait_value")) %>%
  pivot_longer(c("norm_btw", "norm_degree_in", "norm_degree_out", "tl_std", "tl_sw", "oi_std", "oi_sw"), names_to = c("variable"), values_to = c("variable_value")) %>%
  mutate_all(., as.character) %>%
  mutate_at(vars(variable_value, size_select), as.numeric) %>%
  mutate(web = factor(web, levels = websizeorder))
stats_nodehabits$trait_value <- jlookup(stats_nodehabits$trait_value, trait_defs, matching_col = "short_trait", new_values = "full_trait")
stats_nodehabits$variable <- jlookup(stats_nodehabits$variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")

# Size correlations
ggplot(stats_nodehabits) +
  geom_density(aes(x = log(size_select), color = web))

ggplot(
  stats_nodehabits %>% dplyr::select(taxon, web, variable_value, size_select, variable) %>%
    distinct() %>%
    mutate(size_type = ifelse(web %in% c("burgess", "chengjiang"), "length", "body mass")) %>%
    filter(!variable %in% c("OI (sw)", "TL (sw)")) %>%
    mutate(web = factor(web, levels = websizeorder)),
  aes(x = log(size_select), y = variable_value, color = web)
) +
  geom_smooth(method = "lm") +
  ylab("Value") +
  xlab("Log body mass (UNIT)") +
  facet_grid(variable ~ ., scales = "free")

stats_nodehabits_webmean <- stats_nodehabits %>%
  group_by(web, trait, trait_value, variable) %>%
  summarise(
    mean_numeric_vals = mean(variable_value, na.rm = T),
    n_taxon = n_distinct(taxon)
  ) %>%
  group_by(web, trait, variable) %>%
  mutate(prop_taxon = n_taxon / sum(n_taxon)) %>%
  distinct()

stats_nodehabits_allmean <- stats_nodehabits %>%
  group_by(trait, trait_value, variable) %>%
  summarise(mean_numeric_vals = mean(variable_value, na.rm = T)) %>%
  distinct()

ggplot(
  stats_nodehabits_webmean %>% filter(!variable %in% c("OI (sw)", "TL (sw)")),
  aes(x = trait_value, y = mean_numeric_vals)
) +
  xlab("Life habit") +
  ylab("Mean value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(aes(fill = web, color = web, group = web)) +
  geom_line(aes(color = web, group = web)) +  scale_x_discrete() +
  facet_grid(variable ~ trait, scales = "free")

stats_nodehabits$trait_value <- jlookup(stats_nodehabits$trait_value, trait_defs, matching_col = "short_trait", new_values = "full_trait")


##SUPFIGURE7##
p1 <- ggplot(
  stats_nodehabits %>% dplyr::select(taxon, web, variable_value, size_select, variable) %>%
    distinct() %>%
    mutate(size_type = ifelse(web %in% c("burgess", "chengjiang"), "length", "body mass")) %>%
    filter(!variable %in% c("OI (sw)", "TL (sw)")) %>%
    mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
    mutate(variable = factor(variable, levels = fig_names$new)) %>%
    mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
    mutate(web = factor(web, levels = fig_web_names$new)),
  aes(x = log(size_select), y = variable_value, color = web)
) +
  geom_smooth(method = "lm", show.legend = F) +
  ylab("Value") +
  xlab("Log body mass (g)") +
  facet_grid(variable ~ "Size", scales = "free")

p2 <- ggplot(
  stats_nodehabits_webmean %>% filter(!variable %in% c("OI (sw)", "TL (sw)")) %>%
    mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
    mutate(trait = str_to_title(trait)) %>%
    mutate(variable = factor(variable, levels = fig_names$new)) %>%
    mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
    mutate(web = factor(web, levels = fig_web_names$new)),
  aes(x = trait_value, y = mean_numeric_vals)
) +
  xlab("Life habit") +
  ylab("Mean value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  geom_point(aes(fill = web, color = web, group = web)) +
  geom_line(aes(color = web, group = web)) +
  scale_x_discrete() +
  facet_grid(variable ~ trait, scales = "free", space = "free_x")

library(patchwork)
(p1 | p2) + plot_layout(guides = "collect", widths = c(1, 3))




#### PREDICTING INTERACTIONS: interactions between paired traits----


# Generate all possible versus feasible interactions
selected_web_list <- mstr_web_list[["fc2"]]
all_possible_interactions <- c()
for (i in names(selected_web_list)) {
  selected_web <- as.data.frame(selected_web_list[[i]])
  selected_web <- filter_all(selected_web, all_vars(. != "Non-metazoan"))

  selected_list <- as.data.frame(mat2list(as.matrix(as_adjacency_matrix(graph_from_edgelist(as.matrix(selected_web), directed = TRUE)))))
  colnames(selected_list) <- c("consumer", "resource", "interaction")
  selected_list$interaction <- ifelse(selected_list$interaction > 1, 1, selected_list$interaction)
  all_possible_interactions <- rbind(all_possible_interactions, cbind(selected_list, web = i))
}

# Join with metadata
all_possible_interactions_meta <- all_possible_interactions %>%
  left_join(all_metadata, by = c("consumer" = "newest_node_name", "web")) %>%
  distinct() %>%
  left_join(all_metadata, by = c("resource" = "newest_node_name", "web"), suffix = c("_con", "_res")) %>%
  distinct() %>%
  mutate_if(is.character, factor) %>%
  mutate(interaction = as.factor(as.integer(as.logical(interaction))))
# Check if any NA rows
which(is.na(all_possible_interactions_meta), arr.ind = TRUE)


# Paired proportions
stats_pairedprops <- c()
for (i in c("feeding", "motility", "tiering")) {
  for (j in unique(all_possible_interactions_meta$web)) {
    stats_pairedtraits_sub <- all_possible_interactions_meta %>%
      filter(web == j) %>%
      dplyr::select(interaction, contains(i)) %>%
      drop_na(contains(i)) %>% # use this to drop taxa with missing vals
      mutate(interaction = as.numeric(as.character(interaction))) %>%
      group_by_at(setdiff(colnames(.), "interaction")) %>%
      mutate(
        pres = sum(interaction, na.rm = T),
        total = n()
      ) %>%
      mutate(abs = total - pres, proportion = pres / total) %>%
      dplyr::select(-interaction) %>%
      distinct() %>%
      mutate(web = j, trait = i)

    colnames(stats_pairedtraits_sub) <- c("con_trait", "res_trait", "pres", "total", "abs", "proportion", "web", "trait")

    stats_pairedprops <- rbind(stats_pairedprops, stats_pairedtraits_sub)
  }

  print(i)
}


#### PREDICTING INTERACTIONS: Predator prey size ratios----

predprey_sizes <- all_possible_interactions_meta %>%
  dplyr::select(
    consumer, resource, web, interaction,
    size_select_con, size_select_res
  ) %>%
  mutate_at(vars(size_select_con, size_select_res), as.numeric) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  distinct()

ggplot(data = predprey_sizes, aes(x = log(size_select_res), y = log(size_select_con), color = interaction)) +
  geom_point(data = predprey_sizes %>% filter(interaction == 1), alpha = 0.2, size = 0.5) +
  geom_smooth(data = predprey_sizes, method = "lm", se = TRUE) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Log prey size") +
  ylab("Log predator size") +
  facet_wrap(web ~ .)

##SUPFIGURE9##
ggplot(data = predprey_sizes %>% filter(interaction == 1), aes(x = log(size_select_con), y = log(size_select_res))) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color="darkgrey") +
  geom_point(alpha = 0.2, size = 0.5) +
  xlab("Log predator size") +
  ylab("Log prey size") +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  facet_wrap(. ~ (web),nrow=1)

ggplot(predprey_sizes %>% filter(interaction == 1) %>% mutate(size_select_ratio = log(size_select_res) / log(size_select_con)), aes(x = size_select_ratio, color = web)) +
  geom_density() +
  coord_cartesian(xlim = c(-10, 10))

predprey_sizes_ratios <- predprey_sizes %>%
  filter(interaction == 1) %>%
  group_by(web) %>%
  mutate(size_select_res = ifelse(size_select_res == 0, 0.001, size_select_res)) %>%
  mutate(size_select_con = ifelse(size_select_con == 0, 0.001, size_select_con)) %>%
  do(size_ratio = broom::tidy(lm(log(size_select_con) ~ log(size_select_res), data = .))) %>%
  unnest(size_ratio)

# Our range of size ratios is much wider than Brose
shaw_size_rule <- function(res_size, con_size) {
  ifelse(log10(con_size / res_size) > (min(size_vector$stats)) & log10(con_size / res_size) < (max(size_vector$stats)), 1, 0)
}
# As per marine clades in Brose et al 2006
brose_size_rule <- function(res_size, con_size) {
  ifelse(log10(con_size / res_size) > -2 & log10(con_size / res_size) < 5, 1, 0)
}

#### INFER FOOD WEB: Infer FWs based on functional data----


# Choose variables for inference
selected_vars_categorical <- c("motility", "feeding", "tiering")
selected_vars_numerical <- c("size_select") # If size is categorical, then enter NULL

# Use this for a single basal node
basal_node <- data.frame(
  stringsAsFactors = FALSE,
  taxon = c("BASAL NODE", "BASAL NODE", "BASAL NODE"),
  tiering = c("t_pela", "t_surf", "t_shin"),
  motility = c("m_nmun", "m_nmun", "m_nmun"),
  feeding = c("f_prim", "f_prim", "f_prim"),
  size_select = c(1e-12, 1e-12, 1e-12)
)

# Set up lists to store inferred metawebs and species webs
web_list_inferred_mw <- list()
web_list_inferred_species <- list()


# Select number of runs for species webs
webs_to_run <- c("ythanjacob", "stmarks", "kongsfjorden", "weddell","chengjiang","burgess")

runs <- 5
for (i in webs_to_run) {

  # Only need to predict for metazoans
  data <- meta %>%
    filter(web == i & kingdom == "Animalia" & drop_summary == 0) %>%
    dplyr::select(newest_node_name, all_of(selected_vars_categorical), all_of(selected_vars_numerical)) %>%
    mutate(size_select = ifelse(size_select <= 0.0000000001, 0.0000000001, size_select)) %>%
    drop_na() %>%
    group_by(newest_node_name) %>%
    sample_n(1) %>%
    rename("taxon" = newest_node_name) %>%
    mutate_at(all_of(selected_vars_categorical), factor) %>%
    mutate_at(c(selected_vars_numerical), as.numeric) %>%
    drop_na() %>%
    bind_rows(basal_node)

  inf_mw <- infer_edgelist(data, col_taxon = "taxon", col_num_size = selected_vars_numerical, return_full_matrix = FALSE, print_dropped_taxa = FALSE, cat_combo_list = trait_combos)

  # Remove cannibalistic links (necessary if applying the link distributions below)
  inf_mw <- as.matrix(as.data.frame(inf_mw) %>% filter(taxon_resource != taxon_consumer))

  web_list_inferred_mw[[i]] <- inf_mw

  web_list_inferred_species[[i]] <- powerlaw_prey(inf_mw, n_samp = runs)

  print(i)
}




#### COMPARE INFERRED WEBS: Compare metaweb and realized webs----

compare_inferred_list <- list(
  inferred_mw = web_list_inferred_mw,
  fc2 = mstr_web_list[["fc2"]]
)

stats_compare_inferred <- c()
for (i in names(compare_inferred_list)) {
  select_list <- compare_inferred_list[[i]]

  for (j in names(select_list)) {
    select_list2 <- select_list[[j]]

    if (is.matrix(select_list2)) {
      select_fw <- as.matrix(select_list2)
      select_fw <- unique(select_fw[, ])
      ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = 0)
      stats_compare_inferred <- bind_rows(stats_compare_inferred, ccc)
    } else {
      for (k in 1:length(select_list2)) {
        select_list3 <- select_list2[[k]]

        select_fw <- as.matrix(select_list3)
        select_fw <- unique(select_fw[, ])
        ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = k)
        stats_compare_inferred <- bind_rows(stats_compare_inferred, ccc)
      }
    }
  }

  print(i)
}

stats_compare_inferred2 <- stats_compare_inferred %>%
  dplyr::rename("variable" = "ind", "value" = "values") %>%
  filter(!variable %in% c("mean_betweenness", "soi_sw", "mean_oi_std", "mean_oi_sw", "mean_tl_sw", "mot_lin", "mot_omn", "mot_ap_comp", "mot_dir_comp")) %>%
  mutate(web = factor(web, levels = websizeorder))
stats_compare_inferred2$variable <- jlookup(stats_compare_inferred2$variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")
stats_compare_inferred2 <- stats_compare_inferred2 %>% filter(!is.na(variable))

ggplot(stats_compare_inferred2 %>% filter(type %in% c("fc2", "inferred_mw_fd", "inferred_mw_el", "inferred_mw"))) +
  geom_point(aes(x = web, y = value, color = type, shape = type)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free")

#### FIGURE
ggplot(stats_compare_inferred2 %>%
  filter(type %in% c("fc2", "inferred_mw")) %>%
  mutate(type = recode(type, fc2 = "Empirical web", inferred_mw = "Inferred metaweb")) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  group_by(type, variable, web) %>%
  mutate(avg = mean(value, na.rm = TRUE)) %>%
  mutate(variable = factor(variable, levels = fig_names$new))) +
  geom_point(aes(x = web, y = avg, color = type, shape = type)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free", ncol = 5)



#### COMPARE INFERRED WEBS: Compare inferred and FC2 niche model error----

stats_compare_niche_std <- stats_compare_inferred2 %>%
  filter(type %in% c("fc2", "inferred_mw")) %>%
  filter(variable %in% c("Size", "C")) %>%
  dplyr::select(-replicate) %>%
  pivot_wider(names_from = variable, values_from = value)

library(foreach)
library(doParallel)
library(R.utils)

print(Sys.time())
argen <- c()
closeAllConnections()
myCluster <- makeCluster(8)
registerDoParallel(myCluster)
stats_netlevel_niche_stats <- c()
# This really needs to be higher (i.e., 1000)
reps <- 5

argen <- foreach(i = 1:nrow(stats_compare_niche_std), .combine = rbind, .packages = c("igraph", "tidyverse"), .errorhandling = "remove") %dopar% {
  stats_netlevel_niche_stats <- c()
  temp_df <- niche_model_replicates(stats_compare_niche_std[i, c("Size")] %>% pull(),
    stats_compare_niche_std[i, c("C")] %>% pull(),
    reps = reps,
    return_mean = TRUE
  )
  stats_netlevel_niche_stats <- cbind(temp_df,
    web = stats_compare_niche_std[i, c("web")] %>% pull(),
    type = stats_compare_niche_std[i, c("type")] %>% pull()
  )
  as.data.frame(stats_netlevel_niche_stats)
  stats_netlevel_niche_stats
}

stopCluster(myCluster)
closeAllConnections()
registerDoSEQ()
print(Sys.time())

argen$variable_short_name <- argen$variable
argen$variable <- jlookup(argen$variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")
argen$variable_short_name <- NULL

stats_inferred_niche_comparison <- stats_compare_inferred2 %>%
  filter(type %in% c("fc2", "inferred_mw")) %>%
  left_join(argen, by = c("variable", "type", "web"))
stats_inferred_niche_comparison <- as.data.frame(as.matrix(stats_inferred_niche_comparison)) %>% mutate_at(vars(value, min, med, max), as.numeric)

me_vals <- c()
for (i in 1:nrow(stats_inferred_niche_comparison)) {
  temp_me <- calc_ME(
    value = stats_inferred_niche_comparison[i, c("value")],
    model_lower = stats_inferred_niche_comparison[i, c("min")],
    model_median = stats_inferred_niche_comparison[i, c("med")],
    model_upper = stats_inferred_niche_comparison[i, c("max")]
  )
  me_vals <- c(me_vals, temp_me)
}

stats_inferred_niche_comparison_joined <- cbind(stats_inferred_niche_comparison, model_error = me_vals) %>%
  mutate(web = factor(web, levels = websizeorder)) %>%
  filter(!is.na(model_error)) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new)) %>%
  mutate(type = recode(type, fc2 = "Empirical web", inferred_mw = "Inferred metaweb")) %>%
  mutate(type = factor(type, levels = c("Empirical web", "Inferred metaweb", "Hypothetical realized web")))


#### FIGURE
ggplot(
  stats_inferred_niche_comparison_joined,
  aes(x = web, y = model_error, color = type, shape = type)
) +
  geom_hline(yintercept = 1, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkgrey", linetype = "dashed") +
  geom_point() +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free", ncol = 5)




#### COMPARE INFERRED WEBS: Node positions of taxa in inferred and empirical webs----

compare_inferred_list <- list(
  inferred_mw = web_list_inferred_mw,
  fc2 = mstr_web_list[["fc2"]]
)

compare_inferred_nodes <- c()
for (i in names(compare_inferred_list)) {
  select_list <- compare_inferred_list[[i]]

  for (j in names(select_list)) {
    select_list2 <- select_list[[j]]

    if (is.matrix(select_list2)) {
      select_fw <- as.matrix(select_list2)
      select_fw <- unique(select_fw[, ])
      ccc <- cbind(calc_node_stats(graph_from_edgelist(select_fw, directed = TRUE)), web = j, type = i, replicate = 0)
      compare_inferred_nodes <- bind_rows(compare_inferred_nodes, ccc)
    } else {
      for (k in 1:length(select_list2)) {
        select_list3 <- select_list2[[k]]

        select_fw <- as.matrix(select_list3)
        select_fw <- unique(select_fw[, ])
        ccc <- cbind(calc_node_stats(graph_from_edgelist(select_fw, directed = TRUE)), web = j, type = i, replicate = k)
        compare_inferred_nodes <- bind_rows(compare_inferred_nodes, ccc)
      }
    }
  }

  print(i)
}


compare_inferred_nodes2 <- compare_inferred_nodes %>%
  mutate(type = recode(type, fc2 = "Empirical web", inferred_mw = "Inferred metaweb")) %>%
  # filter(type %in% c("inferred_mw","fc2")) %>%
  mutate(type_replicate = paste(type, replicate, sep = "_")) %>%
  dplyr::select(-type, -replicate) %>%
  pivot_longer(cols = c(-taxon, -web, -type_replicate)) %>%
  pivot_wider(names_from = "type_replicate", values_from = "value") %>%
  replace(is.na(.), 0) %>%
  mutate(web = factor(web, levels = websizeorder)) %>%
  left_join(var_names, by = c("name" = "variable_short_name"))

ggplot(compare_inferred_nodes2) +
  geom_point(aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)) +
  xlab("inferred") +
  ylab("fc2") +
  facet_wrap(name ~ web, scales = "free", ncol = length(unique(compare_inferred_nodes2$web)))

ggplot(compare_inferred_nodes2) +
  geom_point(aes(x = inferred_mw_0, y = fc2_0 - inferred_mw_0)) +
  xlab("") +
  ylab("") +
  facet_wrap(name ~ web, scales = "free", ncol = length(unique(compare_inferred_nodes2$web)))

ggplot(compare_inferred_nodes2 %>% filter(!variable_long_name %in% c("Betweenness", "OI (sw)", "TL (sw)", "Generality", "Vulnerability"))) +
  geom_boxplot(aes(x = web, y = inferred_mw_0 - fc2_0)) +
  xlab("web") +
  ylab("mw - fc2 value") +
  facet_wrap(variable_long_name ~ ., scales = "free", ncol = length(unique(compare_inferred_nodes2$web))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### FIGURE
ggplot(
  compare_inferred_nodes2 %>% filter(!name %in% c("betweenness", "degree_in", "degree_out", "oi_sw", "tl_sw")),
  aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)
) +
  geom_point(color = "black", alpha = 0.8) +
  xlab("Inferred") +
  ylab("Empirical") +
  facet_wrap(variable_long_name ~ web, scales = "free", ncol = length(unique(compare_inferred_nodes2$web)))

compare_inferred_nodes_ttt <- compare_inferred_nodes2 %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
  mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
  filter(variable_long_name != "NA")

p1 <- ggplot(compare_inferred_nodes_ttt %>% filter(variable_long_name == "Btw."), aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)) +
  geom_point(color = "black", alpha = 0.8) +
  xlab("Inferred") +
  ylab("Empirical") +
  facet_grid(variable_long_name ~ web) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p2 <- ggplot(compare_inferred_nodes_ttt %>% filter(variable_long_name == "Generality"), aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)) +
  geom_point(color = "black", alpha = 0.8) +
  xlab("Inferred") +
  ylab("Empirical") +
  facet_grid(variable_long_name ~ web) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.x = element_blank())
p3 <- ggplot(compare_inferred_nodes_ttt %>% filter(variable_long_name == "Vulnerability"), aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)) +
  geom_point(color = "black", alpha = 0.8) +
  xlab("Inferred") +
  ylab("Empirical") +
  facet_grid(variable_long_name ~ web) +
  theme(axis.title.x = element_blank(), strip.text.x = element_blank())
p4 <- ggplot(compare_inferred_nodes_ttt %>% filter(variable_long_name == "OI"), aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)) +
  geom_point(color = "black", alpha = 0.8) +
  xlab("Inferred") +
  ylab("Empirical") +
  facet_grid(variable_long_name ~ web) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.x = element_blank())
p5 <- ggplot(compare_inferred_nodes_ttt %>% filter(variable_long_name == "TL"), aes(x = `Inferred metaweb_0`, y = `Empirical web_0`)) +
  geom_point(color = "black", alpha = 0.8) +
  xlab("Inferred") +
  ylab("Empirical") +
  facet_grid(variable_long_name ~ web) +
  theme(axis.title.y = element_blank(), strip.text.x = element_blank())

(p1 / p2 / p3 / p4 / p5)


#### FIGURE
temp <- compare_inferred_nodes2 %>%
  filter(!variable_long_name %in% c("Betweenness", "OI (sw)", "TL (sw)", "Generality", "Vulnerability")) %>%
  pivot_longer(c("Inferred metaweb_0", "Empirical web_0"), names_to = "valuetype", values_to = "value")

ggplot(temp %>%
  mutate(valuetype = gsub("_0", "", valuetype)) %>%
  mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
  mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  filter(variable_long_name != "NA")) +
  geom_boxplot(aes(x = valuetype, y = value, color = valuetype), alpha = 0.5, show.legend = FALSE) +
  xlab("") +
  ylab("Value") +
  facet_grid(variable_long_name ~ web, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p1 <- ggplot(compare_inferred_nodes2 %>% filter(!variable_long_name %in% c("Betweenness", "OI (sw)", "TL (sw)", "Generality", "Vulnerability"))) +
  geom_boxplot(aes(x = web, y = `Inferred metaweb_0` - `Empirical web_0`)) +
  xlab("web") +
  ylab("Difference btw node position in inferred & empirical webs") +
  facet_grid(variable_long_name ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### COMPARE BASAL TYPES: Type1 v type2----

basal_node_type1 <- data.frame(
  stringsAsFactors = FALSE,
  taxon = c("BASAL NODE", "BASAL NODE", "BASAL NODE"),
  tiering = c("t_pela", "t_surf", "t_shin"),
  motility = c("m_nmun", "m_nmun", "m_nmun"),
  feeding = c("f_prim", "f_prim", "f_prim"),
  size_select = c(1e-08, 1e-08, 1e-08)
)

basal_node_type2 <- data.frame(
  stringsAsFactors = FALSE,
  taxon = c("Suspended", "Surficial", "Infaunal"),
  tiering = c("t_pela", "t_surf", "t_shin"),
  motility = c("m_nmun", "m_nmun", "m_nmun"),
  feeding = c("f_prim", "f_prim", "f_prim"),
  size_select = c(1e-08, 1e-08, 1e-08)
)

# Set up lists to store inferred metawebs and species webs
web_list_inferred_mw_basaltype1 <- list()
web_list_inferred_mw_basaltype2 <- list()

# Select number of runs for species webs
for (i in unique(meta$web)) {

  data1 <- meta %>%
    filter(web == i & kingdom == "Animalia" & drop_summary == 0) %>%
    dplyr::select(newest_node_name, all_of(selected_vars_categorical), all_of(selected_vars_numerical)) %>%
    mutate(size_select = ifelse(size_select <= 0.000000000001, 0.000000000001, size_select)) %>%
    drop_na() %>%
    group_by(newest_node_name) %>%
    sample_n(1) %>%
    rename("taxon" = newest_node_name) %>%
    mutate_at(all_of(selected_vars_categorical), factor) %>%
    mutate_at(c(selected_vars_numerical), as.numeric) %>%
    drop_na() %>%
    bind_rows(basal_node_type1)
  inf_mw_type1 <- infer_edgelist(data1, col_taxon = "taxon", col_num_size = selected_vars_numerical, return_full_matrix = FALSE, print_dropped_taxa = FALSE, cat_combo_list = trait_combos)

  # Only need to predict for metazoans
  data2 <- meta %>%
    filter(web == i & kingdom == "Animalia" & drop_summary == 0) %>%
    dplyr::select(newest_node_name, all_of(selected_vars_categorical), all_of(selected_vars_numerical)) %>%
    mutate(size_select = ifelse(size_select <= 0.000000000001, 0.000000000001, size_select)) %>%
   drop_na() %>%
    group_by(newest_node_name) %>%
    sample_n(1) %>%
    rename("taxon" = newest_node_name) %>%
    mutate_at(all_of(selected_vars_categorical), factor) %>%
    mutate_at(c(selected_vars_numerical), as.numeric) %>%
    drop_na() %>%
    bind_rows(basal_node_type2)
  inf_mw_type2 <- infer_edgelist(data2, col_taxon = "taxon", col_num_size = selected_vars_numerical, return_full_matrix = FALSE, print_dropped_taxa = FALSE, cat_combo_list = trait_combos)


  web_list_inferred_mw_basaltype1[[i]] <- inf_mw_type1
  web_list_inferred_mw_basaltype2[[i]] <- inf_mw_type2

  print(i)
}


compare_inferred_list_basal <- list(
  inferred_mw_basaltype1 = web_list_inferred_mw_basaltype1,
  inferred_mw_basaltype2 = web_list_inferred_mw_basaltype2
)

stats_compare_inferred_basaltype <- c()
for (i in names(compare_inferred_list_basal)) {
  select_list <- compare_inferred_list_basal[[i]]

  for (j in names(select_list)) {
    select_list2 <- select_list[[j]]

    select_fw <- as.matrix(select_list2)
    select_fw <- unique(select_fw[, ])
    ccc <- cbind(stack(calc_network_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = 0)
    stats_compare_inferred_basaltype <- bind_rows(stats_compare_inferred_basaltype, ccc)
  }

  print(i)
}


stats_compare_inferred2_basaltype <- stats_compare_inferred_basaltype %>%
  dplyr::rename("variable" = "ind", "value" = "values") %>%
  filter(!variable %in% c("soi_sw", "mean_oi_std", "mean_oi_sw", "mean_tl_sw", "mot_lin", "mot_omn", "mot_ap_comp", "mot_dir_comp")) %>%
  mutate(web = factor(web, levels = websizeorder))
stats_compare_inferred2_basaltype$variable <- jlookup(stats_compare_inferred2_basaltype$variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")
stats_compare_inferred2_basaltype <- stats_compare_inferred2_basaltype %>% filter(!is.na(variable))

ggplot(stats_compare_inferred2_basaltype) +
  geom_point(aes(x = web, y = value, color = type, shape = type)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free")









#### COMPARE LINK DISTRIBUTIONS: Metawebs v species webs----



compare_inferred_list <- list(
  inferred_mw = web_list_inferred_mw,
  inferred_species = web_list_inferred_species,
  fc2 = mstr_web_list[["fc2"]]
)

stats_compare_inferred <- c()
for (i in names(compare_inferred_list)) {
  select_list <- compare_inferred_list[[i]]

  for (j in names(select_list)) {
    select_list2 <- select_list[[j]]

    if (is.matrix(select_list2)) {
      select_fw <- as.matrix(select_list2)
      select_fw <- unique(select_fw[, ])
      ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = 0)
      stats_compare_inferred <- bind_rows(stats_compare_inferred, ccc)
    } else {
      for (k in 1:length(select_list2)) {
        select_list3 <- select_list2[[k]]

        select_fw <- as.matrix(select_list3)
        select_fw <- unique(select_fw[, ])
        ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = k)
        stats_compare_inferred <- bind_rows(stats_compare_inferred, ccc)
      }
    }

    print(j)
  }

  print(i)
}


stats_compare_inferred2 <- stats_compare_inferred %>%
  dplyr::rename("variable" = "ind", "value" = "values") %>%
  mutate(web = factor(web, levels = websizeorder))

stats_compare_inferred2$variable <- jlookup(stats_compare_inferred2$variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")
stats_compare_inferred2 <- stats_compare_inferred2 %>% filter(!is.na(variable))

ggplot(stats_compare_inferred2 %>%
  mutate(type = recode(type, fc2 = "Empirical web", inferred_mw = "Inferred metaweb", inferred_species = "Hypothetical realized web")) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new))) +
  geom_point(aes(x = web, y = value, color = type, shape = type)) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free", nrow = 3)



stats_minmax <- stats_compare_inferred2 %>%
  filter(type %in% c("inferred_species")) %>%
  group_by(web, variable, type) %>%
  summarize(
    avg = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    mid = quantile(value, 0.5, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )
stats_emp <- stats_compare_inferred2 %>% filter(type == "fc2")

ggplot(stats_emp) +
  geom_point(data = stats_emp, aes(x = web, y = value, color = type)) +
  geom_errorbar(data = stats_minmax, aes(x = web, ymin = lower, ymax = upper, color = type)) +
  facet_wrap(. ~ variable, scales = "free")


ggplot(stats_compare_inferred2 %>%
  mutate(type = recode(type, fc2 = "Empirical web", inferred_mw = "Inferred metaweb", inferred_species = "Hypothetical realized web")) %>%
  filter(type != "Hypothetical realized web") %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new))) +
  geom_line(aes(x = web, y = value, color = type, group = type)) +
  geom_line(data = stats_minmax %>% mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
    mutate(web = factor(web, levels = fig_web_names$new)) %>%
    mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
    filter(variable != "NA") %>%
    mutate(variable = factor(variable, levels = fig_names$new)), aes(x = factor(web), y = avg, color = type, group = 1)) +
  facet_wrap(. ~ variable, scales = "free")


##FIGURE6##
dat1 <- stats_compare_inferred2 %>%
  mutate(type = recode(type, fc2 = "Empirical web", inferred_mw = "Inferred metaweb", inferred_species = "Hypothetical realized web")) %>%
  filter(type != "Hypothetical realized web") %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(type = factor(type, levels = c("Empirical web", "Inferred metaweb", "Hypothetical realized web"))) %>%
  mutate(variable = factor(variable, levels = fig_names$new)) %>%
  filter(!variable %in% c("Size", "Mean deg."))
dat2 <- stats_minmax %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(type = recode(type, inferred_species = "Hypothetical realized web")) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(type = factor(type, levels = c("Empirical web", "Inferred metaweb", "Hypothetical realized web"))) %>%
  mutate(variable = factor(variable, levels = fig_names$new)) %>%
  filter(!variable %in% c("Size", "Mean deg."))
temp <- dat1 %>% bind_rows(dat2)

ggplot(temp) +
  geom_errorbar(aes(x = factor(web), ymin = min, ymax = max, color = type, group = 1), alpha = 0.5, width = 0.5, show.legend = F) +
  geom_point(aes(x = factor(web), y = avg, color = type, group = 1), show.legend = T) +
  geom_point(aes(x = web, y = value, color = type, group = type), size = 2, show.legend = T) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(. ~ variable, scales = "free", nrow = 3)





stats_joined <- stats_minmax %>%
  left_join(stats_emp, by = c("web", "variable")) %>%
  mutate(overlap = ifelse(value >= min & value <= max, "OVERLAP", "NO"))


## Model error
mes <- calc_ME(value = stats_joined$value, model_lower = stats_joined$min, model_median = stats_joined$mid, model_upper = stats_joined$max)
stats_joined_me <- cbind(stats_joined, model_error = mes)

ggplot(stats_joined_me) +
  geom_point(aes(x = web, y = model_error, color = overlap), size = 2) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = -1) +
  facet_grid(. ~ variable, scales = "free")

ggplot(stats_joined_me) +
  geom_point(aes(x = web, y = model_error, color = variable), size = 2) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = -1)

library(RColorBrewer)
ggplot(stats_joined_me %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new))) +
  geom_point(aes(x = variable, y = model_error, color = web, shape = web), size = 2) +
  geom_hline(yintercept = 1) +
  ylab("Model error") +
  xlab("Metric") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = -1)


temp <- stats_joined_me %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new)) %>%
  group_by(web) %>%
  summarize(avg_web = mean(abs(model_error), na.rm = TRUE))

ggplot(temp) +
  geom_point(aes(x = web, y = avg_web), size = 2) +
  geom_hline(yintercept = 1) +
  ylab("Model error") +
  xlab("Metric") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = -1)

#### COMPARE LINK DISTRIBUTIONS: TSS of replicate webs----

tss_realized_web_stats <- c()
for (i in webs_to_run) {
  web_fc2 <- as.data.frame(mstr_web_list[["fc2"]][[i]]) %>% distinct()

  web_fc2[, 1] <- gsub("Non-metazoan", "BASAL NODE", web_fc2[, 1])
  web_fc2[, 2] <- gsub("Non-metazoan", "BASAL NODE", web_fc2[, 2])
  colnames(web_fc2) <- c("res_new_node_name", "con_new_node_name")

  tax_names <- c(unique(as.character(as.matrix(web_fc2))))
  tax_list <- as.data.frame(expand.grid(tax_names, tax_names))
  colnames(tax_list) <- c("res_new_node_name", "con_new_node_name")

  realized_web_list <- web_list_inferred_species[[i]]

  for (j in 1:length(realized_web_list)) {
    web_realized <- as.data.frame(realized_web_list[[j]]) %>% distinct()
    colnames(web_realized) <- c("res_new_node_name", "con_new_node_name")

    # First create list of all possible comparisons
    pred_comp <- as.data.frame(tax_list) %>%
      left_join((as.data.frame(web_fc2) %>% mutate(empirical_int = 1)), by = c("res_new_node_name", "con_new_node_name")) %>%
      full_join(as.data.frame(web_realized %>% mutate(inferred_int = 1)), by = c("res_new_node_name", "con_new_node_name"))
    pred_comp[is.na(pred_comp)] <- 0

    cm <- Calc_TSS_new(as.data.frame(cbind(pred_comp[, c("empirical_int")], pred_comp[, c("inferred_int")])))

    tss_realized_web_stats <- rbind(tss_realized_web_stats, cbind(cm, web = i, replicate = j, type = "realized"))

    print(j)
  }

  mw_web_list <- web_list_inferred_mw[[i]]

  web_mw <- as.data.frame(mw_web_list)
  colnames(web_mw) <- c("res_new_node_name", "con_new_node_name")

  # First create list of all possible comparisons
  pred_comp <- tax_list %>%
    left_join(as.data.frame(web_fc2) %>% mutate(empirical_int = 1), by = c("res_new_node_name", "con_new_node_name")) %>%
    full_join(as.data.frame(web_mw %>% mutate(inferred_int = 1)), by = c("res_new_node_name", "con_new_node_name"))
  pred_comp[is.na(pred_comp)] <- 0

  cm <- Calc_TSS_new(as.data.frame(cbind(pred_comp[, c("empirical_int")], pred_comp[, c("inferred_int")])))

  tss_realized_web_stats <- rbind(tss_realized_web_stats, cbind(cm, web = i, replicate = 0, type = "metaweb"))

  print(j)

  print(i)
}

tss_realized_web_stats2 <- as.data.frame(tss_realized_web_stats) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new))

ggplot(tss_realized_web_stats2) +
  geom_density(data = tss_realized_web_stats2 %>% filter(type == "realized"), aes(x = TSS, color = web)) +
  ylab("Density") +
  geom_vline(data = tss_realized_web_stats2 %>% filter(type == "metaweb"), aes(xintercept = TSS, color = web), size = 2, linetype = "dashed") +
  theme(legend.title = element_blank()) +
  scale_color_brewer(palette = "Set1")


#### TABLE FOR PAPER
t1 <- tss_realized_web_stats2 %>%
  pivot_longer(cols = c(a, b, c, d)) %>%
  mutate(name = gsub("a", "True positive", name)) %>%
  mutate(name = gsub("b", "False positive", name)) %>%
  mutate(name = gsub("c", "False negative", name)) %>%
  mutate(name = gsub("d", "True negative", name)) %>%
  group_by(web, type, replicate) %>%
  mutate(total = sum(value, na.rm = TRUE))
t1 <- t1 %>%
  mutate(prop = value / total) %>%
  group_by(web, type, name, total) %>%
  summarize(mean = mean(prop, na.rm = TRUE), sd = sd(prop, na.rm = TRUE)) %>%
  pivot_wider(names_from = c("name"), values_from = c("mean", "sd"))
# TSS values
t2 <- tss_realized_web_stats2 %>%
  group_by(web, type) %>%
  summarise(mean_TSS = mean(TSS, na.rm = TRUE)) %>%
  right_join(t1, by = c("web", "type"))

temp <- tss_realized_web_stats2 %>%
  pivot_longer(cols = c(a, b, c, d)) %>%
  mutate(name = gsub("a", "True positive", name)) %>%
  mutate(name = gsub("b", "False positive", name)) %>%
  mutate(name = gsub("c", "False negative", name)) %>%
  mutate(name = gsub("d", "True negative", name))


ggplot(temp) +
  geom_point(aes(x = name, y = value)) +
  facet_wrap(web ~ type, scales = "free", ncol = 2)

#### NEW FIGURE DEC 2
temp1 <- temp %>% filter(type == "metaweb")

temp2 <- temp %>%
  filter(type == "realized") %>%
  group_by(web, name) %>%
  summarize(
    avg = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.25, na.rm = TRUE),
    mid = quantile(value, 0.5, na.rm = TRUE),
    upper = quantile(value, 0.75, na.rm = TRUE)
  )
p1 <- ggplot(temp1) +
  geom_point(data = temp1, aes(x = name, y = value), color = "#619EFF") +
  geom_errorbar(data = temp2, aes(x = name, ymin = lower, ymax = upper), color = "#43BA38") +
  xlab("Link type") +
  ylab("Number of links") +
  facet_grid(web ~ ., scales = "free")

temp3 <- temp %>%
  select(web, TSS, type) %>%
  distinct()
p2 <- ggplot(temp3) +
  geom_density(data = temp3 %>% filter(type == "realized"), aes(x = TSS), color = "#43BA38") +
  geom_vline(data = temp3 %>% filter(type == "metaweb"), aes(xintercept = TSS), color = "#619EFF") +
  ylab("Density") +
  facet_grid(web ~ ., scales = "free")

p1 + p2 + plot_layout(widths = c(2, 1))



#### Alt metric--overall accurracy
##FIGURE5##
temp <- as.data.frame(tss_realized_web_stats) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  dplyr::select(-TSS) %>%
  mutate(Accuracy = ((a + d) / (a + b + c + d))) %>%
  mutate(Sensitivity = (a / (a + c))) %>%
  mutate(Specificity = (d / (b + d))) %>%
  mutate(TSS = (Sensitivity + Specificity - 1)) %>%
  pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, TSS)) %>%
  filter(name != "Accuracy")
ggplot(temp) +
  geom_density(data = temp %>% filter(type == "realized"), aes(x = value), color = "#619EFF") +
  geom_vline(data = temp %>% filter(type == "metaweb"), aes(xintercept = value), color = "#43BA38") +
  xlab("Value") +
  ylab("Density") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  facet_grid(name ~ web, scales = "free_y")

#### COMPARE LINK DISTRIBUTIONS: Just TSS of inferred metaweb----

tss_inferred_web_stats <- c()
tss_inferred_web_taxa <- c()
for (i in webs_to_run) {
  web_fc2 <- mstr_web_list[["fc2"]][[i]]

  web_fc2[, 1] <- gsub("Non-metazoan", "BASAL NODE", web_fc2[, 1])
  web_fc2[, 2] <- gsub("Non-metazoan", "BASAL NODE", web_fc2[, 2])

  colnames(web_fc2) <- c("res_new_node_name", "con_new_node_name")

  tax_names <- c(unique(as.character(web_fc2)))
  tax_list <- as.data.frame(expand.grid(tax_names, tax_names))
  colnames(tax_list) <- c("res_new_node_name", "con_new_node_name")

  mw_web_list <- web_list_inferred_mw[[i]]

  web_mw <- as.data.frame(mw_web_list)
  colnames(web_mw) <- c("res_new_node_name", "con_new_node_name")

  # First create list of all possible comparisons
  pred_comp <- tax_list %>%
    left_join(as.data.frame(web_fc2) %>% mutate(empirical_int = 1), by = c("res_new_node_name", "con_new_node_name")) %>%
    full_join((web_mw %>% mutate(inferred_int = 1)), by = c("res_new_node_name", "con_new_node_name")) %>%
    distinct()
  pred_comp[is.na(pred_comp)] <- 0

  tss_inferred_web_taxa <- rbind(tss_inferred_web_taxa, cbind(pred_comp, web = i))

  cm <- Calc_TSS_new(as.data.frame(cbind(pred_comp[, c("empirical_int")], pred_comp[, c("inferred_int")])))

  tss_inferred_web_stats <- rbind(tss_inferred_web_stats, cbind(cm, web = i, replicate = 0, type = "metaweb"))

  print(i)
}

tss_inferred_web_stats2 <- as.data.frame(tss_inferred_web_stats) %>%
  dplyr::select(-c(type, replicate)) %>%
  pivot_longer(cols = -c(web, TSS)) %>%
  group_by(web, TSS) %>%
  mutate(total_links = sum(value, na.rm = TRUE)) %>%
  mutate(prop_links = value / total_links) %>%
  mutate(name = gsub("a", "True positive", name)) %>%
  mutate(name = gsub("b", "False positive", name)) %>%
  mutate(name = gsub("c", "False negative", name)) %>%
  mutate(name = gsub("d", "True negative", name)) %>%
  mutate(TSS = paste("TSS =", round(TSS, digits = 2))) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(name = factor(name, levels = c("True negative", "True positive", "False negative", "False positive")))



ggplot(data = tss_inferred_web_stats2, aes(x = web, fill = name, y = prop_links)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Proportion of all possible links") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  ) +
  facet_grid(. ~ web + TSS, scales = "free")



#### COMPARE LINK DISTRIBUTIONS: Hypothetical realized node stats----
hyp_node_stats <- c()
for (i in names(web_list_inferred_species)) {
  for (j in 1:length(web_list_inferred_species[[i]])) {
    sub_graph <- web_list_inferred_species[[i]][[j]]
    temp <- calc_node_stats(graph = graph_from_edgelist(sub_graph, directed = TRUE)) %>% mutate(web = i, replicate = j)

    hyp_node_stats <- rbind(hyp_node_stats, temp)
  }
  print(i)
}

hyp_node_stats2 <- hyp_node_stats %>%
  dplyr::select(taxon, norm_btw, norm_degree_in, norm_degree_out, tl_std, oi_std, web, replicate) %>%
  pivot_longer(cols = c(norm_btw, norm_degree_in, norm_degree_out, tl_std, oi_std), names_to = c("name"), values_to = c("value"))

hyp_node_stats3 <- hyp_node_stats2 %>%
  group_by(taxon, web, name) %>%
  summarise(
    avg = median(value, na.rm = TRUE),
    range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE)
  ) %>%
  left_join(var_names, by = c("name" = "variable_short_name")) %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
  mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
  filter(variable_long_name != "NA") %>%
  mutate(valuetype = "Hypothetical realized web") %>%
  mutate(value = avg)
ggplot(hyp_node_stats3) +
  geom_point(aes(x = avg, y = range)) +
  facet_wrap(web ~ name, scales = "free")

# Compare with empirical and inferred from script 7

all_node_comparison <- compare_inferred_nodes2 %>%
  mutate(web = jlookup(web, fig_web_names, matching_col = "old", new_values = "new")) %>%
  mutate(variable_long_name = jlookup(variable_long_name, fig_names, matching_col = "old", new_values = "new")) %>%
  mutate(variable_long_name = factor(variable_long_name, levels = fig_names$new)) %>%
  filter(variable_long_name != "NA") %>%
  pivot_longer(c("Inferred metaweb_0", "Empirical web_0"), names_to = "valuetype", values_to = "value") %>%
  bind_rows(hyp_node_stats3) %>%
  mutate(web = factor(web, levels = fig_web_names$new)) %>%
  mutate(valuetype = recode(valuetype, `Empirical web_0` = "Empirical web", `Inferred metaweb_0` = "Inferred metaweb")) %>%
  mutate(valuetype = factor(valuetype, levels = c("Empirical web", "Inferred metaweb", "Hypothetical realized web")))

##SUPFIGURE5##
ggplot(all_node_comparison) +
  geom_boxplot(aes(x = valuetype, y = value, color = valuetype), show.legend = FALSE) +
  facet_grid(variable_long_name ~ web, scales = "free") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())




#### FUNDAMENTAL NETWORKS: Setting up PFIM constants----

load("Data/trait_rules.rda")

# List of all used traits
traits_motility <- as.data.frame(example_data_traitrules) %>%
  filter(trait_type_resource == "motility" | trait_type_consumer == "motility") %>%
  dplyr::select(trait_resource, trait_consumer)
traits_motility <- unique(unlist(traits_motility))

traits_tiering <- as.data.frame(example_data_traitrules) %>%
  filter(trait_type_resource == "tiering" | trait_type_consumer == "tiering") %>%
  dplyr::select(trait_resource, trait_consumer)
traits_tiering <- unique(unlist(traits_tiering))

traits_feeding <- as.data.frame(example_data_traitrules) %>%
  filter(trait_type_resource == "feeding" | trait_type_consumer == "feeding") %>%
  dplyr::select(trait_resource, trait_consumer)
traits_feeding <- unique(unlist(traits_feeding))

# All possible trait combinations
traits_all <- as.data.frame(expand.grid(traits_motility, traits_tiering, traits_feeding))

# Specifically assign trait combos to taxa in niche web
taxon_traits <- cbind(traits_all)
colnames(taxon_traits) <- c("motility", "tiering", "feeding")

# Drop Basal taxa to make a faunal list like a fossil list (i.e., PFIM input)
# Set up basal node
basal_node <- data.frame(
  stringsAsFactors = FALSE,
  taxon = c("BASAL NODE", "BASAL NODE", "BASAL NODE"),
  tiering = c("t_pela", "t_surf", "t_shin"),
  motility = c("m_nmun", "m_nmun", "m_nmun"),
  feeding = c("f_prim", "f_prim", "f_prim") # ,size = c(1e-12, 1e-12, 1e-12)
)

taxon_traits <- taxon_traits %>%
  mutate(taxon = paste("taxon", as.character(1:nrow(.)), sep = "")) %>%
  filter(feeding != "f_prim") %>%
  bind_rows(basal_node)


# PFIM inferrence
predicted_web <- infer_edgelist(taxon_traits, col_taxon = "taxon")
# 42 taxa dropped = of all 168 possible trait combinations, only 126 are valid combos (excluding size)

predicted_graph <- graph_from_edgelist(predicted_web)


# Calc trophic level of different taxa
pg_tl <- as.data.frame(stack(calc_node_tl_std(igraph::as_adjacency_matrix(predicted_graph))))
colnames(pg_tl) <- c("trophic_level", "taxon")

trait_list_tl <- pg_tl %>%
  left_join(taxon_traits, by = "taxon")


nodes_meta <- trait_list_tl %>%
  filter(taxon != "BASAL NODE")
nodes_basal <- trait_list_tl %>%
  filter(taxon == "BASAL NODE")

nodes_extra <- nodes_meta %>%
  mutate(taxon = gsub(x = taxon, pattern = "taxon", replacement = "extra"))

# Std web
fw_std <- nodes_meta %>%
  bind_rows(nodes_extra) %>%
  bind_rows(nodes_basal)

meta_size <- nrow(nodes_meta)

# Top heavy web
nodes_meta_upper <- nodes_meta %>%
  slice_max(trophic_level, n = (meta_size / 2), with_ties = FALSE) %>%
  mutate(taxon = gsub(x = taxon, pattern = "taxon", replacement = "extra_high"))
nodes_meta_lower <- nodes_meta %>%
  slice_min(trophic_level, n = (meta_size / 2), with_ties = FALSE) %>%
  mutate(taxon = gsub(x = taxon, pattern = "taxon", replacement = "extra_low"))

fw_top <- (nodes_meta_upper %>% mutate(taxon = paste(taxon, "_1", sep = ""))) %>%
  bind_rows(nodes_meta_upper %>% mutate(taxon = paste(taxon, "_2", sep = ""))) %>%
  bind_rows(nodes_meta_upper %>% mutate(taxon = paste(taxon, "_3", sep = ""))) %>%
  bind_rows(nodes_meta_lower) %>%
  bind_rows(nodes_basal)

fw_bottom <- (nodes_meta_lower %>% mutate(taxon = paste(taxon, "_1", sep = ""))) %>%
  bind_rows(nodes_meta_lower %>% mutate(taxon = paste(taxon, "_2", sep = ""))) %>%
  bind_rows(nodes_meta_lower %>% mutate(taxon = paste(taxon, "_3", sep = ""))) %>%
  bind_rows(nodes_meta_upper) %>%
  bind_rows(nodes_basal)


#### FUNDAMENTAL NETWORKS: Compare stats of unsampled food webs----

fw_fl_std <- fw_std %>% dplyr::select(-trophic_level)
fw_fl_top <- fw_top %>% dplyr::select(-trophic_level)
fw_fl_bottom <- fw_bottom %>% dplyr::select(-trophic_level)

runs <- 5

meta_std <- infer_edgelist(fw_fl_std, col_taxon = "taxon")
meta_std <- as.matrix(as.data.frame(meta_std) %>% filter(taxon_resource != taxon_consumer))
realized_std <- powerlaw_prey(meta_std, n_samp = runs)

meta_top <- infer_edgelist(fw_fl_top, col_taxon = "taxon")
meta_top <- as.matrix(as.data.frame(meta_top) %>% filter(taxon_resource != taxon_consumer))
realized_top <- powerlaw_prey(meta_top, n_samp = runs)

meta_bottom <- infer_edgelist(fw_fl_bottom, col_taxon = "taxon")
meta_bottom <- as.matrix(as.data.frame(meta_bottom) %>% filter(taxon_resource != taxon_consumer))
realized_bottom <- powerlaw_prey(meta_bottom, n_samp = runs)

# Make list of webs
compare_list <- list(
  std = list(metaweb = meta_std, realized = realized_std),
  bottom = list(metaweb = meta_bottom, realized = realized_bottom),
  top = list(metaweb = meta_top, realized = realized_top)
)

stats_compare_inferred_int <- c()
stats_node_inferred <- c()
for (i in names(compare_list)) {
  select_list <- compare_list[[i]]

  for (j in names(select_list)) {
    select_list2 <- select_list[[j]]

    if (is.matrix(select_list2)) {
      select_fw <- as.matrix(select_list2)
      select_fw <- unique(select_fw[, ])
      ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = 0)
      stats_compare_inferred_int <- bind_rows(stats_compare_inferred_int, ccc)

      bbb <- cbind(calc_node_stats(graph_from_edgelist(select_fw, directed = T)), web = j, type = i, replicate = 0)
      stats_node_inferred <- bind_rows(stats_node_inferred, bbb)
    } else {
      for (k in 1:length(select_list2)) {
        select_list3 <- select_list2[[k]]

        select_fw <- as.matrix(select_list3)
        select_fw <- unique(select_fw[, ])
        ccc <- cbind(stack(calc_select_stats(graph_from_edgelist(select_fw, directed = T))), web = j, type = i, replicate = k)
        stats_compare_inferred_int <- bind_rows(stats_compare_inferred_int, ccc)

        bbb <- cbind(calc_node_stats(graph_from_edgelist(select_fw, directed = T)), web = j, type = i, replicate = k)
        stats_node_inferred <- bind_rows(stats_node_inferred, bbb)
      }
    }
  }

  print(i)
}


stats_compare_inferred2_int <- stats_compare_inferred_int %>%
  dplyr::rename("variable" = "ind", "value" = "values") %>%
  filter(!variable %in% c("size", "mean_normalized_degree"))


#### FIGURE
ggplot(stats_compare_inferred2_int) +
  geom_point(aes(x = type, y = value, color = web, size = web), alpha = 0.5) +
  scale_size_manual(values = c(3, 1)) +
  xlab("Web shape") +
  ylab("Metric value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free")


temp <- stats_compare_inferred2_int %>%
  mutate(type = recode(type, bottom = "Bottom-heavy", std = "Standard", top = "Top-heavy")) %>%
  mutate(web = recode(web, realized = "Realized webs", metaweb = "Metaweb")) %>%
  mutate(variable = jlookup(variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new))

temp_points <- temp %>% filter(web == "Metaweb")
temp_bars <- temp %>%
  filter(web != "Metaweb") %>%
  group_by(variable, web, type) %>%
  dplyr::summarize(
    avg = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    mid = quantile(value, 0.5, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )

ggplot(temp_points) +
  geom_errorbar(data = temp_bars, aes(x = type, ymin = lower, ymax = upper), color = "#609DFF") +
  geom_point(data = temp_points, aes(x = type, y = value), color = "#00BA39", alpha = 0.5) +
  xlab("Web shape") +
  ylab("Metric value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free")


ggplot(temp_points) +
  geom_errorbar(data = temp_bars, aes(x = type, ymin = lower, ymax = upper), color = "#609DFF") +
  geom_point(data = temp_points, aes(x = type, y = value), color = "#00BA39", alpha = 0.5) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(variable ~ ., scales = "free")


temp_std <- stats_compare_inferred2_int %>% filter(web == "metaweb")
temp_other <- stats_compare_inferred2_int %>%
  filter(web != "metaweb") %>%
  group_by(variable, web, type) %>%
  dplyr::summarize(
    avg = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.25, na.rm = TRUE),
    mid = quantile(value, 0.5, na.rm = TRUE),
    upper = quantile(value, 0.75, na.rm = TRUE)
  )

temp_join <- temp_std %>% full_join(temp_other, by = c("variable", "type"))

me_vals <- c()
for (i in 1:nrow(temp_join)) {
  temp_me <- calc_ME(
    value = temp_join[i, c("value")],
    model_lower = temp_join[i, c("lower")],
    model_median = temp_join[i, c("mid")],
    model_upper = temp_join[i, c("upper")]
  ) %>% as.numeric()
  me_vals <- c(me_vals, temp_me)
}

temp_join2 <- as.data.frame(cbind(temp_join, model_error = me_vals))
temp_join2$model_error[is.nan(temp_join2$model_error)] <- 0

ggplot(temp_join2 %>%
  mutate(type = recode(type, bottom = "Bottom-heavy", std = "Standard", top = "Top-heavy")) %>%
  # mutate(web=recode(web, realized="Realized webs",metaweb="Metaweb")) %>%
  mutate(variable = jlookup(variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new))) +
  geom_hline(yintercept = 1, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "darkgrey", linetype = "dashed") +
  geom_point(aes(x = type, y = model_error), color = "#00BA39", alpha = 0.8) +
  xlab("") +
  ylab("Effect size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_wrap(variable ~ ., scales = "free")



#### FUNDAMENTAL NETWORKS: Analyze node level stats----

# Separated plot of metaweb and realized
stats_nodelevel <- as.data.frame(stats_node_inferred) %>%
  pivot_longer(cols = !c("taxon", "web", "type", "replicate"), names_to = "variable", values_to = "value") %>%
  mutate(value = as.numeric(as.character(value))) %>%
  filter(variable %in% c("norm_btw", "norm_degree_in", "norm_degree_out", "oi_std", "tl_std"))

stats_nodelevel_realized <- stats_nodelevel %>%
  filter(web == "realized") %>%
  group_by(variable, taxon, type, web) %>%
  dplyr::summarize(
    avg = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )

ggplot(stats_nodelevel_realized, aes(x = type, y = sd)) +
  geom_boxplot() +
  xlab("Metaweb replicate") +
  ylab("SD of taxon-specific metric values across realized webs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free", space = "free_x")

ggplot(stats_nodelevel_realized, aes(x = type, y = avg)) +
  geom_boxplot() +
  xlab("Metaweb replicate") +
  ylab("SD of taxon-specific metric values across realized webs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free", space = "free_x")

#### FIGURE
ggplot(stats_nodelevel_realized %>%
  mutate(type = recode(type, bottom = "Bottom-heavy", std = "Standard", top = "Top-heavy")) %>%
  # mutate(web=recode(web, realized="Realized webs",metaweb="Metaweb")) %>%
  mutate(variable = jlookup(variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new)), aes(x = type, y = sd)) +
  geom_boxplot(color = "#609DFF") +
  xlab("Web shape") +
  ylab("SD of taxon-specific metric values across realized webs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free", space = "free_x")



# Avg and SD plot
temp <- stats_nodelevel_realized %>%
  mutate(type = recode(type, bottom = "Bottom-heavy", std = "Standard", top = "Top-heavy")) %>%
  mutate(variable = jlookup(variable, var_names, matching_col = "variable_short_name", new_values = "variable_long_name")) %>%
  mutate(variable = jlookup(variable, fig_names, matching_col = "old", new_values = "new")) %>%
  filter(variable != "NA") %>%
  mutate(variable = factor(variable, levels = fig_names$new))
temp2 <- temp %>% pivot_longer(cols = c("avg", "sd"))
ggplot(temp2, aes(x = type, y = value)) +
  geom_boxplot(color = "#609DFF") +
  xlab("Web shape") +
  ylab("SD of taxon-specific metric values across realized webs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ name, scales = "free")

p1 <- ggplot(temp, aes(x = type, y = avg)) +
  stat_boxplot(geom = "errorbar", color = "#609DFF", width = 0.5) +
  geom_boxplot(color = "#609DFF", fill = "white", outlier.alpha = 0.5, outlier.shape = 21) +
  xlab("") +
  ylab("Average taxon-specific metric values across realized webs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free")
p2 <- ggplot(temp, aes(x = type, y = sd)) +
  stat_boxplot(geom = "errorbar", color = "#609DFF", width = 0.5) +
  geom_boxplot(color = "#609DFF", fill = "white", outlier.alpha = 0.5, outlier.shape = 21) +
  xlab("") +
  ylab("SD of taxon-specific metric values across realized webs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  facet_grid(variable ~ ., scales = "free")

p1 + p2


#### GLOBI/ADBM COMPARISON----

# Set up ABDM model


NM.RH.web <- readRDS("Data/ythan_adbm.RDS")

abdm_fw <- NM.RH.web$web

sp <- read.csv("Data/cirtwill_ythan_taxa.csv")

bses <- sp %>%
  arrange(BodyWeight) %>%
  pull(Species)

g2 <- abdm_fw
rownames(g2) <- bses
colnames(g2) <- bses
g2 <- g2[bses, bses]


plot_data2 <- as.data.frame(as_edgelist(graph_from_adjacency_matrix(as.matrix(g2)))) %>%
  rename(taxon_consumer = V2) %>%
  rename(taxon_resource = V1) %>%
  mutate(
    taxon_consumer = factor(taxon_consumer, levels = (bses)),
    taxon_resource = factor(taxon_resource, levels = rev(bses))
  ) %>%
  mutate(interaction_abdm = 1)





# Set up globi
inferred_interactions_df <- read.csv("Data/inferred_interactions_list.csv") %>% filter(web == "ythanjacob")

temp <- inferred_interactions_df %>% # This is from script 8
  filter(web == "ythanjacob") %>%
  mutate(taxon_resource2 = word(taxon_resource, 1)) %>%
  mutate(taxon_consumer2 = word(taxon_consumer, 1))

ythan_taxa <- unique(c(temp$taxon_resource, temp$taxon_consumer))

interaction_list <- read.csv("Data/ythan_globi.csv")

interaction_list_clean <- interaction_list %>%
  dplyr::select(source, target) %>%
  mutate(interaction_globi = 1) %>%
  distinct()

temp3 <- temp %>%
  # filter(interaction_empirical==0 & interaction_inferred==1) %>% #look at interactions inferred but not observed
  left_join(interaction_list_clean, by = c("taxon_resource2" = "source", "taxon_consumer2" = "target")) %>%
  filter(web == "ythanjacob") %>%
  replace(is.na(.), 0)

temp4 <- temp3 %>%
  left_join(plot_data2 %>% filter(taxon_resource %in% ythan_taxa) %>%
    filter(taxon_consumer %in% ythan_taxa), by = c("taxon_consumer", "taxon_resource")) %>%
  replace(is.na(.), 0)

# Make long
bses <- meta %>%
  filter(web == "ythanjacob" & drop_summary == 0 & kingdom == "Animalia") %>%
  arrange(size_select) %>%
  pull(newest_node_name)
bses <- unique(bses)

temp5 <- temp4 %>%
  select(taxon_resource, taxon_consumer, interaction_empirical, interaction_inferred, interaction_globi, interaction_abdm) %>%
  mutate_at(c("interaction_empirical", "interaction_inferred", "interaction_globi", "interaction_abdm"), as.numeric) %>%
  pivot_longer(c(interaction_empirical, interaction_inferred, interaction_globi, interaction_abdm)) %>%
  filter(value == 1) %>%
  mutate(
    taxon_consumer = factor(taxon_consumer, levels = (bses)),
    taxon_resource = factor(taxon_resource, levels = rev(bses))
  )

temp6 <- temp5 %>%
  filter(taxon_resource %in% bses) %>%
  filter(taxon_consumer %in% bses) %>%
  mutate(
    taxon_consumer = factor(taxon_consumer, levels = (bses)),
    taxon_resource = factor(taxon_resource, levels = rev(bses))
  )

length(unique(temp6$taxon_resource))
length(unique(temp6$taxon_consumer))

ggplot(temp6 %>% filter(name == "interaction_empirical"), aes(x = taxon_consumer, y = taxon_resource, fill = name)) +
  geom_tile(alpha = 0.5)

# Make blank matrix in order to make the axes equal
blank <- expand.grid(bses, bses) %>%
  rename(taxon_consumer = Var2) %>%
  rename(taxon_resource = Var1) %>%
  mutate(
    taxon_consumer = factor(taxon_consumer, levels = (bses)),
    taxon_resource = factor(taxon_resource, levels = rev(bses))
  ) %>%
  mutate(interaction_all = 1)

ggplot(temp6) +
  geom_tile(data = blank, aes(x = taxon_consumer, y = taxon_resource), fill = "white") +
  geom_tile(data = temp6, aes(x = taxon_consumer, y = taxon_resource, fill = name), alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0),
    aspect.ratio = 1
  ) +
  facet_grid(. ~ name)

# Plot with all other webs compared to inferred
ggplot(temp6) +
  geom_tile(data = blank, aes(x = taxon_consumer, y = taxon_resource), fill = "white") +
  geom_tile(data = temp6 %>% filter(name != "interaction_inferred") %>% mutate(name = "all"), aes(x = taxon_consumer, y = taxon_resource, fill = name), alpha = 0.5) +
  geom_tile(data = temp6 %>% filter(name == "interaction_inferred"), aes(x = taxon_consumer, y = taxon_resource, fill = name), alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0),
    aspect.ratio = 1
  )

temp6$name <- factor(temp6$name, levels = c("interaction_empirical", "interaction_inferred", "interaction_globi", "interaction_abdm"))

# Create the adjacency matrix plot
ggplot() +
  geom_tile(data = blank, aes(x = taxon_consumer, y = taxon_resource), fill = "white") +
  geom_tile(data = temp6 %>%
    mutate(name = gsub("interaction_empirical", "Empirical", name)) %>%
    mutate(name = gsub("interaction_inferred", "Inferred", name)) %>%
    mutate(name = gsub("interaction_globi", "GloBI", name)) %>%
    mutate(name = gsub("interaction_abdm", "ABDM", name)), aes(x = taxon_consumer, y = taxon_resource)) +
  theme_bw() +
  scale_x_discrete(drop = FALSE, position = "top") +
  scale_y_discrete(drop = FALSE, position = "left") +
  geom_tile(data = as.data.frame(bses), aes(x = bses, y = bses), alpha = 0.3) +
  theme(
    # Rotate the x-axis labels
    axis.text.x = element_text(angle = 270, hjust = 0),
    # Force  into a square aspect ratio
    aspect.ratio = 1,
    # Hide the legend
    legend.position = "none"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  ) +
  facet_grid(. ~ name)



temp7 <- temp4 %>%
  select(taxon_resource, taxon_consumer, interaction_empirical, interaction_inferred, interaction_globi, interaction_abdm) %>%
  mutate_at(c("interaction_empirical", "interaction_inferred", "interaction_globi", "interaction_abdm"), as.numeric) %>%
  filter(taxon_resource %in% bses) %>%
  filter(taxon_consumer %in% bses) %>%
  mutate(
    taxon_consumer = factor(taxon_consumer, levels = (bses)),
    taxon_resource = factor(taxon_resource, levels = rev(bses))
  ) %>%
  mutate(empirical_inferred = paste(interaction_empirical, interaction_inferred, sep = ",")) %>%
  mutate(empirical_inferred = str_replace(empirical_inferred, "0,0", "True negative")) %>%
  mutate(empirical_inferred = str_replace(empirical_inferred, "0,1", "False positive")) %>%
  mutate(empirical_inferred = str_replace(empirical_inferred, "1,0", "False negative")) %>%
  mutate(empirical_inferred = str_replace(empirical_inferred, "1,1", "True positive")) %>%
  mutate(empirical_inferred = factor(empirical_inferred, levels = c("True positive", "False positive", "False negative", "True negative"))) %>%
  filter(empirical_inferred == "False positive") %>%
  mutate(interaction_any = ifelse((interaction_inferred == 1 & interaction_globi == 0 & interaction_abdm == 0), "PFWIM",
    ifelse((interaction_inferred == 1 & interaction_globi == 1 & interaction_abdm == 0), "PFWIM + GloBI",
      ifelse(interaction_inferred == 1 & interaction_abdm == 1 & interaction_globi == 0, "PFWIM + ABDM",
        ifelse((interaction_globi == 1 & interaction_abdm == 1 & interaction_inferred == 1), "PFWIM + GloBI + ABDM", "NONE")
      )
    )
  ))


ttt <- temp7 %>%
  group_by(empirical_inferred, interaction_any) %>%
  summarise(n = n()) %>%
  group_by(empirical_inferred) %>%
  mutate(freq = n / sum(n))

ggplot(ttt) +
  geom_col(aes(x = empirical_inferred, y = n, fill = interaction_any))
ggplot(ttt) +
  geom_col(aes(x = empirical_inferred, y = freq, fill = interaction_any))

#### FIGURE
ggplot(ttt) +
  geom_col(aes(x = empirical_inferred, y = n, fill = interaction_any)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("") +
  ylab("Number of false positive interactions") +
  theme(
    legend.title = element_blank()
  )
