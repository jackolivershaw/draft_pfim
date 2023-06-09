library(igraph)
library(tidyverse)


####Functions----

##Make perfect trophic species web from adjacency matrix
make.trophic.species.web <- function(fw) {
  if(length(fw)>1){
    
    dupl_cols <- duplicated(fw)
    dupl_rows <- duplicated(t(fw))
    dupl <- which(dupl_rows+dupl_cols == 2)
    
    newnw <- as.data.frame(fw[-dupl,-dupl])
    
  }else{
    newnw<-NULL
  }
  return(newnw)
}

##Calculate short weighted trophic level - via Susanne Kortsch
calc.swtl<-function(web){
  
  if(length(web)>1){
    
    library(NetIndices)
    library(igraph)  
    
    web<-as.matrix(web)
    rownames(web)<-colnames(web)
    
    basal <- rownames(subset(web, apply(web, 2, sum)==0) & apply(web, 1, sum)!= 0)
    edge.list_web <- graph.adjacency(web, mode = "directed");
    paths_prey <- shortest.paths(graph = edge.list_web, v= V(edge.list_web),
                                 to = V(edge.list_web)[basal], mode = "in", weights = NULL, algorithm = "unweighted")
    paths_prey[is.infinite(paths_prey)] <- NA
    shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm=TRUE)))
    #longest_paths <- as.matrix(apply(paths_prey, 1, max, na.rm=TRUE))
    in_deg <- apply(web, 2, sum) 				## Take each columns, sum  the rows
    out_deg <- apply(web, 1, sum) 			## Take each rows, sum the columns
    
    
    # Shortest TL
    sTL <- 1 + shortest_paths  # Commonly, detritus have a TL value of 1. (Shortest path to basal = 0)
    # Longest TL
    #lTL <- 1 + longest_paths
    
    S<-dim(web)[1]
    
    # Creating the matrix  
    short_TL_matrix <- matrix(NA, S, S)
    #long_TL_matrix <- matrix(NA, S, S)
    prey_ave <- rep(NA, S) # averaged prey
    chain_ave<-rep(NA, S)
    #
    for(j in 1:S){
      for(i in 1:S){
        lij <- web[i, j] 					# is the interaction matrix
        prey_ave[j] <- 1 / in_deg[j] 			# if Bas species in, no in-degrees ; re-attribute a value of 0
        short_TL_matrix[i,j] <- lij * sTL[i]   	# Shortest path in matrix
        #long_TL_matrix[i,j] <- lij * lTL[i]
      }  
    }
    
    prey_ave[which(prey_ave == Inf)] <- 0
    prey_ave[which(prey_ave == -Inf)] <- 0
    
    short_TL_matrix[is.nan(short_TL_matrix)] <- 0
    
    
    short_TL_matrix[which(short_TL_matrix == Inf)] <- 0
    short_TL_matrix[which(short_TL_matrix == -Inf)] <- 0
    
    
    
    sumShortTL <- as.matrix(apply(short_TL_matrix, 2, sum)) # sum all shortest path
    
    
    # Short-weighted TL
    # weigth by the number of prey
    
    SWTL <- matrix(data = NA, nrow = S, ncol = 1)
    for(i in 1:S){
      SWTL[i] <- 1 + (prey_ave[i] * sumShortTL[i])  
    }
    
    
    TL_NI<-TrophInd(web)
    TL_NetInd<- TL_NI[,1]
    
    TLS<-cbind(SWTL, TL_NetInd)
    colnames(TLS)<-c("SWTL","NetInd_TL")
    rownames(TLS)<-rownames(web)
    
    return(TLS)
  }else{
    return(NA)
  }
  
}



####Burgess from supplemental
sup_dat<-read.csv("/Users/jackshaw/Downloads/testing_dunne_08/raw_dunne_08_files_from_supplement/t7_burgess_fw_data_tbl1.csv")
temp<-sup_dat[,c("res_number","con_number")]
temp<-na.omit(as.matrix(temp))
temp[] <- sapply(temp, function(x) paste('ID', x,sep=""))
temp_el<-as.data.frame(temp)
t2<-graph_from_edgelist(temp,directed=TRUE)
vcount(t2)
# food web connectance
(ecount(t2))/(vcount(t2))^2
t3<-as.matrix(as_adjacency_matrix(t2))
t4<-as.data.frame(calc.swtl(t3))
mean(t4$SWTL,na.rm=TRUE)
max(t4$SWTL,na.rm=TRUE)

#ID-name list
sup_cons<-sup_dat[,c("con_number","consumer_sp")]
sup_cons<-unique(sup_cons[,])
colnames(sup_cons)<-c("id_num","id_name")
sup_cons$id_num<-paste("ID",sup_cons$id_num,sep="")
sup_res<-sup_dat[,c("res_number","resource_sp")]
sup_res<-unique(sup_res[,])
colnames(sup_res)<-c("id_num","id_name")
sup_res$id_num<-paste("ID",sup_res$id_num,sep="")
sup_join<-rbind(sup_cons,sup_res)
sup_join<-unique(sup_join[,])

sup_join$id_name<-tolower(gsub("[[:punct:][:blank:]]+", " ", sup_join$id_name))
sup_join<-as.data.frame(na.omit(as.matrix(sup_join)))


#Import new names from meta (post processing in "compiled_scripts.R")
temp_meta<-meta %>% filter(web=="burgess") %>%
  dplyr::select(oldest_node_name,newest_node_name)
temp_meta$oldest_node_name<-tolower(gsub("[[:punct:][:blank:]]+", " ", temp_meta$oldest_node_name))
temp_meta$newest_node_name<-tolower(gsub("[[:punct:][:blank:]]+", " ", temp_meta$newest_node_name))
temp_meta<- temp_meta %>%
  left_join(sup_join,by=c("oldest_node_name"="id_name")) %>%
  distinct() %>%
  dplyr::select(-oldest_node_name)

#Make new edgelist
final_el<-temp_el %>%  left_join(temp_meta,by=c("res_number"="id_num"))
names(final_el)[names(final_el) == 'newest_node_name'] <- 'res_name'
final_el<-final_el %>%  left_join(temp_meta,by=c("con_number"="id_num"))
names(final_el)[names(final_el) == 'newest_node_name'] <- 'con_name'
burgess_el <- final_el[,c("res_name","con_name")] %>% distinct()



####Chengjiang from supplemental
sup_dat<-read.csv("/Users/jackshaw/Downloads/testing_dunne_08/raw_dunne_08_files_from_supplement/t6_chengjiang_fw_data_tbl1.csv")
temp<-sup_dat[,c("res_number","con_number")]
temp<-na.omit(as.matrix(temp))
temp[] <- sapply(temp, function(x) paste('ID', x,sep=""))
temp_el<-as.data.frame(temp)
t2<-graph_from_edgelist(temp,directed=TRUE)
vcount(t2)
# food web connectance
(ecount(t2))/(vcount(t2))^2
t3<-as.matrix(as_adjacency_matrix(t2))
t4<-as.data.frame(calc.swtl(t3))
mean(t4$SWTL,na.rm=TRUE)
max(t4$SWTL,na.rm=TRUE)

#ID-name list
sup_cons<-sup_dat[,c("con_number","consumer_sp")]
sup_cons<-unique(sup_cons[,])
colnames(sup_cons)<-c("id_num","id_name")
sup_cons$id_num<-paste("ID",sup_cons$id_num,sep="")
sup_res<-sup_dat[,c("res_number","resource_sp")]
sup_res<-unique(sup_res[,])
colnames(sup_res)<-c("id_num","id_name")
sup_res$id_num<-paste("ID",sup_res$id_num,sep="")
sup_join<-rbind(sup_cons,sup_res)
sup_join<-unique(sup_join[,])

sup_join$id_name<-tolower(gsub("[[:punct:][:blank:]]+", " ", sup_join$id_name))
sup_join<-as.data.frame(na.omit(as.matrix(sup_join)))

#Import new names from meta (post processing in "compiled_scripts.R")
temp_meta<-meta %>% filter(web=="chengjiang") %>%
  dplyr::select(oldest_node_name,newest_node_name)
temp_meta$oldest_node_name<-tolower(gsub("[[:punct:][:blank:]]+", " ", temp_meta$oldest_node_name))
temp_meta$newest_node_name<-tolower(gsub("[[:punct:][:blank:]]+", " ", temp_meta$newest_node_name))
temp_meta<- temp_meta %>%
  left_join(sup_join,by=c("oldest_node_name"="id_name")) %>%
  distinct() %>%
  dplyr::select(-oldest_node_name)

#Make new edgelist
final_el<-temp_el %>%  left_join(temp_meta,by=c("res_number"="id_num"))
names(final_el)[names(final_el) == 'newest_node_name'] <- 'res_name'
final_el<-final_el %>%  left_join(temp_meta,by=c("con_number"="id_num"))
names(final_el)[names(final_el) == 'newest_node_name'] <- 'con_name'
chengjiang_el <- final_el[,c("res_name","con_name")] %>% distinct()

####Export new files----
og_che<-chengjiang_el
og_bur<-burgess_el

write.csv(og_che, "chengjiang_dunne08_exported_221128.csv")
write.csv(og_bur, "burgess_dunne08_exported_221128.csv")
