#load packages
library(mediation)
library(lme4)
library(lmerTest)
library(reshape2)
library(dplyr)
library(tidyr)
library(parallel)

#合并表格(plasma_meta)
species_table <- as.data.frame(t(species_otu_table))
genus_table <- as.data.frame(t(genus_otu_table))
meta.tab <- merge(species_table,genus_table,by='row.names')
meta.tab <- merge(metadata[c(1:8,23:37)],meta.tab,by.x='SampleID',by.y="Row.names",all=T)
meta.tab <- merge(meta.tab,plasma_metabolites_meta[,c(1:2,9:219)],by=c("Subject","Time"),all=T)
meta.tab[,9:ncol(meta.tab)] <- lapply(meta.tab[,9:ncol(meta.tab)],scale)

species <- c(sig_species,sig_genus)
marker <- colnames(metadata)[c(23:27,29:37)]
metabolite <- plasma_metabolite_sig
covar <- c("Age","Gender")

params_list <- expand.grid(species, marker, metabolite)
colnames(params_list) <- c("spe","ser","m")
#Group1
g1.mid.tab <- subset(meta.tab,Treatment=="Group1")

# forward mediation --------------------
g1_Mediation.forward_plasma <- NULL
#创建函数
run_forward_mediation_analysis <- function(spe, ser, m){
  spe<-as.character(spe)
  ser<-as.character(ser)
  m<-as.character(m)
  
  tab1 <- g1.mid.tab[,c("Age","Gender","Subject",spe,ser,m)]
  tab1 <- na.omit(tab1)
  id = paste(spe, m, ser, sep = "|")
  
  model.m = lme4::lmer(as.formula( paste("`",m,"`", " ~ ", spe,"+ (1 | Subject)",sep = "")), data = tab1)
  model.y = lme4::lmer(as.formula( paste(ser, " ~ ", spe," + ", "`",m,"`","+ (1 | Subject)",sep = "")), data = tab1)
  summary = summary(mediate(model.m,model.y,treat=spe,mediator=m,boot=F,sims=500))
  res <- capture.output(summary,append=FALSE)
  
  tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" & tmp!="."]
  tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
  ACME.p <- tmp[length(tmp)]
  ACME.Es <- tmp[3]
  tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" & tmp!="."]
  tmp <- tmp[!grepl("*",tmp,fixed = T) ]
  ADE.p <- tmp[length(tmp)]
  ADE.Es <- tmp[3]
  tmp <- base::strsplit(res[grep("Total",res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" & tmp!="."]
  tmp <- tmp[!grepl("*",tmp,fixed = T) ]
  Total.p <- tmp[length(tmp)]
  Total.Es <- tmp[3]
  tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" ]
  tmp <- tmp[!grepl("*",tmp,fixed = T)]
  i_str = which(grepl("Mediated", tmp))
  prop.mediated <- tmp[(i_str + 2)]
  prop.p <- tmp[length(tmp)]
  return(c(id, ACME.p, ACME.Es, ADE.p, ADE.Es, Total.p, Total.Es, prop.mediated, prop.p))
}

#并行循环
results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_forward_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96) # 使用所有可用核心 - 1

#合并结果
g1_Mediation.forward_plasma <- do.call(rbind, results)


# reverse mediation --------------------
g1_Mediation.reverse_plasma <- NULL

#创建函数
run_reverse_mediation_analysis <- function(spe, ser, m){
  spe<-as.character(spe)
  ser<-as.character(ser)
  m<-as.character(m)
  
  tab1 <- g1.mid.tab[,c("Age","Gender","Subject",spe,ser,m)]
  tab1 <- na.omit(tab1)
  id = paste(spe, ser, m, sep = "|")
  
  model.m = lme4::lmer(as.formula( paste(ser, " ~ ", spe,"+ (1 | Subject)", sep = "")), data = tab1)
  model.y = lme4::lmer(as.formula( paste("`",m,"`", " ~ ", spe," + ", ser,"+ (1 | Subject)",sep = "")), data = tab1)
  summary = summary(mediate(model.m,model.y,treat=spe,mediator=ser,boot=F,sims=500))
  res <- capture.output(summary,append=FALSE)
  
  tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" & tmp!="."]
  tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
  ACME.p <- tmp[length(tmp)]
  ACME.Es <- tmp[3]
  
  tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" & tmp!="."]
  tmp <- tmp[!grepl("*",tmp,fixed = T) ]
  ADE.p <- tmp[length(tmp)]
  ADE.Es <- tmp[3]
  
  tmp <- base::strsplit(res[grep("Total",res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" & tmp!="."]
  tmp <- tmp[!grepl("*",tmp,fixed = T) ]
  Total.p <- tmp[length(tmp)]
  Total.Es <- tmp[3]
  
  tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
  tmp <- tmp[tmp != "" ]
  tmp <- tmp[!grepl("*",tmp,fixed = T)]
  i_str = which(grepl("Mediated", tmp))
  prop.mediated <- tmp[(i_str + 2)]
  prop.p <- tmp[length(tmp)]
  
  return(c(id, ACME.p, ACME.Es, ADE.p, ADE.Es, Total.p, Total.Es, prop.mediated, prop.p))
}

#并行循环
results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_reverse_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96) # 使用96核心

#合并结果
g1_Mediation.reverse_plasma <- do.call(rbind, results)


###group2
g1.mid.tab <- subset(meta.tab,Treatment=="Group2")

# forward mediation --------------------
g2_Mediation.forward_plasma <- NULL

results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_forward_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96)

g2_Mediation.forward_plasma <- do.call(rbind, results)


# reverse mediation --------------------
g2_Mediation.reverse_plasma <- NULL

results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_reverse_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96)

g2_Mediation.reverse_plasma <- do.call(rbind, results)



####合并表格(fecal_meta)
species <- c(sig_species,sig_genus)
metabolite <- fecal_metabolites_sig
marker <- colnames(metadata)[c(23:27,29:37)]

species_table <- as.data.frame(t(species_otu_table))
genus_table <- as.data.frame(t(genus_otu_table))
meta.tab <- merge(species_table,genus_table,by='row.names')
meta.tab <- merge(metadata[c(1:8,23:37)],meta.tab,by.x='SampleID',by.y="Row.names",all=T)
meta.tab <- merge(meta.tab,fecal_metabolites_meta[,c(1:2,9:219)],by=c("Subject","Time"),all=T)
meta.tab[,9:ncol(meta.tab)] <- lapply(meta.tab[,9:ncol(meta.tab)],scale)

params_list <- expand.grid(species, marker, metabolite)
colnames(params_list) <- c("spe","ser","m")

#Group1
g1.mid.tab <- subset(meta.tab,Treatment=="Group1")

# forward mediation --------------------
g1_Mediation.forward_fecal <- NULL

results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_forward_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96)

g1_Mediation.forward_fecal <- do.call(rbind, results)

# reverse mediation --------------------
g1_Mediation.reverse_fecal <- NULL

results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_reverse_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96)

g1_Mediation.reverse_fecal <- do.call(rbind, results)


#Group2
g1.mid.tab <- subset(meta.tab,Treatment=="Group2")

# forward mediation --------------------
g2_Mediation.forward_fecal <- NULL

results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_forward_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96)

g2_Mediation.forward_fecal <- do.call(rbind, results)

# reverse mediation --------------------
g2_Mediation.reverse_fecal <- NULL

results <- mclapply(1:nrow(params_list), function(i) {
  cat(i,"\n")
  run_reverse_mediation_analysis(params_list$spe[i], params_list$ser[i], params_list$m[i])
}, mc.cores = 96)

g2_Mediation.reverse_fecal <- do.call(rbind, results)


