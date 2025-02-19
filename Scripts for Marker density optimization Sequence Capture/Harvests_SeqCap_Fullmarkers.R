
rm(list = ls(all.names = TRUE))
library(ASRgenomics)
library(AGHmatrix)
require(asreml)
require(MASS)
require(dplyr)

setwd("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/Sequence Capture/Full")
load("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Xr_matrix_MAF0.05.RData")
######Load marker data from Sequence Capture and from DArTag 3K and keep only genotypic data from samples genotyped by both platforms
Parents<-paste0("F", c(1072:1111))
DART<-read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/Haplotypes DArT/GMatrixHaplo.csv")
Samples<-paste0("F", DART$X)
NewXr<-Xr[rownames(Xr) %in% Samples,]
Family<-rownames(NewXr)
###Load phenotypic data from trial and curate data to match between genotypic files and phenotypic files
data=read.csv('/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/z.csv',header=TRUE)
data$Family<-gsub("_","",as.factor(data$Family))
data<-transform(data, Harvest=factor(Harvest), Row=factor(Row), Column=factor(Column),Family=factor(Family), FamC=factor(FamC), Type=factor(Type)) 
male.info <- read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/New_AlfalfaPed_info_ERiosCurated2.csv", h=T)
colnames(male.info)[1]<-"Sample"
colnames(male.info)[3]<-"Family"
male.info$Family<-gsub("_","",as.factor(male.info$Family))
male.info$Family<-as.factor(male.info$Family)
male.info$Juliana.Code2<-as.factor(male.info$Juliana.Code2)
newmale<-male.info[male.info$Juliana.Code2 %in% Family,]
newdata<-merge(data, newmale, by = "Family")
newdata$Juliana.Code2<-as.factor(newdata$Juliana.Code2)
info2 <- read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/info2.csv")
info2<-info2[,-9]
info2$Family<-as.factor(info2$Family)
newinfo<-info2[info2$Juliana.Code2 %in% rownames(NewXr),]

Gmatrix<-Gmatrix(NewXr, method="VanRaden", ratio = TRUE, ploidy = 4, ploidy.correction = TRUE)
GINV <- G.inverse(G=Gmatrix, sparseform=TRUE,blend=T)$Ginv.sparse  

####Cross-validations
####The k loop will retrieve phenotypic data from each harvest from 1 to 11
####The i loop will generate 100 iterations of cross-validation
for (k in 1:11) {
  hn <- k
data_byharv <- subset(newdata,Harvest==hn)
eblueH1<-read.csv(paste0("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/BLUES/BLUE_H",hn,'.csv'), header=TRUE,row.names = 1)
eblueH1$Family<-gsub("_","",as.factor(eblueH1$Family))
eblueH1$Family<-as.factor(eblueH1$Family)
eblueH1b<-merge(eblueH1, male.info, by='Family', all.x = TRUE)
neweblueH1b<-eblueH1b[eblueH1b$Juliana.Code2 %in% Family,]

Maleinfo<-male.info$Unique.Code
male.info$Unique.Code<-paste0("F",Maleinfo)
data2<-merge(male.info, data_byharv, by='Family')

colnames(eblueH1b)[2:3]<-c('eblue', 'std.error_eblue')
colnames(neweblueH1b)[2:3]<-c('eblue', 'std.error_eblue')

val<-list()
validation<-list()
predicted<-list()
metrics<-matrix(0,nrow = 1000, ncol = 6)
colnames(metrics)<-c('harvest','pcor','pcor_spearman','mse','slope','trn_fam')

for (i in 1:100) {
  set.seed(i)
  val[[i]]<-data.frame(newinfo %>% group_by(groups_dapc) %>% dplyr::slice_sample(prop = 0.117, replace = F))
  val[[i]]<-merge(val[[i]],neweblueH1b,by="Juliana.Code2")
  data10fold<-newdata
  data10fold$DM_kg_ha[data10fold$Juliana.Code2 %in% val[[i]]$Juliana.Code2]<-NA 
  genot_validation<-unique(val[[i]]$Juliana.Code2)
    modelharvest1<-asreml(fixed = DM_kg_ha ~ 1,
                          random = ~ Row + Column + vm(Juliana.Code2,GINV), 
                          na.action = na.method(y='include',x='include'),
                          data=data10fold)
    modelharvest1<-update.asreml(modelharvest1)
    modelharvest1<-update.asreml(modelharvest1)
    
    predicted<-data.frame(predict(modelharvest1,classify = "Juliana.Code2")$pvals)
    eblues_validation<-neweblueH1b[neweblueH1b$Juliana.Code2 %in% genot_validation,]
    predicted_validation<-data.frame(predicted[predicted$Juliana.Code2 %in% genot_validation,])
    all<-merge(predicted_validation,eblues_validation, by.x = "Juliana.Code2", 
               by.y = "Juliana.Code2")
    validation[[i]]<-all
    pcor<-cor(all$eblue,all$predicted.value)
    mse<-mean((all$eblue-all$predicted.value)^2)
    slope<-lm(all$eblue~all$predicted.value)$coef[2]
    pcor_spearman<-cor(all$eblue,all$predicted.value,method = 'spearman')
    trn_fam_i<-i
    metrics[i,]<-c(paste0("H ",hn),pcor,pcor_spearman,mse,slope,trn_fam_i)
    metrics<-as.data.frame(metrics)
    metrics$pcor<-as.numeric(metrics$pcor)
    metrics$pcor_spearman<-as.numeric(metrics$pcor_spearman)
    metrics$mse<-as.numeric(metrics$mse)
    metrics$slope<-as.numeric(metrics$slope)
    mean(metrics[,3],na.rm=T)
    median(metrics[,3],na.rm=T)
    sd(metrics[,3],na.rm=T)
    setwd("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/Sequence Capture/Fullmarkers/Harvests")
    write.csv(metrics,paste0('./metrics10fold_SeqCapFullmarkers_DMY_',hn, '.csv'),quote=FALSE)
  }
}

