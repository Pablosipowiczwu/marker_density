rm(list = ls(all.names = TRUE))
library(ASRgenomics)
library(AGHmatrix)
require(asreml)
require(MASS)
require(dplyr)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
######Load marker data from Sequence Capture and from DArTag 3K and make a vector with families genotyped by both platforms to later extract information from each  harvest
load("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/Xr_matrix_MAF0.05.RData")
Data<-read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/DaRT/AlleleFreqParents.csv")
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
new<-merge(data, newmale, by = "Family")
new$Juliana.Code2<-as.factor(new$Juliana.Code2)
info2 <- read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/info2.csv")
info2<-info2[,-9]
info2$Family<-as.factor(info2$Family)
newinfo<-info2[info2$Juliana.Code2 %in% rownames(NewXr),]

#######Randomly select 500 markers (10 times) to estimate G matrix and create G matrix with AGHmatrix package
A<-list()
B<-list()
C<-list()
D<-list()
#E<-list()
#G<-list()
Total<-list()
newData<-list()
Check<-list()
Gmatrix<-list()
for (i in 1:10) {
  set.seed(1051+i)  
  Tot<-c(1:2190)
  A[[i]] <- sample(Tot,10, replace = FALSE)
  n <- Tot[c(Tot %not in% A[[i]])]
  B[[i]] <- sample(n,40, replace = FALSE)
  n <- Tot[c(Tot %not in% c(A[[i]],B[[i]]))]
  C[[i]] <- sample(n,50, replace = FALSE)
  n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]]))]
  D[[i]] <- sample(n,400, replace = FALSE)
  n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]],D[[i]]))]
#  E[[i]] <- sample(n,500, replace = FALSE)
#  n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]],D[[i]],E[[i]]))]
#  G[[i]] <- sample(n,500, replace = FALSE)
#  n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]],D[[i]],E[[i]],G[[i]]))]
  Total[[i]]<-c(A[[i]]
                ,B[[i]]
                ,C[[i]]
                ,D[[i]]
                #,E[[i]]
                #,G[[i]]
  )
  Check[[i]]<-length(unique(Total[[i]]))###This line is to check if there are actually 500 markers in each iteration
  newData[[i]]<-Data[Data$Marker %in% Total[[i]],]
  newData[[i]]<-newData[[i]][,-c(1:2)]
  newData[[i]]<-t(newData[[i]])
  rownames(newData[[i]])<-gsub("X","F",rownames(newData[[i]]))
  newData[[i]]<-newData[[i]][rownames(newData[[i]]) %in% Family,]
  Gmatrix[[i]]<-Gmatrix(newData[[i]], method="VanRaden", ratio = TRUE, ploidy = 4)
}

##Estimate the inverse of each Gmatrix to run models

GINV<-list()

for (i in 1:10) {
  GINV[[i]] <- G.inverse(G=Gmatrix[[i]], sparseform=TRUE,blend=T)$Ginv.sparse 
}

####Cross-validations
####The k loop will retrieve phenotypic data from each harvest from 1 to 11
####The i loop will generate 100 iterations of cross-validation
####The loop j will use a different Gmatrix built from 500 markers (10 different Gmatrices by randomly selecting 50 markers at)
for (k in 1:11) {
  hn <- k
  data_byharv <- subset(new,Harvest==hn)
  eblueH1<-read.csv(paste0("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/BLUES/DMY/BLUE_H",hn,'.csv'), header=TRUE,row.names = 1)
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
metrics<-matrix(0,nrow = 1000, ncol = 7)
colnames(metrics)<-c('harvest','GINV','pcor','pcor_spearman','mse','slope','trn_fam')
require(asreml)
for (i in 1:100) {
  set.seed(1)
  val[[i]]<-data.frame(newinfo %>% group_by(groups_dapc) %>% 
   dplyr::slice_sample(prop = 0.117, replace = F))
  val[[i]]<-merge(val[[i]],neweblueH1b,by="Juliana.Code2")
  data2$Juliana.Code2<-data2$Juliana.Code2.x
  data10fold<-data2
  data10fold$DM_kg_ha[data10fold$Juliana.Code2 %in% val[[1]]$Juliana.Code2]<-NA 
  genot_validation<-unique(val[[i]]$Juliana.Code2)
  for (j in 1:10) {
    modelharvest1<-asreml(fixed = DM_kg_ha ~ 1,
                          random = ~ Row + Column +     vm(Juliana.Code2,GINV[[j]]), 
                          na.action = na.method(y='include',x='include'),
                          data=data10fold)
    modelharvest1<-update.asreml(modelharvest1)
    modelharvest1<-update.asreml(modelharvest1)
    
    predicted<-data.frame(predict(modelharvest1,classify = "Juliana.Code2")$pvals)
    eblues_validation<-neweblueH1b[neweblueH1b$Juliana.Code2 %in% genot_validation,]
    predicted_validation<-data.frame(predicted[predicted$Juliana.Code2 %in% genot_validation,])
    all<-merge(predicted_validation,eblues_validation, by.x = "Juliana.Code2", 
               by.y = "Juliana.Code2")
    validation[[j]]<-all
    pcor<-cor(all$eblue,all$predicted.value)
    mse<-mean((all$eblue-all$predicted.value)^2)
    slope<-lm(all$eblue~all$predicted.value)$coef[2]
    pcor_spearman<-cor(all$eblue,all$predicted.value,method = 'spearman')
    GINV_j<-j
    trn_fam_i<-i
    metrics[(j*100)+i-100,]<-c(paste0("H ",hn),GINV_j,pcor,pcor_spearman,mse,slope,trn_fam_i)
    print(paste('All done',j))
    metrics<-as.data.frame(metrics)
    metrics$pcor<-as.numeric(metrics$pcor)
    metrics$pcor_spearman<-as.numeric(metrics$pcor_spearman)
    metrics$mse<-as.numeric(metrics$mse)
    metrics$slope<-as.numeric(metrics$slope)
    mean(metrics[,3],na.rm=T)
    median(metrics[,3],na.rm=T)
    sd(metrics[,3],na.rm=T)
    
    setwd("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/DaRT/AlleleFreq/GBLUP/New/500markers/Harvests")
    write.csv(metrics,paste0('./metrics10fold_DaRT_500markers_DMY_GBLUP',hn, '.csv'),quote=FALSE)
  }
}
}
