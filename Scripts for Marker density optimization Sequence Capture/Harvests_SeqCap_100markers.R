
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

#Create Gmatrix from random sampling of 100 markers ten times

Columns<-list()
Data<-list()
Gmatrix<-list()
Check<-list()
A<-list()
B<-list()
C<-list()
#D<-list()
#E<-list()
#FF<-list()
#G<-list()
#H<-list()
#I<-list()
#J<-list()
Total<-list()
Check<-list()
for (i in 1:10) {
  set.seed(1051+i)
  A[[i]] <- assign(paste0("Col", i), runif(10, min = 1, max = 114945))
  B[[i]] <- assign(paste0("Col", i), runif(40, min = 1, max = 114945))
  C[[i]] <- assign(paste0("Col", i), runif(50, min = 1, max = 114945))
  #D[[i]] <- assign(paste0("Col", i), runif(400, min = 1, max = 114945))
  #E[[i]] <- assign(paste0("Col", i), runif(500, min = 1, max = 114945))
  #FF[[i]] <- assign(paste0("Col", i), runif(4000, min = 1, max = 114945))
  #G[[i]] <- assign(paste0("Col", i), runif(5000, min = 1, max = 114945))
  #H[[i]] <- assign(paste0("Col", i), runif(10000, min = 1, max = 114945))
  #I[[i]] <- assign(paste0("Col", i), runif(20000, min = 1, max = 114945))
  #J[[i]] <- assign(paste0("Col", i), runif(30000, min = 1, max = 114945))
  Total[[i]]<-c(A[[i]]
                ,B[[i]]
                ,C[[i]]
                #,D[[i]]
                #,E[[i]]
                #,FF[[i]]
                #,G[[i]]
                #,H[[i]]
                #,I[[i]]
                #,J[[i]]
  )
  Check[[i]]<-length(unique(Total[[i]]))
  Columns[[i]]<- NewXr[,c(Total[[i]])]
  Data[[i]]<-as.matrix(Columns[[i]], nrow = 160, ncol = length(Total[[i]]), dimnames(NewXr))
  Gmatrix[[i]]<-Gmatrix(Data[[i]], method="VanRaden", ratio = TRUE, ploidy = 4, ploidy.correction = TRUE)
}

GINV<-list()
for (i in 1:10) {
  GINV[[i]] <- G.inverse(G=Gmatrix[[i]], sparseform=TRUE,blend=T)$Ginv.sparse  
}

####Cross-validations
####The k loop will retrieve phenotypic data from each harvest from 1 to 11
####The i loop will generate 100 iterations of cross-validation
####The j loop will use a different Gmatrix built from 100 markers (10 different Gmatrices by randomly selecting 50 markers at)
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
metrics<-matrix(0,nrow = 1000, ncol = 7)
colnames(metrics)<-c('harvest','GINV','pcor','pcor_spearman','mse','slope','trn_fam')
require(asreml)

setwd("~/Documents/Forage Breeding Team/AGP17/Genomic Selection/Sequence Capture/Validation")

for (i in 1:100) {
  set.seed(i)
  val[[i]]<-data.frame(newinfo %>% group_by(groups_dapc) %>% dplyr::slice_sample(prop = 0.117, replace = F))
  val[[i]]<-merge(val[[i]],neweblueH1b,by="Juliana.Code2")
  data10fold<-newdata
  data10fold$DM_kg_ha[data10fold$Juliana.Code2 %in% val[[i]]$Juliana.Code2]<-NA 
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
    
    setwd("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/Sequence Capture/100markers/Harvests")
    write.csv(metrics,paste0('./metrics10fold_SeqCap100markers_DMY_',hn, '.csv'),quote=FALSE)
  }
}
}
