rm(list = ls(all.names = TRUE))
library(ASRgenomics)
library(AGHmatrix)
require(asreml)
require(MASS)
require(dplyr)
######Load marker data from Sequence Capture and from DArTag 3K and keep only genotypic data from samples genotyped by both platforms
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
load("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/Xr_matrix_MAF0.05.RData")
Data<-read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/DaRT/AlleleFreqParents.csv")
DomData<-read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Genomic Selection/DaRT/DomAlleleFreqParents.csv")
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
#######Randomly select 100 markers (10 times) to estimate G matrix and create G matrix with AGHmatrix package
A<-list()
B<-list()
C<-list()
#D<-list()
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
  #D[[i]] <- sample(n,400, replace = FALSE)
  #n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]],D[[i]]))]
  #E[[i]] <- sample(n,500, replace = FALSE)
  #n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]],D[[i]],E[[i]]))]
  #  G[[i]] <- sample(n,500, replace = FALSE)
   # n <- Tot[c(Tot %not in% c(A[[i]],B[[i]],C[[i]],D[[i]],E[[i]],G[[i]]))]
  Total[[i]]<-c(A[[i]]
                ,B[[i]]
                ,C[[i]]
             #   ,D[[i]]
              #  ,E[[i]]
               # ,G[[i]]
  )
  Check[[i]]<-length(unique(Total[[i]]))
  newData[[i]]<-Data[Data$Marker %in% Total[[i]],]
  newData[[i]]<-newData[[i]][,-c(1:2)]
  newData[[i]]<-t(newData[[i]])
  rownames(newData[[i]])<-gsub("X","F",rownames(newData[[i]]))
  newData[[i]]<-newData[[i]][rownames(newData[[i]]) %in% Family,]
  Gmatrix[[i]]<-Gmatrix(newData[[i]], method="VanRaden", ratio = TRUE, ploidy = 4)
}

GINV<-list()
for (i in 1:10) {
  GINV[[i]] <- G.inverse(G=Gmatrix[[i]], sparseform=TRUE,blend=T)$Ginv.sparse 
}

#### Cross-validations (to test prediction of totally unrelated testing and training populations 
#### we will predict a full-sib family and any family including the male or female from that specific 
#### full-sib family will be eliminated from the training population)


####The k loop will retrieve phenotypic data from each harvest from 1 to 11
####The i loop will select a specific male 
####The j loop will select a specific female 
####The l loop will use a different Gmatrix built from 100 markers (10 different Gmatrices by randomly selecting 50 markers)

males<-unique(male.info$Male.name)
males<-males[1:6]
data2<-data2[data2$Male.name.x %in% males,]
females<-unique(male.info$Female.name)
females<-females[1:27]


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
  all_n<-data.frame()
  r1 = data.frame("Harvest","Male","Female","Ginv","pcor","mse")
  r1 = r1[-1,]
  for (l in 1:10) {
    all_n<-data.frame()
  for (i in 1:6) {
  for (j in 1:27) {
  m<-males[i]
  fe<-females[j]
  index_males<-which(data2$Male.name.x == m)
  index_females<-which(data2$Female.name.x == fe)
  data10fold<-data2
  data10fold$DM_kg_ha[c(index_males,index_females)]<-NA ###This works!!
  
  val <- data10fold[data10fold$Male.name.x == m & data10fold$Female.name.x == fe, ]
  genot_validation<-unique(val$Juliana.Code2.x)
  
  colnames(eblueH1b)[2:3]<-c('eblue', 'std.error_eblue')
  colnames(neweblueH1b)[2:3]<-c('eblue', 'std.error_eblue')
  

  validation<-list()
  predicted<-list()

  colnames(r1)<-c("Harvest","Male","Female","Ginv","pcor","mse")
    
      modelharvest1<-asreml(fixed = DM_kg_ha ~ 1,
                            random = ~ Row + Column + vm(Juliana.Code2.x,GINV[[l]])
                            #      + vm(Juliana.Code2,DINV[[j]])
                            , 
                            na.action = na.method(y='include',x='include'),
                            data=data10fold)
      modelharvest1<-update.asreml(modelharvest1)
      modelharvest1<-update.asreml(modelharvest1)
      
      predicted[[l]]<-data.frame(predict(modelharvest1,classify = "Juliana.Code2.x")$pvals)
      eblues_validation<-neweblueH1b[neweblueH1b$Juliana.Code2 %in% genot_validation,]
      predicted_validation<-data.frame(predicted[[l]][predicted[[l]]$Juliana.Code2.x %in% genot_validation,])
      all<-merge(predicted_validation,eblues_validation, by.x = "Juliana.Code2.x", 
                 by.y = "Juliana.Code2")
      all_n<-rbind(all_n,all)
  }
  }
    pcor<-cor(all_n$eblue,all_n$predicted.value, use = "complete.obs")
    mse<-mean((all_n$eblue-all_n$predicted.value)^2, use = "complete.obs")
    metrics_i<-c(hn,j,i,l,pcor,mse)
    r1 = rbind(r1,metrics_i)
  }
  setwd("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Unrelated/DaRT")
  write.csv(r1,paste0(hn,'_AGP17_unrelated_100', '.csv'),quote=FALSE)
}

###Load pedigree file, estimate Amatrix with AGHmatrix package and then curate to retain same families as analyzed with genotypic data
PedigreeAGP17<-read.csv("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/PedigreeAGP17.csv")
Amat<-Amatrix(PedigreeAGP17, ploidy = 4)
head(PedigreeAGP17)
PedigreeAGP17<-PedigreeAGP17[,-1]

Amat<-Amatrix(PedigreeAGP17, ploidy = 4)

Amat2<-Amat
colnames(Amat2)<-paste0("F",colnames(Amat2))
rownames(Amat2)<-paste0("F",rownames(Amat2))
index<-which(rownames(Amat2) %in% rownames(NewXr))
Amat3<-Amat2[index,index]
sorted_row_names <- rownames(Gmatrix)[order(rownames(Gmatrix))]
Gmatrix <- Gmatrix[sorted_row_names, ]
sorted_col_names <- colnames(Gmatrix)[order(colnames(Gmatrix))]
Gmatrix <- Gmatrix[,sorted_col_names ]
sorted_row_names <- rownames(Amat3)[order(rownames(Amat3))]
Amat3 <- Amat3[sorted_row_names, ]
sorted_col_names <- colnames(Amat3)[order(colnames(Amat3))]
Amat3 <- Amat3[,sorted_col_names ]

#### Cross-validations (to test prediction of totally unrelated testing and training populations 
#### we will predict a full-sib family and any family including the male or female from that specific 
#### full-sib family will be eliminated from the training population)


####The k loop will retrieve phenotypic data from each harvest from 1 to 11
####The i loop will select a specific male 
####The j loop will select a specific female 


males<-unique(male.info$Male.name)
males<-males[1:6]
data2<-data2[data2$Male.name.x %in% males,]
females<-unique(male.info$Female.name)
females<-females[1:27]

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
  r1 = data.frame("Harvest","Male","Female","Ginv","pcor","mse")
  r1 = r1[-1,]
  all_n<-data.frame()
    for (i in 1:6) {
      for (j in 1:27) {
        m<-males[i]
        fe<-females[j]
        index_males<-which(data2$Male.name.x == m)
        index_females<-which(data2$Female.name.x == fe)
        data10fold<-data2
        data10fold$DM_kg_ha[c(index_males,index_females)]<-NA ###This works!!
        
        val <- data10fold[data10fold$Male.name.x == m & data10fold$Female.name.x == fe, ]
        genot_validation<-unique(val$Juliana.Code2.x)
        
        colnames(eblueH1b)[2:3]<-c('eblue', 'std.error_eblue')
        colnames(neweblueH1b)[2:3]<-c('eblue', 'std.error_eblue')
   
        modelharvest1<-asreml(fixed = DM_kg_ha ~ 1,
                              random = ~ Row + Column + vm(Juliana.Code2.x,Amat3)
                              , 
                              na.action = na.method(y='include',x='include'),
                              data=data10fold)
        modelharvest1<-update.asreml(modelharvest1)
        modelharvest1<-update.asreml(modelharvest1)
        
        predicted<-data.frame(predict(modelharvest1,classify = "Juliana.Code2.x")$pvals)
        eblues_validation<-neweblueH1b[neweblueH1b$Juliana.Code2 %in% genot_validation,]
        predicted_validation<-data.frame(predicted[predicted$Juliana.Code2.x %in% genot_validation,])
        all<-merge(predicted_validation,eblues_validation, by.x = "Juliana.Code2.x", 
                   by.y = "Juliana.Code2")
        all_n<-rbind(all_n,all)
     }
    }
    pcor<-cor(all_n$eblue,all_n$predicted.value, use = "complete.obs")
    mse<-mean((all_n$eblue-all_n$predicted.value)^2, use = "complete.obs")
    metrics_i<-c(paste0("H ", hn),i,j,pcor,mse)
    r1 = rbind(r1,metrics_i)
    colnames(r1)<-c("Harvest","Male","Female","pcor","mse")
  setwd("/Users/pablosipowicz/Documents/Forage Breeding Team/AGP17/Unrelated")
  write.csv(r1,paste0(hn,'_AGP17_unrelated_pedigree', '.csv'),quote=FALSE)
}




