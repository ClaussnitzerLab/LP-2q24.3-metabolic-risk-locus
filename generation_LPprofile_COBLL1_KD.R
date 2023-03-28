#####Generation of LP profile of KD COBLL1####
#####set working directory###
setwd("/Users/sophiestrobel/Desktop/LP/LP data import")

####load LP data set#####
KD_data<-read.csv(file = "KD_COBLL1_PAC_data.csv")
rownames(KD_data)<-KD_data$X
data<-as.data.frame(KD_data[,-1])
############
#####subset individuals#####
PAC1<-subset(data, data$Metadata_patient_number =="603")
PAC2<-subset(data, data$Metadata_patient_number =="589")
#####subset days of differentiation#####
#day14<-subset(PAC1, PAC1$Metadata_Day == "14" )
#day9<-subset(PAC1, PAC1$Metadata_Day == "9" )
#day3<-subset(PAC1, PAC1$Metadata_Day == "3" )
#day0<-subset(PAC1, PAC1$Metadata_Day == "0" )

day14<-subset(PAC2, PAC2$Metadata_Day == "14" )
#day9<-subset(PAC2, PAC2$Metadata_Day == "9" )
#day3<-subset(PAC2, PAC2$Metadata_Day == "3" )
#day0<-subset(PAC2, PAC2$Metadata_Day == "0" )

########input########
###########test##########
AP_subset<-day14
#####subset in meta and feature data#############
meta<-AP_subset[, 1:10]
input<-AP_subset[, 11:1164]
#####filter data, remove blocklisted features, SmallBODIPY objects, DNA features and features with 0´s or NA´s across all samples#####
#####see methods LP manuscripte Laber and Strobel et al.#####
variables<-input

variables[,] <- sapply(sapply(variables[,], as.character),as.numeric) 
back<-as.data.frame(t(variables))

mean<-rowMeans(back)
var<-cbind(mean, back)

var.1<-subset(var, var$mean != "0")
var.1a<- var.1[complete.cases(var.1), ]

variables<-as.data.frame(t(var.1a[,2:7]))

##########combine meta and feature data#######
df<-cbind(meta, variables)

#########split into case and control########
control<- subset(df, si_type=='siNT')

case<- subset(df, si_type=='siCOBLL')

#######numeric matrix control########
mcluster<-control[,11:1164]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcontrol<-as.matrix(mcluster)

###############numeric matrix case ########
mcluster<-case[,11:1164]

mcluster[,] <- sapply(sapply(mcluster[,], as.character),as.numeric) 

matrixcase<-as.matrix(mcluster)

#####t-test############
####define features#######

colvector<-as.character(colnames(matrixcase))

################t-test
outtest<-list()

for(i in colvector){
  
  x<-matrixcase[,i]
  y<-matrixcontrol[,i]
  
  outtest[[i]]<-t.test(x,y)
}


pvalue <- data.frame(matrix(unlist(outtest), nrow=length(outtest), byrow=T))
rownames(pvalue)<-colnames(colvector)
colnames(pvalue)<-c("t-test", "df", "pvalue", "confin1", "confin2", "meanx", "meany", "diffmean", "stderr", "alternative", "method", "dataname")

#######adjust for multiple testing using FDR######
library(qvalue)
p<- as.numeric(as.character(pvalue$pvalue))

qvalues<-qvalue(p)
hist(qvalues)

#######adjust for multiple testing using bonferroni#####
q<-qvalues[["qvalues"]]

padj<--log10(p.adjust(p, method="bonferroni", n=length(p)))

######combine adjusted data#####
ppadj<-as.data.frame(cbind(p,padj, pvalue, q))

####rearrange data for plotting######
volcanomatrix<-ppadj[, c("t-test", "p", "padj", "q")]

colnames(volcanomatrix)<-c("t-test", "p", "nglog10ajp", "qvalue")

volcanomatrix[,] <- sapply(sapply(volcanomatrix[,], as.character),as.numeric) 

feature<-colvector
volcano<-as.data.frame(cbind(volcanomatrix, feature))
rownames(volcano)<-feature

###########plot  volcano######
volcano <- volcano %>%
  dplyr::mutate(feature_color = 
                  ifelse(grepl("BODIPY", rownames(volcano)),'bodipy_feature', 
                         ifelse(grepl("AGP", rownames(volcano)),'AGP_feature', 
                                ifelse(grepl("Mito", rownames(volcano)),'Mito_feature', 
                                              'other_feature'))))

volcano$feature_color<-factor(volcano$feature_color,
                              levels = c("Mito_feature","AGP_feature","bodipy_feature", "other_feature"))

volcano<- volcano[complete.cases(volcano), ]


volcano <- volcano %>%
  dplyr::mutate(feature_sig = 
                  ifelse(volcano$q <= 0.05 & 
                         volcano$p<= 0.05,'5%FDR',
                         'not'))

library (ggplot2)

plot <- ggplot(volcano, aes(x=`t-test`, y=-log10(p))) +
  geom_point(aes(color=feature_color,
                 size = feature_sig,
                 alpha = feature_sig ))+
  xlab("t-statistics") +
  ylab("-log10 p-value") +
  ylim(0,5)+
  xlim(-30,50)+
  scale_size_manual(name = "",
                    values = c("5%FDR" = 1.5, "not" = 1))+
  scale_alpha_manual(name = "",
                     values = c("5%FDR" = 0.7, "not" = 0.05))+
  scale_color_manual(name = "", 
                     values = c("Mito_feature" = "#f56464" ,
                                "AGP_feature" = "#f2cf41",
                                "bodipy_feature" = "#4dac26",
                                "other_feature" = "#959ca3")) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black",
                                        fill = "#fdfff4"))

#########subset significant features##########
significant_features<-subset(volcano, volcano$SNP.pvalue < 0.05 & volcano$q < 0.05)

##########save file#########
write.table(volcano, file = "siCOBll1 vs siNT profile subq D14.csv", sep = ",", col.names = NA,
            qmethod = "double")




