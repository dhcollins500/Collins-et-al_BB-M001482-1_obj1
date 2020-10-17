#DavidHCollins
#03 08 2020

##This script is relevant to the following publication:

#David H. Collins1*, Anders Wirén1, 2*, Marjorie Labédan1, 3, Michael Smith1, David C. Prince1, Irina Mohorianu1, 4, Tamas Dalmay1, and Andrew F. G. Bourke1. Gene expression in larval caste determination of bumblebees and the transition to advanced eusociality.

#The 1st purpose of this script is to establish if number of novel HDEGs is enriched in worker destined-larvae relative to queen destined larvae

#The 2nd purpose of this script is to establish if the number of unannotated genes is enriched among the DEGs and HDEGs compared to the total geneset

#First you must alwys reset R
rm(list=ls())

#packages
library(tidyverse)

#Set working directory
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/RAnalysis/mRNAseqdata")


######Task 1: establish if number of novel HDEGs is enriched in worker destined-larvae relative to queen destined larvae

#Make a contingency table of the mid larvae 


Mid_instar <-
  matrix(c(14, 92, 7, 54),
         nrow = 2,
         dimnames = list(Gene = c("Novel", "Not Novel"),
                         Caste = c("Worker", "Not Worker")))
fisher.test(Mid_instar)

#There is no enrichment of novel genes in the worker pathway compared to the queen pathway in mid-instar larvae

#Make a contingency table of the late larvae

late_instar <-
  matrix(c(7, 14, 5, 40),
         nrow = 2,
         dimnames = list(Gene = c("Novel", "Not Novel"),
                         Caste = c("Worker", "Not Worker")))
fisher.test(late_instar)


##There is significant enrichment of novel genes in the worker pathway compared to the queen pathway in late-instar larvae


#####Task 2: establish if number of unannotated DEGs/HDEGs is enriched relative to the total numbers of genes in each phenotype


### DEGs ###

unannotated_analysis <- read_csv("unannotated_analysis.csv")

#change regions in the type column to genes instead
unannotated_analysis <- unannotated_analysis %>% 
mutate(type=replace(type, type=="region", "gene"))

#How many unnannotated DEGs, annotated DEGs, unannotated non-DEGs, unannotate non-DEGs there are for each phenotype

summary <- unannotated_analysis %>% 
  group_by(Phenotype,type,DEG) %>% 
  summarise(n=n())


## EQ Fishers test##

(EQ <-  matrix(c(409, 408, 4403, 9597),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))
fisher.test(EQ) 

#i.e unannotated DEGs significantly over represented compared to annotated DEGs

(EW <- matrix(c(200, 1137, 4612, 8868),
                    nrow = 2,
                    dimnames = list(Gene = c("Unannotated", "Annotated"),
                                    Caste = c("DEG", "Not DEG"))))

fisher.test(EW)

#i.e unannotated DEGs significantly under represented compared to annotated DEGs
 

(MQ <- matrix(c(691, 3029, 4121, 6976),
                 nrow = 2,
                 dimnames = list(Gene = c("Unannotated", "Annotated"),
                                 Caste = c("DEG", "Not DEG"))))
fisher.test(MQ)

#i.e unannotated DEGs significantly under represented compared to annotated DEGs


(MW <-  matrix(c(1348, 1224, 3464, 8781),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))
fisher.test(MW)

#i.e unannotated DEGs significantly over represented compared to annotated DEGs


(LQ <-  matrix(c(403, 1573, 4409, 8432),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))
fisher.test(EQ)

#i.e unannotated DEGs significantly under represented compared to annotated DEGs


(LW <-  matrix(c(427, 1005, 4385, 9000),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))
fisher.test(LW)

#i.e unannotated DEGs significantly (borderline) under represented compared to annotated DEGs

#Overall these results provide little or no evidence for the hypothesis that unannoted genes are more likely to be differentially expressed


### HDEGs ####

#How many unnannotated HDEGs, annotated HDEGs, unannotated non-HDEGs, unannotate non-HDEGs are there are for each phenotype

summary_HDEG <- unannotated_analysis %>% 
  group_by(Phenotype,type,HDEG) %>% 
  summarise(n=n())



(MQ_HDEG <- matrix(c(19, 39, 4793, 9966),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))
fisher.test(MQ_HDEG)

#i.e unannotated HDEGs not significantly over/under represented compared to annotated DEGs



(MW_HDEG <- matrix(c(49, 45, 4763, 9960),
                   nrow = 2,
                   dimnames = list(Gene = c("Unannotated", "Annotated"),
                                   Caste = c("DEG", "Not DEG"))))

fisher.test(MW_HDEG)

#i.e unannotated HDEGs significantly over represented compared to annotated DEGs


(L_HDEG <- matrix(c(7, 33, 4805, 9972),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))

fisher.test(LQ_HDEG)

#i.e unannotated HDEGs significantly (borderline) under represented compared to annotated DEGs


(LW_HDEG <- matrix(c(10, 4, 4802, 10001),
              nrow = 2,
              dimnames = list(Gene = c("Unannotated", "Annotated"),
                              Caste = c("DEG", "Not DEG"))))

fisher.test(LW_HDEG)

#i.e unannotated HDEGs significantly (borderline) over represented compared to annotated DEGs


#Overall these results show little support for the hypothesis that unnannoted genes are more likely to be highly differentially expressed than annotated genes

#However they provide some support for the hypothesis that unannotated genes are more likely to be highly upregulated in the worker pathway compared to the queen pathway.

