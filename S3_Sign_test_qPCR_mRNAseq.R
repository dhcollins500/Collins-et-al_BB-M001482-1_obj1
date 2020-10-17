#DavidHCollins
#30 07 2020

#This script is relevant to the following publication:

#David H. Collins1*, Anders Wirén1, 2*, Marjorie Labédan1, 3, Michael Smith1, David C. Prince1, Irina Mohorianu1, 4, Tamas Dalmay1, and Andrew F. G. Bourke1. Gene expression in larval caste determination of bumblebees and the transition to advanced eusociality.

#This script is to comparison tests to show that the qRT-PCR data shows the same results as the mRNA-seq data. We do this in two ways. The first way is to directly identify whether the average expression is greater in queens or in workers for each gene in each developmental phenotype, and using each method (qPCR or mRNAseq). This first analysis therefore, does not account for significant and non-significant differences in expression between queens and workers.

#The second comparison looks at significance. For both comparisons we use signed rank and other non-parametric tests to test for differences between the two methods.


###############################################################################################

#First you must alwys reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr)
library(tidyverse)
library(gridExtra)

#Set working directory
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Academic writing/My papers/Collins et al_BB-M001482-1_obj1 paper/R scripts/")

#Import the RNAseq data for the target genes
Targetgenes <- read_csv("Targetgenes.csv")
View(Targetgenes)

#Import the qRT-PCR data for the target genes
qPCR <- read_csv("qpcrdatafifthtry171201.csv")
View(qpcrdatafifthtry171201)

#############################################################################################

#Create a summary of the gene expression values for each phenotype within each gene in the mRNA-seq data
df <- Targetgenes %>% 
  group_by(GeneCode, Phenotype) %>% 
  summarise(mean_ReadNumber= mean(ReadNumber)) %>% 
  pivot_wider(names_from = Phenotype, values_from = mean_ReadNumber)

#Create a pivot table so that for each comparison (queen vs worker) at each developmental stage, we know whether the gene is greater in the queens than in the workers in the mRNA-seq data

df2 <- df %>% 
  mutate(EQvsEW = ifelse(EQ>EW,T,F)) %>% 
  mutate(MQvsMW = ifelse(MQ>MW,T,F)) %>% 
  mutate(LQvsLW = ifelse(LQ>LW,T,F)) %>% 
  select(-(EQ:MW))%>% 
  pivot_longer(-GeneCode, names_to = "Comparison", values_to = "Q_greater_mRNAseq")
         
#Do the same thing as above for the qRT-PCR data
df3 <- qPCR %>% 
  group_by(GeneCode, Phenotype) %>% 
  summarise(mean_Quant= mean(RelQuan))%>% 
  pivot_wider(names_from = Phenotype, values_from = mean_Quant)

df4 <- df3 %>% 
  mutate(EQvsEW = ifelse(EQ>EW,T,F)) %>% 
  mutate(MQvsMW = ifelse(MQ>MW,T,F)) %>% 
  mutate(LQvsLW = ifelse(LQ>LW,T,F)) %>% 
  select(-(EQ:MW)) %>% 
  pivot_longer(-GeneCode, names_to = "Comparison", values_to = "Q_greater_qPCR")

df5 <- df2 %>% 
  left_join(df4,by= c("GeneCode","Comparison")) %>% 
  mutate(Methods_match = ifelse(Q_greater_mRNAseq==Q_greater_qPCR,T,F))
        
df5 %>% 
  ungroup() %>% 
count(Methods_match)

wilcox.test(df5$Q_greater_mRNAseq, df5$Q_greater_qPCR, paired=TRUE) 

#Signed rank test, 32 successes out of 48 trials, null hypothesis is equal successes and failures, alternative is that the successes are greater than the failures

binom.test(32, n = 48, alternative = "greater")

#McNemar's test - another test I could run


df6 <- df5 %>% 
  ungroup() %>% 
  mutate(mRNAseq = as.numeric(Q_greater_mRNAseq)) %>% 
  mutate(qPCR=as.numeric(Q_greater_qPCR)) %>% 
  select(Q_greater_mRNAseq,Q_greater_qPCR)

cont.tab <- xtabs(data = df6)

mcnemar.test(cont.tab)

###  Neither mRNAseq nor qPCR are more likely to show that a gene is more expressed in queens than workers, specifically because neither the 8 nor the 8 in the table are large relative to the other. This isn't really the right test for this hypothesis (that the mRNAseq and the qPCR are related to each other)

### Instead we want a test of significance
chisq.test(cont.tab)
# Bordeline non-signficant - we should not reject the null hypothesis that qPCR and mRNAseq are not associated with each other - this could be a sample size issue, and is partly explained by the non-significant differences affecting the predictions. Indeed, it's still quite close to 0.05 considering we didn't take that into account.



#########################################################################################################################
#New data, these data are modified from table S13 (the results of the qPCR and mRNAseq together which shows they produce similar results). Data are arranged by GeneCode, Developmental stage, mRNAseq result, qPCR result. For the results greater means queen-destined larvae expression is significantly greater than worke-destined larvae expression, lesser means queen-destined llarvae expression is significantly smaller than worker-destined larvae expression, equal means expression differences between castes are not significantly different from zero. Significance for mRNA is log offset fold change differences greater than 1, and non-overlapping biological replicates.


Both <- read_csv("mRNA_qPCR_comparison.csv")
View(Both)

df7 <- Both %>% 
  mutate(Methods_match = ifelse(mRNA_seq==qRT_PCR,T,F))

df7 %>% 
  ungroup() %>% 
  count(Methods_match)


binom.test(38, n = 48, alternative = "greater")
#binomial test, 38 successes out of 48 trials, null hypothesis is equal successes and failures, alternative is that the successes are greater than the failures

#Andrew prefers to look at differentially expressed and non-differentially expressed results seperately 

#DE results
binom.test(x = c(11,6), p = 0.35 , alternative = "greater")

#Non DE results
binom.test(27, n = 31, alternative = "greater")

?binom.test()




df8 <- df7 %>% 
  select(mRNA_seq,qRT_PCR)

cont.tab2 <- xtabs(data = df8)

Matrix = as.matrix(cont.tab2,
                   header=TRUE,
                   row.names=1)

#Here is an example from the R book on how to analyse these kinds of data

Input =("
Before         Yes.after   No.after   Maybe.after
Yes.before      6          0           1
No.before       5          3           7
Maybe.before   11          1          12
")

Matrix.2 = as.matrix(read.table(textConnection(Input),
                                header=TRUE,
                                row.names=1))

Matrix.2

sum(Matrix.2)

mcnemar.test(Matrix.2)

#Try the same with cont.tab2



Input2 =("
Before         Yes.PCR   No.PCR   Maybe.PCR
Yes.mRNA      9          0           4
No.mRNA       0          2           2
Maybe.mRNA    1          3          27
")


Matrix.3 = as.matrix(read.table(textConnection(Input2),
                                header=TRUE,
                                row.names=1))
mcnemar.test(Matrix.3)
#For some reason this won't work even though it works for the example. Something funny about the numbers. I checked it, and it doesn't like there being two zeros in the table. It can deal with 1 zero but not two.

chisq.test(Matrix.3)
#The chisq test gives a slightly worrying warning so I used the simulate p value argument

chisq.test(Matrix.3, simulate.p.value=TRUE)


#Could try using a Global test for symmetry
library(rcompanion)

nominalSymmetryTest(Matrix.3,
                    method="fdr",
                    digits = 3,
                    MonteCarlo = TRUE,
                    ntrial = 100000)
#This didn't run all values, but it effectively conludes something similar to before, that there is no difference in the method used (qPCR or mRNAseq) in their likelihood of showing whether a given gene is up or down regulated.

#####################################################################################################################
#Another method used by Duncan et al. 2020 was to use correlation coefficients


#Make a table of all the mRNAseq values and all of the qRT_PCR values
df9 <- Targetgenes %>% 
  group_by(GeneCode, Phenotype) %>% 
  summarise(mean_ReadNumber= mean(ReadNumber))

df10 <- qPCR %>% 
  group_by(GeneCode, Phenotype) %>% 
  summarise(mean_Quant= mean(RelQuan)) %>% 
  inner_join(df9,by= c("GeneCode","Phenotype"))
  

df11 <- df10 %>% 
  filter(GeneCode=="plrp2") 


cor_fun <- function(df) cor.test(df$mean_Quant, df$mean_ReadNumber, method = "pearson")


cor_fun(df10)
cor_fun(df11)

#Use paste0 to produce a comma seperated string of gene names
paste0(df$GeneCode, collapse=",")


df10$GeneCode <- factor(df10$GeneCode,levels = c("plrp2","p17.29C", "unch.G0271606","Hex", "Kruppel", "Takeout", "Chymotry2","Yellow","Cyt6A1.649469", "Cyt6k1.648995","Cyt305A1.647578","Cyt6k1.642936","Nos.res.640031","Nos.res.645614","dnmt3","VHDL"))

correlation <- df10 %>% 
  group_by (GeneCode) %>% 
  summarize(cor = cor.test(mean_Quant, mean_ReadNumber, method = "pearson")$estimate,
            df = cor.test(mean_Quant, mean_ReadNumber, method = "pearson")$parameter,
            p = cor.test(mean_Quant, mean_ReadNumber, method = "pearson")$p.value)

correlation_sp <- correlation %>% 
  ggplot(aes(x=GeneCode ,y=cor)) +
  geom_point(width=0.2,size=2)+
  scale_colour_manual(values=c("black"))+
  scale_x_discrete(name = "Gene Name", labels=c("plrp-2", "p17.29C-like", "Uncharacterised (XM_003399878.3)","hex", "Kr-h1", "takeout","chymotrypsin-2-like", "yellow", "Cyp6A1","Cyp6k1 (XM_012314831.2)", "Cyp305a1", "Cyp6k1 (XM_012314843.2)", "Nrf-6", "Nrf-like", "Dnmt3", "Vitellogenin-like"))+
  scale_y_continuous(name="Correlation coefficient")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position="none")


ggsave("Correlation_coefficients.pdf", correlation_sp, width = 8, height = 6)


#As we can see the correlation is not very good for some genes, and the results are similar if we use spearmans or kendalls. Another way might be to correlate the mRNAseq logfold changes with the qPCR differences


foldchange<- read_csv("mRNAseq_foldchanges.csv")

df12 <- df3 %>% 
  mutate(early = EQ-EW) %>% 
  mutate(mid = MQ-MW) %>% 
  mutate(late = LQ-LW) %>% 
  select(-(EQ:MW)) %>% 
  pivot_longer(-GeneCode, names_to = "Phenotype", values_to = "qPCR")%>% 
  inner_join(foldchange,by= c("GeneCode","Phenotype")) %>% 
  select(-(GeneID))

plrp2 <- df12 %>% 
  filter(GeneCode=="plrp2") 

cor_fun2 <- function(df) cor.test(df$fold, df$qPCR, method = "pearson")

cor_fun2(df12)

cor_fun2(plrp2)

correlation2 <- df12 %>% 
  group_by (GeneCode) %>% 
  summarize(cor = cor.test(qPCR, fold, method = "pearson")$estimate,
            df = cor.test(qPCR, fold, method = "pearson")$parameter,
            p = cor.test(qPCR, fold, method = "pearson")$p.value)



correlation2 %>% 
  mutate(low=ifelse(cor<0.7,T,F)) %>% 
  ggplot(aes(x=GeneCode ,y=cor, colour=low)) +
  geom_jitter(width=0.2,size=2)+
  scale_colour_manual(values=c("black", "red"))+
  scale_x_discrete(name="Gene Name")+
  scale_y_continuous(name="Correlation coefficient")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position="none")

correlation2 %>% 
  summarise(mean=mean(cor),
            SD = sd(cor))
#This method shows a tighter correlation for some genes, but not for others, for the paper I will use the first correlation method.

#overall I'm more leaning towards using the discrete measurements for comparing gene expression as these are less likely to be affected by highly variable gene expression in the phenotypes that aren't significantly different from each other.
