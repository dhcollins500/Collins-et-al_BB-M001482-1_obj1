#DavidHCollins
#30 07 2020

#This script is relevant to the following publication:

#David H. Collins1*, Anders Wirén1, 2*, Marjorie Labédan1, 3, Michael Smith1, David C. Prince1, Irina Mohorianu1, 4, Tamas Dalmay1, and Andrew F. G. Bourke1. Gene expression in larval caste determination of bumblebees and the transition to advanced eusociality.

#This script is to analyse the qRT-PCR data using MWU tests and to produce the a figure showing the qRT-PCR and mRNA-seq for the genes of interest side by side.

#Files required:
#Targetgenes.csv - the mRNAseq data for each of the 'genes of interest' (see manuscript for details)

#qpcrdatafifthtry171201.csv - the qRT-PCR data for each of the genes of interest. The analysis which determined the 'relative quantification' values for each gene is included in the excel file. In addition, the information on the selected target genes is available in "DC_gnorm_script"


###############################################################################################

#First you must alwys reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr)
library(tidyverse)
library(grid)
library(gridExtra)

#Set working directory
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Academic writing/My papers/Collins et al_BB-M001482-1_obj1 paper/R scripts/")


###############################################################################################
###############################################################################################

#Import the RNAseq data for the target genes
Targetgenes <- read_csv("Targetgenes.csv")
View(Targetgenes)

#Import the qRT-PCR data for the target genes
qpcrdatafifthtry171201 <- read_csv("qpcrdatafifthtry171201.csv")
View(qpcrdatafifthtry171201)


####################### Use Mann Whitney U tests to test significance of qRT-PCR Data#######################

#The following shows individual MWU tests for each pair of caste-pathway phenotypes in each developmental stage.  This information is required for plotting the qRT-PCR results.

# Note we have had some debate about the need to correct for multiple comparisons with these data. For practical reasons we can not correct for multiple comparisons because there are so many comparisons and the sample sizes for each qRT-PCR are small (n=6), so it would be impossible to establish significance this way. Scaling up the qRT-PCR sample numbers to achieve the requisite power would have been unfeasible with these data.

# In theory each comparison is an independent test of a hypothesis generated from the mRNA-seq results. I.e. if a gene is differentially expressed between caste phenotypes in a developmental phenotype in mRNA-seq, then the hypothesis is that the same gene will be differentially expressed in the same phenotype according to qRT-PCR, and that it will not be significcantly differentially expressed when the mRNA-seq results are not significantlyy differentially expressed.

#This way we can accept that each qRT-PCR MWU test is independent from every other

####################plrp2###############################################

Early.plrp2 <- filter (qpcrdatafifthtry171201, GeneCode== "plrp2", Development == "Early")
Medium.plrp2 <- filter (qpcrdatafifthtry171201, GeneCode== "plrp2",  Development == "Medium")
Late.plrp2 <- filter (qpcrdatafifthtry171201, GeneCode== "plrp2",  Development == "Late")

Early.plrp2.MWUTest <-  wilcox.test(RelQuan ~ Phenotype, data=Early.plrp2)
Medium.plrp2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.plrp2)
Late.plrp2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.plrp2)

Early.plrp2.MWUTest
Medium.plrp2.MWUTest
Late.plrp2.MWUTest

####################p17.29C###############################################

Early.p17.29C <- filter (qpcrdatafifthtry171201, GeneCode== "p17.29C",  Development == "Early")
Medium.p17.29C <- filter (qpcrdatafifthtry171201, GeneCode== "p17.29C", Development == "Medium")
Late.p17.29C <- filter (qpcrdatafifthtry171201, GeneCode== "p17.29C", Development == "Late")

Early.p17.29C.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.p17.29C)
Medium.p17.29C.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.p17.29C)
Late.p17.29C.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.p17.29C)

Early.p17.29C.MWUTest
Medium.p17.29C.MWUTest
Late.p17.29C.MWUTest

####################unch.G0271606###############################################

Early.unch.G0271606 <- filter (qpcrdatafifthtry171201, GeneCode== "unch.G0271606", Development == "Early")
Medium.unch.G0271606 <- filter (qpcrdatafifthtry171201, GeneCode== "unch.G0271606", Development == "Medium")
Late.unch.G0271606 <- filter (qpcrdatafifthtry171201, GeneCode== "unch.G0271606", Development == "Late")

Early.unch.G0271606.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.unch.G0271606)
Medium.unch.G0271606.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.unch.G0271606)
Late.unch.G0271606.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.unch.G0271606)

Early.unch.G0271606.MWUTest
Medium.unch.G0271606.MWUTest
Late.unch.G0271606.MWUTest

####################Hex###############################################

Early.Hex <- filter (qpcrdatafifthtry171201, GeneCode== "Hex", Development == "Early")
Medium.Hex <- filter (qpcrdatafifthtry171201, GeneCode== "Hex", Development == "Medium")
Late.Hex <- filter (qpcrdatafifthtry171201, GeneCode== "Hex", Development == "Late")

Early.Hex.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Hex)
Medium.Hex.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Hex)
Late.Hex.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Hex)

Early.Hex.MWUTest
Medium.Hex.MWUTest
Late.Hex.MWUTest

####################Kruppel###############################################

Early.Kruppel <- filter (qpcrdatafifthtry171201, GeneCode== "Kruppel", Development == "Early")
Medium.Kruppel <- filter (qpcrdatafifthtry171201, GeneCode== "Kruppel", Development == "Medium")
Late.Kruppel <- filter (qpcrdatafifthtry171201, GeneCode== "Kruppel", Development == "Late")

Early.Kruppel.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Kruppel)
Medium.Kruppel.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Kruppel)
Late.Kruppel.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Kruppel)

Early.Kruppel.MWUTest
Medium.Kruppel.MWUTest
Late.Kruppel.MWUTest

####################Takeout###############################################

Early.Takeout <- filter (qpcrdatafifthtry171201, GeneCode== "Takeout", Development == "Early")
Medium.Takeout <- filter (qpcrdatafifthtry171201, GeneCode== "Takeout", Development == "Medium")
Late.Takeout <- filter (qpcrdatafifthtry171201, GeneCode== "Takeout", Development == "Late")

Early.Takeout.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Takeout)
Medium.Takeout.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Takeout)
Late.Takeout.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Takeout)

Early.Takeout.MWUTest
Medium.Takeout.MWUTest
Late.Takeout.MWUTest

###################Chymotrypsin 2#################################################

Early.Chymotry2 <- filter (qpcrdatafifthtry171201, GeneCode== "Chymotry2", Development == "Early")
Medium.Chymotry2 <- filter (qpcrdatafifthtry171201, GeneCode== "Chymotry2", Development == "Medium")
Late.Chymotry2 <- filter (qpcrdatafifthtry171201, GeneCode== "Chymotry2", Development == "Late")

#Then do the MWU test on each pair of phenotypes and assign them to a seperate object (for easy reference)

Early.Chymotry2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Chymotry2)
Medium.Chymotry2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Chymotry2)
Late.Chymotry2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Chymotry2)

Early.Chymotry2.MWUTest 
Medium.Chymotry2.MWUTest
Late.Chymotry2.MWUTest


####################Yellow###############################################

Early.Yellow <- filter (qpcrdatafifthtry171201, GeneCode== "Yellow", Development == "Early")
Medium.Yellow <- filter (qpcrdatafifthtry171201, GeneCode== "Yellow", Development == "Medium")
Late.Yellow <- filter (qpcrdatafifthtry171201, GeneCode== "Yellow", Development == "Late")

Early.Yellow.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Yellow)
Medium.Yellow.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Yellow)
Late.Yellow.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Yellow)

Early.Yellow.MWUTest
Medium.Yellow.MWUTest
Late.Yellow.MWUTest

####################Cyt6A1.649469###############################################

Early.Cyt6A1.649469 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6A1.649469", Development == "Early")
Medium.Cyt6A1.649469 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6A1.649469", Development == "Medium")
Late.Cyt6A1.649469 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6A1.649469", Development == "Late")

Early.Cyt6A1.649469.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt6A1.649469)
Medium.Cyt6A1.649469.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt6A1.649469)
Late.Cyt6A1.649469.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt6A1.649469)

Early.Cyt6A1.649469.MWUTest 
Medium.Cyt6A1.649469.MWUTest
Late.Cyt6A1.649469.MWUTest

####################Cyt6k1.648995###############################################

Early.Cyt6k1.648995 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6k1.648995", Development == "Early")
Medium.Cyt6k1.648995 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6k1.648995", Development == "Medium")
Late.Cyt6k1.648995 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6k1.648995", Development == "Late")

Early.Cyt6k1.648995.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt6k1.648995)
Medium.Cyt6k1.648995.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt6k1.648995)
Late.Cyt6k1.648995.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt6k1.648995)

Early.Cyt6k1.648995.MWUTest 
Medium.Cyt6k1.648995.MWUTest
Late.Cyt6k1.648995.MWUTest

####################Cyt305A1.647578###############################################

Early.Cyt305A1.647578 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt305A1.647578", Development == "Early")
Medium.Cyt305A1.647578 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt305A1.647578", Development == "Medium")
Late.Cyt305A1.647578 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt305A1.647578", Development == "Late")

Early.Cyt305A1.647578.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt305A1.647578)
Medium.Cyt305A1.647578.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt305A1.647578)
Late.Cyt305A1.647578.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt305A1.647578)

Early.Cyt305A1.647578.MWUTest
Medium.Cyt305A1.647578.MWUTest
Late.Cyt305A1.647578.MWUTest


####################Cyt6k1.642936###############################################

Early.Cyt6k1.642936 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6k1.642936", Development == "Early")
Medium.Cyt6k1.642936 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6k1.642936", Development == "Medium")
Late.Cyt6k1.642936 <- filter (qpcrdatafifthtry171201, GeneCode== "Cyt6k1.642936", Development == "Late")

Early.Cyt6k1.642936.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt6k1.642936)
Medium.Cyt6k1.642936.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt6k1.642936)
Late.Cyt6k1.642936.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt6k1.642936)

Early.Cyt6k1.642936.MWUTest
Medium.Cyt6k1.642936.MWUTest
Late.Cyt6k1.642936.MWUTest 



####################Nos.res.640031###############################################

Early.Nos.res.640031 <- filter (qpcrdatafifthtry171201, GeneCode== "Nos.res.640031", Development == "Early")
Medium.Nos.res.640031 <- filter (qpcrdatafifthtry171201, GeneCode== "Nos.res.640031", Development == "Medium")
Late.Nos.res.640031 <- filter (qpcrdatafifthtry171201, GeneCode== "Nos.res.640031", Development == "Late")

Early.Nos.res.640031.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Nos.res.640031)
Medium.Nos.res.640031.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Nos.res.640031)
Late.Nos.res.640031.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Nos.res.640031)

Early.Nos.res.640031.MWUTest 
Medium.Nos.res.640031.MWUTest
Late.Nos.res.640031.MWUTest

####################Nos.res.645614###############################################

Early.Nos.res.645614 <- filter (qpcrdatafifthtry171201, GeneCode== "Nos.res.645614", Development == "Early")
Medium.Nos.res.645614 <- filter (qpcrdatafifthtry171201, GeneCode== "Nos.res.645614", Development == "Medium")
Late.Nos.res.645614 <- filter (qpcrdatafifthtry171201, GeneCode== "Nos.res.645614", Development == "Late")

Early.Nos.res.645614.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Nos.res.645614)
Medium.Nos.res.645614.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Nos.res.645614)
Late.Nos.res.645614.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Nos.res.645614)

Early.Nos.res.645614.MWUTest
Medium.Nos.res.645614.MWUTest
Late.Nos.res.645614.MWUTest


####################dnmt3###############################################

Early.dnmt3 <- filter (qpcrdatafifthtry171201, GeneCode== "dnmt3", Development == "Early")
Medium.dnmt3 <- filter (qpcrdatafifthtry171201, GeneCode== "dnmt3", Development == "Medium")
Late.dnmt3 <- filter (qpcrdatafifthtry171201, GeneCode== "dnmt3", Development == "Late")

Early.dnmt3.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.dnmt3)
Medium.dnmt3.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.dnmt3)
Late.dnmt3.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.dnmt3)

Early.dnmt3.MWUTest
Medium.dnmt3.MWUTest
Late.dnmt3.MWUTest

####################VHDL###############################################

Early.VHDL <- filter (qpcrdatafifthtry171201, GeneCode== "VHDL", Development == "Early")
Medium.VHDL <- filter (qpcrdatafifthtry171201, GeneCode== "VHDL", Development == "Medium")
Late.VHDL <- filter (qpcrdatafifthtry171201, GeneCode== "VHDL", Development == "Late")

Early.VHDL.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.VHDL)
Medium.VHDL.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.VHDL)
Late.VHDL.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.VHDL)

Early.VHDL.MWUTest
Medium.VHDL.MWUTest
Late.VHDL.MWUTest


######################### Make bar plot grid for presentation of mRNA-seq data and qRT-PCR ################


####Plrp2###


#filter out the rows for a gene of interest, plrp2 
#Make univariate scatter plots for mRNA-seq results. Use the genelist (Figure s9 to determine significance bars)

plrp2.us <- filter(Targetgenes,GeneCode == "plrp2") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 24000, yend = 24000, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=25500, label=">1"), color="dim grey",size=3 )+
  annotate("segment", x = 5, xend = 6, y = 7250, yend = 7250, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=8750, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results - (use the output of the mwu tests to deterrmine the significance bars)
plrp2.bp <- filter(qpcrdatafifthtry171201,GeneCode == "plrp2") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 3, xend = 4, y = 1.87, yend = 1.87, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=2.42, label="**"), color="black",size=5)+
  annotate("segment", x = 5, xend = 6, y = 0.7, yend = 0.7, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=1.25, label="**"), color="black",size=5)
  

plrp2.A <- grid.arrange(plrp2.us,plrp2.bp, nrow=1,
                        top = textGrob(expression(""~bold(A)~""~italic(plrp-2)~""), x = 0, hjust = 0))


#Do the same for the rest

####p17.29C####

p17.29C.us <- filter(Targetgenes,GeneCode == "p17.29C") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 7000, yend = 7000, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=9857, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
p17.29C.bp <- filter(qpcrdatafifthtry171201,GeneCode == "p17.29C") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 3, xend = 4, y = -5.2, yend = -5.2, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=-4.2, label="**"), color="black",size=5)


p17.29C.B <- grid.arrange(p17.29C.us,p17.29C.bp, nrow=1,
                          top = textGrob(expression(""~bold(B)~""~italic(p17.29C-like)~""), x = 0, hjust = 0))



###putative uncharacterized protein DDB_G0271606###


unch.G0271606.us <- filter(Targetgenes,GeneCode == "unch.G0271606") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 5, xend = 6, y = 774, yend = 774, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=832.5, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
unch.G0271606.bp <- filter(qpcrdatafifthtry171201,GeneCode == "unch.G0271606") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 5, xend = 6, y = 9.2, yend = 9.2, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=9.9, label="*"), color="black",size=5)


unch.G0271606.C <- grid.arrange(unch.G0271606.us,unch.G0271606.bp, nrow=1,
                                top = textGrob(expression(""~bold(C)~""~italic(Uncharacterised)~"(XM_003399878.3)"), x = 0, hjust = 0))




####hexamerin####


hex.us <- filter(Targetgenes,GeneCode == "Hex") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 12360, yend = 12360, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 13060, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
hex.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Hex") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 5, xend = 6, y = 2.54, yend = 2.54, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=2.76, label="*"), color="black",size=5)


hex.D <- grid.arrange(hex.us,hex.bp, nrow=1,
                      top = textGrob(expression(""~bold(D)~""~italic(hex)~""), x = 0, hjust = 0))




####Kruppel homolog-1####



kruppel.us <- filter(Targetgenes,GeneCode == "Kruppel") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 740, yend = 740, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 782, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
kruppel.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Kruppel") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 1, xend = 2, y = 3.4, yend = 3.4, colour = "black", size=0.5)+
  geom_text(aes( x=1.5, y=3.72, label="**"), color="black",size=5)+
  annotate("segment", x = 3, xend = 4, y = 3.7, yend = 3.7, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=4.02, label="*"), color="black",size=5)


kruppel.E <- grid.arrange(kruppel.us,kruppel.bp, nrow=1,
                          top = textGrob(expression(""~bold(E)~""~italic(Kr-h1)~""), x = 0, hjust = 0))


####Takeout####


takeout.us <- filter(Targetgenes,GeneCode == "Takeout") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 816, yend = 816, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 916, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
takeout.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Takeout") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))



takeout.F <- grid.arrange(takeout.us,takeout.bp, nrow=1,
                          top = textGrob(expression(""~bold(F)~""~italic(takeout)~""), x = 0, hjust = 0))


####Chymotry2####

Chymotry2.us <- filter(Targetgenes,GeneCode == "Chymotry2") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 5, xend = 6, y = 9500, yend = 9500, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y= 10150, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
Chymotry2.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Chymotry2") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 5, xend = 6, y = 2.09, yend = 2.09, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=2.5, label="*"), color="black",size=5)


chymotry2.G <- grid.arrange(Chymotry2.us,Chymotry2.bp, nrow=1,
                            top = textGrob(expression(""~bold(G)~""~italic(chymotrypsin-2-like)~""), x = 0, hjust = 0))

####Yellow####



yellow.us <- filter(Targetgenes,GeneCode == "Yellow") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 75.2, yend = 75.2, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 81.4, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
yellow.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Yellow") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))



yellow.H <- grid.arrange(yellow.us,yellow.bp, nrow=1,
                         top = textGrob(expression(""~bold(H)~""~italic(yellow)~""), x = 0, hjust = 0))

####Cyt6A1.649469####

Cyt6A1.649469.us <- filter(Targetgenes,GeneCode == "Cyt6A1.649469") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 5, xend = 6, y = 9400, yend = 9400, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y= 10000, label=">1"), color="dim grey",size=3 )


#Make boxplots scatter plots for qRT-PCR results
Cyt6A1.649469.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Cyt6A1.649469") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 1, xend = 2, y = 0.4, yend = 0.4, colour = "black", size=0.5)+
  geom_text(aes( x=1.5, y=0.65, label="*"), color="black",size=5)



Cyt6A1.649469.I <- grid.arrange(Cyt6A1.649469.us,Cyt6A1.649469.bp, nrow=1,
                                top = textGrob(expression(""~bold(I)~""~italic(Cyp6A1)~""), x = 0, hjust = 0))


####Cyt6k1.648995####


#Make scatterplot for mRNA-seq results
Cyt6k1.648995.us <- filter(Targetgenes,GeneCode == "Cyt6k1.648995") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 5, xend = 6, y = 444, yend = 444, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y= 472, label=">1"), color="dim grey",size=3 )


#Make boxplots for qRT-PCR results
Cyt6k1.648995.bp <- filter(qpcrdatafifthtry171201,GeneCode == "Cyt6k1.648995") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 5, xend = 6, y = 4.86, yend = 4.86, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=5.65, label="**"), color="black",size=5)



Cyt6A1.649469.J <- grid.arrange(Cyt6k1.648995.us,Cyt6k1.648995.bp, nrow=1,
                                top = textGrob(expression(""~bold(J)~""~italic(Cyp6k1)~"(XM_012314831.2)"), x = 0, hjust = 0))


####Cyt305A1.647578####


#Make scatterplot for mRNA-seq results
Cyt305A1.647578.us <- filter(Targetgenes,GeneCode == "Cyt305A1.647578") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 1390, yend = 1390, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 1490, label=">1"), color="dim grey",size=3 )


#Make boxplots for qRT-PCR results
Cyt305A1.647578.bp<- filter(qpcrdatafifthtry171201,GeneCode == "Cyt305A1.647578") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 3, xend = 4, y = 1.87, yend = 1.87, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=2.17, label="**"), color="black",size=5)




Cyt305A1.647578.K <- grid.arrange(Cyt305A1.647578.us,Cyt305A1.647578.bp, nrow=1,
                                  top = textGrob(expression(""~bold(K)~""~italic(Cyp305a1)~""), x = 0, hjust = 0))



####Cyt6k1.642936####


#Make scatterplot for mRNA-seq results
Cyt6k1.642936.us <- filter(Targetgenes,GeneCode == "Cyt6k1.642936") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 406, yend = 406, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 448, label=">1"), color="dim grey",size=3 )


#Make boxplots for qRT-PCR results
Cyt6k1.642936.bp<- filter(qpcrdatafifthtry171201,GeneCode == "Cyt6k1.642936") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 5, xend = 6, y = 1.13, yend = 1.13, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=1.505, label="**"), color="black",size=5)





Cyt6k1.642936.L <- grid.arrange(Cyt6k1.642936.us,Cyt6k1.642936.bp, nrow=1,
                                top = textGrob(expression(""~bold(L)~""~italic(Cyp6k1)~"(XM_012314843.2)"), x = 0, hjust = 0))

####Nos.res.640031####


#Make scatterplot for mRNA-seq results
Nos.res.640031.us <- filter(Targetgenes,GeneCode == "Nos.res.640031") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 365, yend = 365, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 390, label=">1"), color="dim grey",size=3 )+
  annotate("segment", x = 5, xend = 6, y = 280, yend = 280, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y= 305, label=">1"), color="dim grey",size=3 )



#Make boxplots for qRT-PCR results
Nos.res.640031.bp<- filter(qpcrdatafifthtry171201,GeneCode == "Nos.res.640031") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 5, xend = 6, y = 3.0, yend = 3.0, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=3.57, label="**"), color="black",size=5)



Nos.res.640031.M <- grid.arrange(Nos.res.640031.us,Nos.res.640031.bp, nrow=1,
                                 top = textGrob(expression(""~bold(M)~""~italic(Nrf-6)~""), x = 0, hjust = 0))


####Nos.res.645614####


#Make scatterplot for mRNA-seq results
Nos.res.645614.us<- filter(Targetgenes,GeneCode == "Nos.res.645614") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")+
  annotate("segment", x = 3, xend = 4, y = 390, yend = 390, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y= 423, label=">1"), color="dim grey",size=3 )+
  annotate("segment", x = 5, xend = 6, y = 550, yend = 550, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y= 583, label=">1"), color="dim grey",size=3 )

#Make boxplots for qRT-PCR results
Nos.res.645614.bp<- filter(qpcrdatafifthtry171201,GeneCode == "Nos.res.645614") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  annotate("segment", x = 3, xend = 4, y = 2.1, yend = 2.1, colour = "black", size=0.5)+
  geom_text(aes( x=3.5, y=2.8, label="**"), color="black",size=5)+
  annotate("segment", x = 5, xend = 6, y = 2.27, yend = 2.27, colour = "black", size=0.5)+
  geom_text(aes( x=5.5, y=2.97, label="**"), color="black",size=5)



Nos.res.645614.N <- grid.arrange(Nos.res.645614.us,Nos.res.645614.bp, nrow=1,
                                 top = textGrob(expression(""~bold(N)~""~italic(Nrf-like)~""), x = 0, hjust = 0))


####dnmt3####

#Make scatterplot for mRNA-seq results
dnmt3.us<- filter(Targetgenes,GeneCode == "dnmt3") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")

#Make boxplots for qRT-PCR results
dnmt3.bp <- filter(qpcrdatafifthtry171201,GeneCode == "dnmt3") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))

Dnmt3.O <- grid.arrange(dnmt3.us,dnmt3.bp, nrow=1,
                        top = textGrob(expression(""~bold(O)~""~italic(Dnmt3)~""), x = 0, hjust = 0))


####VHDL####

#Make scatterplot for mRNA-seq results
Vitellogenin.us<- filter(Targetgenes,GeneCode == "VHDL") %>% 
  ggplot(aes(x=Phenotype,y=ReadNumber, shape=factor(Order))) +
  geom_jitter(width=0.2,size=2,colour="black")+
  scale_shape_manual(values=c(15, 0, 16, 1, 17, 2))+
  theme_classic()+
  theme(plot.title=element_text(size=10),
        legend.position="none", 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title = "mRNA-seq", x = "Phenotype", y = "Normalized read count")+
  scale_y_continuous(position = "left")


#Make boxplots for qRT-PCR results
Vitellogenin.bp <- filter(qpcrdatafifthtry171201,GeneCode == "VHDL") %>% 
  ggplot(aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot(outlier.shape=4) +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  labs(title = "qRT-PCR",x = "Phenotype", y = "Expression Ratio")+
  theme(plot.title=element_text(size=10),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x = element_text(size = 8, angle=45,hjust=1),
        axis.text.y = element_text(size=8)) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))

Vitellogenin.P <- grid.arrange(Vitellogenin.us,Vitellogenin.bp, nrow=1,
                               top = textGrob(expression(""~bold(P)~""~italic(Vitellogenin-like)~""), x = 0, hjust = 0))

#####Plot them all as a Grid################

Gene_results_panel <- grid.arrange(plrp2.A,p17.29C.B,unch.G0271606.C,hex.D,
                  kruppel.E, takeout.F,chymotry2.G,yellow.H,
                  Cyt6A1.649469.I,Cyt6A1.649469.J,Cyt305A1.647578.K,Cyt6k1.642936.L,
                  Nos.res.640031.M,Nos.res.645614.N, Dnmt3.O,Vitellogenin.P, nrow=4)

Gene_results_panel

ggsave("Gene_results_panel.pdf", Gene_results_panel, width = 16, height = 10)

