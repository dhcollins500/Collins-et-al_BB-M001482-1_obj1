#DavidHCollins
#01 December 2017
#This script is my fifth attempt to use r to analyse my qRT-PCR data, it is different from my other attempts because I have redone
#many of the qPCR biological replicates and included a interplate calibrator which is the same gene for each new qPCR
#This includes the qRTPCR data for PCRs that were run in November
#See 'qRTPCR plate design and results_BB-M001482-1_obj1_v55_2017_11_10' for summary of plate design
#See 'BB-MM001482-1_obj1_qRTPCR_rawdata_v1_20181201 for summary of results

########################################################################################################################

#First you must alwys reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(dplyr)
library(ggplot2)
library("cowplot")


#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor

#Find file using one of two methods (actually four methods in the book but two are beyond me)

library(readr)
qpcrdatafifthtry171201 <- read_csv("U:/Documents/RAnalysis/qRTPCRdata/data/qpcrdatafifthtry171201.csv")
View(qpcrdatafifthtry171201)

#do you have the right data?

glimpse(qpcrdatafifthtry171201)

#filter out the rows for a gene of interest, chymotrypsin2 in this case, and assign it to its own object

Chymotry2<-filter(qpcrdatafifthtry171201,GeneCode == "Chymotry2")

#Is the object correct?

glimpse(Chymotry2)

#make a very simple boxplot to view some of the data and assign it to a fucntion

Chymotry2_bp <- ggplot(Chymotry2, aes(x=Phenotype,y=RelQuan)) + #basic boxplot plotting commands
  geom_boxplot() + #boxplotfunction
  theme_bw() + #change the theme from the default
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) + #change the x-axis order (if it has discreet variables) using the scale_x_discrete() function
  labs(title= "Chymotrypsin 2", x = "Phenotype", y = "Expression Ratio") # add title and labels                                          



#Use these functions for each other gene on the dataframe

Cyt305A1.647578 <- filter(qpcrdatafifthtry171201,GeneCode == "Cyt305A1.647578")
Cyt6A1.649469 <- filter(qpcrdatafifthtry171201,GeneCode == "Cyt6A1.649469")
Cyt6k1.642936 <- filter(qpcrdatafifthtry171201,GeneCode == "Cyt6k1.642936")
Cyt6k1.648995 <- filter(qpcrdatafifthtry171201,GeneCode == "Cyt6k1.648995")
Hex <- filter(qpcrdatafifthtry171201,GeneCode == "Hex")
Kruppel <- filter(qpcrdatafifthtry171201,GeneCode == "Kruppel")
Nos.res.640031 <- filter(qpcrdatafifthtry171201,GeneCode == "Nos.res.640031")
Nos.res.645614 <- filter(qpcrdatafifthtry171201,GeneCode == "Nos.res.645614")
p17.29C <- filter(qpcrdatafifthtry171201,GeneCode == "p17.29C")
plrp2 <- filter(qpcrdatafifthtry171201,GeneCode == "plrp2")
unch.G0271606 <- filter(qpcrdatafifthtry171201,GeneCode == "unch.G0271606")
Takeout <- filter(qpcrdatafifthtry171201,GeneCode == "Takeout")
dnmt3 <- filter(qpcrdatafifthtry171201,GeneCode == "dnmt3")
VHDL <- filter(qpcrdatafifthtry171201,GeneCode == "VHDL")
Yellow <- filter(qpcrdatafifthtry171201,GeneCode == "Yellow")

#then plot each of them as a boxplot

Cyt305A1.647578_bp <- ggplot(Cyt305A1.647578, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Cytochrome p450 305A1 (100647578)", x = "Phenotype", y = "Expression Ratio")

Cyt6A1.649469_bp <- ggplot(Cyt6A1.649469, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Cytochrome p450 6A1 (10064969)", x = "Phenotype", y = "Expression Ratio")


Cyt6k1.642936_bp <- ggplot(Cyt6k1.642936, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Cytochrome p450 6k1 (100642936)", x = "Phenotype", y = "Expression Ratio")

Cyt6k1.648995_bp <- ggplot(Cyt6k1.648995, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) + 
  labs(title= "Cytochrome p450 6k1 (100648995)", x = "Phenotype", y = "Expression Ratio")

Hex_bp <- ggplot(Hex, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Hexamerin", x = "Phenotype", y = "Expression Ratio")

Kruppel_bp <- ggplot(Kruppel, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Kruppel", x = "Phenotype", y = "Expression Ratio")

Nos.res.640031_bp <- ggplot(Nos.res.640031, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Nose resistant to fluoxetine (100640031)", x = "Phenotype", y = "Expression Ratio")

Nos.res.645614_bp <- ggplot(Nos.res.645614, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Nose resistant to fluoxetine (100645614)", x = "Phenotype", y = "Expression Ratio")

p17.29C_bp <- ggplot(p17.29C, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title= "p17/29C protein", x = "Phenotype", y = "Expression Ratio")

plrp2_bp <- ggplot(plrp2, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Pancreatic lipase related protein 2", x = "Phenotype", y = "Expression Ratio")

unch.G0271606_bp <- ggplot(unch.G0271606, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Uncharacterized protein (G0271606)", x = "Phenotype", y = "Expression Ratio")

Takeout_bp <- ggplot(Takeout, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Takeout", x = "Phenotype", y = "Expression Ratio")

dnmt3_bp <- ggplot(dnmt3, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "DNA methyltransferase 3", x = "Phenotype", y = "Expression Ratio")

VHDL_bp <- ggplot(VHDL, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Very High Density Lipoprotein", x = "Phenotype", y = "Expression Ratio")

Yellow_bp <- ggplot(Yellow, aes(x=Phenotype,y=RelQuan)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Yellow", x = "Phenotype", y = "Expression Ratio")

Chymotry2_bp
Cyt305A1.647578_bp
Cyt6A1.649469_bp
Cyt6k1.642936_bp
Cyt6k1.648995_bp
Hex_bp
Kruppel_bp
Nos.res.640031_bp
Nos.res.645614_bp
p17.29C_bp
plrp2_bp 
unch.G0271606_bp
Takeout_bp 
dnmt3_bp 
VHDL_bp
Yellow_bp




#Next we want to do a Mann-Whitney-U test between the two caste phenotypes for each gene, at each Development stage. 
#First use the filter function for each developmental stage and assign them each to a new object

####################plrp2###############################################

Early.plrp2 <- filter (plrp2, Development == "Early")
Medium.plrp2 <- filter (plrp2, Development == "Medium")
Late.plrp2 <- filter (plrp2, Development == "Late")

Early.plrp2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.plrp2)
Medium.plrp2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.plrp2)
Late.plrp2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.plrp2)

Early.plrp2.MWUTest
Medium.plrp2.MWUTest
Late.plrp2.MWUTest

####################p17.29C###############################################

Early.p17.29C <- filter (p17.29C, Development == "Early")
Medium.p17.29C <- filter (p17.29C, Development == "Medium")
Late.p17.29C <- filter (p17.29C, Development == "Late")

Early.p17.29C.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.p17.29C)
Medium.p17.29C.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.p17.29C)
Late.p17.29C.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.p17.29C)

Early.p17.29C.MWUTest
Medium.p17.29C.MWUTest
Late.p17.29C.MWUTest

####################unch.G0271606###############################################

Early.unch.G0271606 <- filter (unch.G0271606, Development == "Early")
Medium.unch.G0271606 <- filter (unch.G0271606, Development == "Medium")
Late.unch.G0271606 <- filter (unch.G0271606, Development == "Late")

Early.unch.G0271606.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.unch.G0271606)
Medium.unch.G0271606.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.unch.G0271606)
Late.unch.G0271606.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.unch.G0271606)

Early.unch.G0271606.MWUTest
Medium.unch.G0271606.MWUTest
Late.unch.G0271606.MWUTest

####################Hex###############################################

Early.Hex <- filter (Hex, Development == "Early")
Medium.Hex <- filter (Hex, Development == "Medium")
Late.Hex <- filter (Hex, Development == "Late")

Early.Hex.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Hex)
Medium.Hex.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Hex)
Late.Hex.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Hex)

Early.Hex.MWUTest
Medium.Hex.MWUTest
Late.Hex.MWUTest

####################Kruppel###############################################

Early.Kruppel <- filter (Kruppel, Development == "Early")
Medium.Kruppel <- filter (Kruppel, Development == "Medium")
Late.Kruppel <- filter (Kruppel, Development == "Late")

Early.Kruppel.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Kruppel)
Medium.Kruppel.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Kruppel)
Late.Kruppel.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Kruppel)

Early.Kruppel.MWUTest
Medium.Kruppel.MWUTest
Late.Kruppel.MWUTest

####################Takeout###############################################

Early.Takeout <- filter (Takeout, Development == "Early")
Medium.Takeout <- filter (Takeout, Development == "Medium")
Late.Takeout <- filter (Takeout, Development == "Late")

Early.Takeout.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Takeout)
Medium.Takeout.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Takeout)
Late.Takeout.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Takeout)

Early.Takeout.MWUTest
Medium.Takeout.MWUTest
Late.Takeout.MWUTest

###################Chymotrypsin 2#################################################

Early.Chymotry2 <- filter (Chymotry2, Development == "Early")
Medium.Chymotry2 <- filter (Chymotry2, Development == "Medium")
Late.Chymotry2 <- filter (Chymotry2, Development == "Late")

#Then do the MWU test on each pair of phenotypes and assign them to a seperate object (for easy reference)

Early.Chymotry2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Chymotry2)
Medium.Chymotry2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Chymotry2)
Late.Chymotry2.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Chymotry2)

Early.Chymotry2.MWUTest 
Medium.Chymotry2.MWUTest
Late.Chymotry2.MWUTest


####################Yellow###############################################

Early.Yellow <- filter (Yellow, Development == "Early")
Medium.Yellow <- filter (Yellow, Development == "Medium")
Late.Yellow <- filter (Yellow, Development == "Late")

Early.Yellow.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Yellow)
Medium.Yellow.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Yellow)
Late.Yellow.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Yellow)

Early.Yellow.MWUTest
Medium.Yellow.MWUTest
Late.Yellow.MWUTest

####################Cyt6A1.649469###############################################

Early.Cyt6A1.649469 <- filter (Cyt6A1.649469, Development == "Early")
Medium.Cyt6A1.649469 <- filter (Cyt6A1.649469, Development == "Medium")
Late.Cyt6A1.649469 <- filter (Cyt6A1.649469, Development == "Late")

Early.Cyt6A1.649469.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt6A1.649469)
Medium.Cyt6A1.649469.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt6A1.649469)
Late.Cyt6A1.649469.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt6A1.649469)

Early.Cyt6A1.649469.MWUTest 
Medium.Cyt6A1.649469.MWUTest
Late.Cyt6A1.649469.MWUTest

####################Cyt6k1.648995###############################################

Early.Cyt6k1.648995 <- filter (Cyt6k1.648995, Development == "Early")
Medium.Cyt6k1.648995 <- filter (Cyt6k1.648995, Development == "Medium")
Late.Cyt6k1.648995 <- filter (Cyt6k1.648995, Development == "Late")

Early.Cyt6k1.648995.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt6k1.648995)
Medium.Cyt6k1.648995.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt6k1.648995)
Late.Cyt6k1.648995.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt6k1.648995)

Early.Cyt6k1.648995.MWUTest 
Medium.Cyt6k1.648995.MWUTest
Late.Cyt6k1.648995.MWUTest

####################Cyt305A1.647578###############################################

Early.Cyt305A1.647578 <- filter (Cyt305A1.647578, Development == "Early")
Medium.Cyt305A1.647578 <- filter (Cyt305A1.647578, Development == "Medium")
Late.Cyt305A1.647578 <- filter (Cyt305A1.647578, Development == "Late")

Early.Cyt305A1.647578.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt305A1.647578)
Medium.Cyt305A1.647578.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt305A1.647578)
Late.Cyt305A1.647578.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt305A1.647578)

Early.Cyt305A1.647578.MWUTest
Medium.Cyt305A1.647578.MWUTest
Late.Cyt305A1.647578.MWUTest


####################Cyt6k1.642936###############################################

Early.Cyt6k1.642936 <- filter (Cyt6k1.642936, Development == "Early")
Medium.Cyt6k1.642936 <- filter (Cyt6k1.642936, Development == "Medium")
Late.Cyt6k1.642936 <- filter (Cyt6k1.642936, Development == "Late")

Early.Cyt6k1.642936.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Cyt6k1.642936)
Medium.Cyt6k1.642936.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Cyt6k1.642936)
Late.Cyt6k1.642936.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Cyt6k1.642936)

Early.Cyt6k1.642936.MWUTest
Medium.Cyt6k1.642936.MWUTest
Late.Cyt6k1.642936.MWUTest 



####################Nos.res.640031###############################################

Early.Nos.res.640031 <- filter (Nos.res.640031, Development == "Early")
Medium.Nos.res.640031 <- filter (Nos.res.640031, Development == "Medium")
Late.Nos.res.640031 <- filter (Nos.res.640031, Development == "Late")

Early.Nos.res.640031.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Nos.res.640031)
Medium.Nos.res.640031.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Nos.res.640031)
Late.Nos.res.640031.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Nos.res.640031)

Early.Nos.res.640031.MWUTest 
Medium.Nos.res.640031.MWUTest
Late.Nos.res.640031.MWUTest

####################Nos.res.645614###############################################

Early.Nos.res.645614 <- filter (Nos.res.645614, Development == "Early")
Medium.Nos.res.645614 <- filter (Nos.res.645614, Development == "Medium")
Late.Nos.res.645614 <- filter (Nos.res.645614, Development == "Late")

Early.Nos.res.645614.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.Nos.res.645614)
Medium.Nos.res.645614.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.Nos.res.645614)
Late.Nos.res.645614.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.Nos.res.645614)

Early.Nos.res.645614.MWUTest
Medium.Nos.res.645614.MWUTest
Late.Nos.res.645614.MWUTest


####################dnmt3###############################################

Early.dnmt3 <- filter (dnmt3, Development == "Early")
Medium.dnmt3 <- filter (dnmt3, Development == "Medium")
Late.dnmt3 <- filter (dnmt3, Development == "Late")

Early.dnmt3.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.dnmt3)
Medium.dnmt3.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.dnmt3)
Late.dnmt3.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.dnmt3)

Early.dnmt3.MWUTest
Medium.dnmt3.MWUTest
Late.dnmt3.MWUTest

####################VHDL###############################################

Early.VHDL <- filter (VHDL, Development == "Early")
Medium.VHDL <- filter (VHDL, Development == "Medium")
Late.VHDL <- filter (VHDL, Development == "Late")

Early.VHDL.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Early.VHDL)
Medium.VHDL.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Medium.VHDL)
Late.VHDL.MWUTest <- wilcox.test(RelQuan ~Phenotype, data = Late.VHDL)

Early.VHDL.MWUTest
Medium.VHDL.MWUTest
Late.VHDL.MWUTest


###################Poster Figure#########################################

#We want a four graph panel, Takeout then p1729.c on the upper, and cytochrome 6k1 and uncharacterised on the lower

Takout_bp2 <- ggplot(Takeout, aes(x=Phenotype,y=RelQuan, fill = Development)) + #aesthetics, fill = Development
  #to fill different colours by Developmental stage
  #next plot the boxplot
  geom_boxplot() +
  #Change to 'classic theme'
  theme_classic() +
  #fill in the boxes using the orange fill
  scale_fill_brewer(palette = "Oranges") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank(), #remove x axis elements, and 
        axis.text.x=element_blank(), # remove y axis title
        axis.text=element_text(size=20)) + #change the font size of the axis to 20
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Takeout")

#Do same fo p17 plot

p17.29C_bp2 <- ggplot(p17.29C, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Oranges") +
  theme(legend.position="none",
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x=element_blank(),
        axis.text=element_text(size=20)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title= "p17/29C protein", x = "Phenotype", y = "Expression Ratio")

#and the same for cytochrome 6k1

Cyt6k1.648995_bp2 <- ggplot(Cyt6k1.648995, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Oranges") +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.x=element_blank(), 
        axis.text=element_text(size=20)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title= "Cytochrome p450 6k1 (100648995)", x = "Phenotype", y = "Expression Ratio")


#and the same for uncharacterised

unch.G0271606_bp2 <- ggplot(unch.G0271606, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Oranges") +
  theme(legend.position="none",
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x=element_blank(), 
        axis.text=element_text(size=20)) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW"))+
  labs(title= "Uncharacterized protein (G0271606)", x = "Phenotype", y = "Expression Ratio")

#And for plrp2

plrp2_bp2 <- ggplot(plrp2, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Oranges") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(x = "Phenotype", y = "Expression Ratio")


#If you have loaded cowplot you can use the following function to arrange them side by side

plot_grid(Takout_bp2,p17.29C_bp2, Cyt6k1.648995_bp2,unch.G0271606_bp2,
          ncol = 2, nrow = 2, labels = "AUTO")

#Here is another way to do it with the gridextra function in cowplot


FourPanelWithLabels <- gridExtra::arrangeGrob(Takout_bp2,p17.29C_bp2, Cyt6k1.648995_bp2,unch.G0271606_bp2,
                                              ncol=2,
                                              bottom=grid::textGrob("Phenotype",gp=gpar(fontsize=20)), #place a label at the bottom, and make the font size 20 points
                                              left=grid::textGrob("Relative Quantification", rot=90,gp=gpar(fontsize=20))) #place a font 20 label at side and rotate it
grid::grid.newpage() #clear current graphic and start a new page (not sure why this is necessary)
grid::grid.draw(FourPanelWithLabels) #draw your plot on the newpage

#############################Paper Figure###########################################################################

### I quite like the aethetics of the poster figure, and I would like to put all the barplots in a single 4x4 panel
### for the paper, therefore I'm going to replicate the poster figure and use it for each gene


plrp2_bp3 <- ggplot(plrp2, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "A", x = "Phenotype", y = "Expression Ratio")


p17.29C_bp3 <- ggplot(p17.29C, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "B", x = "Phenotype", y = "Expression Ratio")



unch.G0271606_bp3 <- ggplot(unch.G0271606, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "C", x = "Phenotype", y = "Expression Ratio")


Hex_bp3 <- ggplot(Hex, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "D", x = "Phenotype", y = "Expression Ratio")

Kruppel_bp3 <- ggplot(Kruppel, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) + 
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "E", x = "Phenotype", y = "Expression Ratio")


Takeout_bp3 <- ggplot(Takeout, aes(x=Phenotype,y=RelQuan, fill = Development)) + #aesthetics, fill = Development
  #to fill different colours by Developmental stage
  #next plot the boxplot
  geom_boxplot() +
  #Change to 'classic theme'
  theme_classic() +
  #fill in the boxes using the orange fill
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank()) + #remove x axis title, 
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "F", x = "Phenotype", y = "Expression Ratio")


Chymotry2_bp3 <- ggplot(Chymotry2, aes(x=Phenotype,y=RelQuan, fill = Development)) + #basic boxplot plotting commands
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) + 
  labs(title= "G", x = "Phenotype", y = "Expression Ratio") # add title and labels

Yellow_bp3 <- ggplot(Yellow, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "H", x = "Phenotype", y = "Expression Ratio")


Cyt6A1.649469_bp3 <- ggplot(Cyt6A1.649469, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "I", x = "Phenotype", y = "Expression Ratio")

Cyt6k1.648995_bp3 <- ggplot(Cyt6k1.648995, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "J", x = "Phenotype", y = "Expression Ratio")

Cyt305A1.647578_bp3 <- ggplot(Cyt305A1.647578, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "K", x = "Phenotype", y = "Expression Ratio")


Cyt6k1.642936_bp3 <- ggplot(Cyt6k1.642936, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "L", x = "Phenotype", y = "Expression Ratio")


Nos.res.640031_bp3 <- ggplot(Nos.res.640031, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) + 
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "M", x = "Phenotype", y = "Expression Ratio")

Nos.res.645614_bp3 <- ggplot(Nos.res.645614, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "N", x = "Phenotype", y = "Expression Ratio")


dnmt3_bp3 <- ggplot(dnmt3, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "O", x = "Phenotype", y = "Expression Ratio")

VHDL_bp3 <- ggplot(VHDL, aes(x=Phenotype,y=RelQuan, fill = Development)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_brewer(palette = "Greys") +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank()) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +  labs(title= "P", x = "Phenotype", y = "Expression Ratio")


plrp2_bp3
p17.29C_bp3
unch.G0271606_bp3
Hex_bp3
Kruppel_bp3
Takeout_bp3
Chymotry2_bp3
Yellow_bp3
Cyt6A1.649469_bp3
Cyt6k1.648995_bp3
Cyt305A1.647578_bp3
Cyt6k1.642936_bp3
Nos.res.640031_bp3
Nos.res.645614_bp3
dnmt3_bp3 
VHDL_bp3



# Simple way to draw it is using plot_grid() function


qRTPCR.Simple.Figure<-plot_grid(plrp2_bp3,unch.G0271606_bp3,p17.29C_bp3,Hex_bp3,
                 Kruppel_bp3,Takeout_bp3,Chymotry2_bp3,Yellow_bp3,
                 Cyt6A1.649469_bp3,Cyt6k1.648995_bp3,Cyt305A1.647578_bp3,Cyt6k1.642936_bp3,
                 Nos.res.640031_bp3, Nos.res.645614_bp3,dnmt3_bp3,VHDL_bp3,
                 ncol = 4, nrow = 4, labels = "AUTO")

qRTPCR.Simple.Figure


# More complicated way allows us to add x and y labels

SixteenPanelWithLabels <- gridExtra::arrangeGrob(plrp2_bp3,p17.29C_bp3,unch.G0271606_bp3,Hex_bp3,
                                                 Kruppel_bp3,Takeout_bp3,Chymotry2_bp3,Yellow_bp3,
                                                 Cyt6A1.649469_bp3,Cyt6k1.648995_bp3,Cyt305A1.647578_bp3,Cyt6k1.642936_bp3,
                                                 Nos.res.640031_bp3, Nos.res.645614_bp3,dnmt3_bp3,VHDL_bp3,
                                                 ncol = 4,
                                                 bottom=grid::textGrob("Phenotype"), #place a label at the bottom, and make the font size 20 points
                                                 left=grid::textGrob("Relative Quantification", rot=90)) #place a font 20 label at side and rotate it
grid::grid.newpage() #clear current graphic and start a new page (not sure why this is necessary)
grid::grid.draw(SixteenPanelWithLabels) #draw your plot on the newpage




