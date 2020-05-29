#DavidHCollins
#01 January 2019
#This script is to add all the figures and statistical (but not bioinformatics) analysis for the planned manuscript: Collins, D. H., Wirén, A., Labédan, M., Smith, M., Prince, D.C., Mohorianu, I., Dalmay, T., & Bourke, A. F. G. (2020). Gene expression in larval caste determination of bumblebees and the transition to advanced eusociality. 


###############################################################################################

#First you must alwys reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr)
library(dplyr)
library(ggplot2)
library("cowplot")
library(ReadqPCR)
library(NormqPCR)
library(eulerr) 
library(ggplotify)  

#source("https://bioconductor.org/biocLite.R") - ReadqPCR and NormqPCR may require installation usinng these packages rather than the base R package
#BiocManager::install(c("ReadqPCRa"))


###############################################################################################
###############################################################################################
######################### Make stacked bar plot for figure 1 ##################################

#This script was used to make a stacked bar blot summarising the numbers of each type of annotated and unannotated gene in the lists of Bombus terrestris HDEGs

#Enter table - this is transcriptstatus, which is a summary of all of the genes that were either uncharacterised or unannotated (See manuscript for further details). It contains collumns for:
#1)Gene - the ID of the uncharacterised/unannotated transcript
#2)Phenotype - which of the phenotypes (MQ,MW,LQ,LW - see manuscript for definitions) the gene was upregulated in
#3)Status - the transcripts status as either a Novel coding transcript, Novel non-coding transcript, Unannotated not novel transcript, Annotated transcript - see manuscript for further details
#4)Order - a number given to order the transcripts depending on their status Novel coding transcript = 1, Novel non-coding transcript = 2, Unannotated not novel transcript = 3, Annotated transcript = 4
#5)Age - whether the gene was differentially expressed in an early-, mid- , or late-instar larva


transcriptstatus <- read_csv("Gene,Phenotype,Status,Order,Age
            NW_003568470.1_7_877,LQ,Novel coding transcript,1,late
            NC_015769.1_4820898_4822627,LQ,Novel non-coding transcript,2,late
            NC_015771.1_3152030_3152362,LQ,Novel non-coding transcript,2,late
            NW_003568631.1_1_612,LQ,Novel non-coding transcript,2,late
            NC_015778.1_482243_482837,LQ,Novel non-coding transcript,2,late
            NC_015764.1_10238131_10241059,LQ,Unannotated not novel transcript,3,late
            NC_015763.1_11235629_11239608,LQ,Unannotated not novel transcript,3,late
            NC_015763.1_6338074_6346801,LQ,Unannotated not novel transcript,3,late
            NC_015770.1_12905435_12907153,LQ,Unannotated not novel transcript,3,late
            NW_003566526.1_1902_3818,LQ,Unannotated not novel transcript,3,late
            NC_015764.1_10165373_10181599,LQ,Unannotated not novel transcript,3,late
            NC_015774.1_5145354_5146192,LQ,Unannotated not novel transcript,3,late
            NC_015777.1_4620388_4620973,LQ,Unannotated not novel transcript,3,late
            NW_003568489.1_1_1068,LQ,Unannotated not novel transcript,3,late
            NW_003565676.1_5485_6014,LQ,Unannotated not novel transcript,3,late
            NW_003565557.1_44006_46418,LQ,Annotated transcript,4,late
            NC_015767.1_10428321_10672051,LQ,Annotated transcript,4,late
            NC_015769.1_4790925_4844239,LQ,Annotated transcript,4,late
            NC_015771.1_3137226_3141484,LQ,Annotated transcript,4,late
            NC_015772.1_3772353_3775563,LQ,Annotated transcript,4,late
            NC_015764.1_10254752_10258688,LQ,Annotated transcript,4,late
            NC_015764.1_10241255_10243195,LQ,Annotated transcript,4,late
            NC_015765.1_1212762_1214427,LQ,Annotated transcript,4,late
            NW_003567385.1_551_1034,LQ,Annotated transcript,4,late
            NW_003566949.1_1081_1564,LQ,Annotated transcript,4,late
            NC_015772.1_724225_733475,LQ,Annotated transcript,4,late
            NC_015772.1_5395367_5397934,LQ,Annotated transcript,4,late
            NW_003565557.1_26160_29185,LQ,Annotated transcript,4,late
            NC_015769.1_7559514_7562808,LQ,Annotated transcript,4,late
            NC_015769.1_7563046_7568254,LQ,Annotated transcript,4,late
            NC_015773.1_10266169_10268122,LQ,Annotated transcript,4,late
            NC_015773.1_10270623_10273991,LQ,Annotated transcript,4,late
            NW_003568878.1_25_931,LQ,Annotated transcript,4,late
            NW_003569632.1_26_1122,LQ,Annotated transcript,4,late
            NC_015774.1_8595443_8598918,LQ,Annotated transcript,4,late
            NC_015767.1_10624627_10627639,LQ,Annotated transcript,4,late
            NC_015765.1_12833920_12836958,LQ,Annotated transcript,4,late
            NW_003566548.1_26374_46056,LQ,Annotated transcript,4,late
            NC_015767.1_1054083_1068174,LQ,Annotated transcript,4,late
            NC_015771.1_3045554_3144288,LQ,Annotated transcript,4,late
            NC_015771.1_9108635_9108849,LW,Novel non-coding transcript,2,late
            NW_003565415.1_68987_69283 ,LW,Novel non-coding transcript,2,late
            NW_003565415.1_73539_74047,LW,Novel non-coding transcript,2,late
            NW_003570464.1_639_925,LW,Novel non-coding transcript,2,late
            NC_015764.1_6688167_6688610 ,LW,Novel non-coding transcript,2,late
            NC_015773.1_6836828_6837290 ,LW,Novel non-coding transcript,2,late
            NC_015764.1_12810697_12811120 ,LW,Novel non-coding transcript,2,late
            NC_015773.1_6839121_6841995,LW,Unannotated not novel transcript,3,late
            NW_003566045.1_104220_104817 ,LW,Unannotated not novel transcript,3,late
            NW_003566082.1_14813_15171,LW,Unannotated not novel transcript,3,late
            NW_003570107.1_1_628 ,LW,Unannotated not novel transcript,3,late
            NC_015771.1_3163864_3167718,LW,Annotated transcript,4,late
            NC_015768.1_6580879_6587965,LW,Annotated transcript,4,late
            NW_003566591.1_148357_151275,LW,Annotated transcript,4,late
            NW_003568470.1_7_877,MQ,Novel coding transcript,1,mid
            NC_015763.1_6359773_6362415,MQ,Novel non-coding transcript,2,mid
            NC_015771.1_3152030_3152362 ,MQ,Novel non-coding transcript,2,mid
            NC_015764.1_11089645_11090729 ,MQ,Novel non-coding transcript,2,mid
            NC_015768.1_17011118_17014312 ,MQ,Novel non-coding transcript,2,mid
            NW_003569707.1_46_541 ,MQ,Novel non-coding transcript,2,mid
            NW_003569946.1_6_673 ,MQ,Novel non-coding transcript,2,mid
            NC_015764.1_10238131_10241059,MQ,Unannotated not novel transcript,3,mid
            NC_015776.1_3627415_3629700,MQ,Unannotated not novel transcript,3,mid
            NC_015771.1_1527014_1528641,MQ,Unannotated not novel transcript,3,mid
            NW_003566526.1_1902_3818,MQ,Unannotated not novel transcript,3,mid
            NC_015764.1_10165373_10181599,MQ,Unannotated not novel transcript,3,mid
            NC_015764.1_8677570_8678513 ,MQ,Unannotated not novel transcript,3,mid
            NC_015765.1_2486623_2487371 ,MQ,Unannotated not novel transcript,3,mid
            NC_015765.1_2495099_2495655 ,MQ,Unannotated not novel transcript,3,mid
            NC_015770.1_3659628_3660348,MQ,Unannotated not novel transcript,3,mid
            NC_015774.1_7481321_7481757 ,MQ,Unannotated not novel transcript,3,mid
            NW_003565777.1_231_968 ,MQ,Unannotated not novel transcript,3,mid
            NW_003569368.1_70_910,MQ,Unannotated not novel transcript,3,mid
            Combined_NW_003567937.1_2_424_ NW_003567937.1_449_775,MQ,Unannotated not novel transcript,3,mid
            Combined_NW_003568514.1_1_229_NW_003568514.1_236_572 ,MQ,Unannotated not novel transcript,3,mid
            Combined_NW_003570889.1_1_297_ NW_003570889.1_356_696,MQ,Unannotated not novel transcript,3,mid
            NC_015777.1_2971753_2985774,MQ,Annotated transcript,4,mid
            NW_003566135.1_1387093_1390168,MQ,Annotated transcript,4,mid
            NC_015770.1_7042380_7047282,MQ,Annotated transcript,4,mid
            NC_015764.1_4416927_4424286,MQ,Annotated transcript,4,mid
            NC_015768.1_13169118_13184778,MQ,Annotated transcript,4,mid
            NC_015772.1_14946690_14950657,MQ,Annotated transcript,4,mid
            NC_015769.1_4473080_4479231,MQ,Annotated transcript,4,mid
            NW_003566337.1_93_4609,MQ,Annotated transcript,4,mid
            NC_015768.1_4736229_4741671,MQ,Annotated transcript,4,mid
            NC_015768.1_4742166_4745530,MQ,Annotated transcript,4,mid
            NC_015762.1_7539982_7540770,MQ,Annotated transcript,4,mid
            NC_015771.1_3137226_3141484,MQ,Annotated transcript,4,mid
            NC_015764.1_10241255_10243195,MQ,Annotated transcript,4,mid
            NC_015764.1_10150985_10155037,MQ,Annotated transcript,4,mid
            NC_015764.1_10254752_10258688,MQ,Annotated transcript,4,mid
            NC_015765.1_2490563_2492637,MQ,Annotated transcript,4,mid
            NC_015770.1_9450268_9451697,MQ,Annotated transcript,4,mid
            NC_015777.1_4103831_4107629,MQ,Annotated transcript,4,mid
            NC_015772.1_2431712_2441570,MQ,Annotated transcript,4,mid
            NC_015773.1_5908624_5913137,MQ,Annotated transcript,4,mid
            NW_003565702.1_4661_4969,MQ,Annotated transcript,4,mid
            NC_015769.1_5235044_5245163,MQ,Annotated transcript,4,mid
            NC_015776.1_3633957_3636554,MQ,Annotated transcript,4,mid
            NC_015766.1_8097203_8103084,MQ,Annotated transcript,4,mid
            NC_015772.1_16038359_16041849,MQ,Annotated transcript,4,mid
            NC_015775.1_3270154_3271584,MQ,Annotated transcript,4,mid
            NC_015767.1_2698404_2700314,MQ,Annotated transcript,4,mid
            NC_015770.1_13488361_13493334,MQ,Annotated transcript,4,mid
            NC_015772.1_1091520_1107559,MQ,Annotated transcript,4,mid
            NC_015768.1_15383038_15385558,MQ,Annotated transcript,4,mid
            NC_015770.1_1089710_1091504,MQ,Annotated transcript,4,mid
            NC_015771.1_3045554_3144288,MQ,Annotated transcript,4,mid
            NW_003566532.1_6054_6308 ,MW,Novel non-coding transcript,2,mid
            NW_003566567.1_1_207 ,MW,Novel non-coding transcript,2,mid
            NW_003565415.1_147383_147647,MW,Novel non-coding transcript,2,mid
            NW_003570323.1_319_600,MW,Novel non-coding transcript,2,mid
            NW_003566392.1_12836_13581 ,MW,Novel non-coding transcript,2,mid
            NW_003566392.1_10084_10717 ,MW,Novel non-coding transcript,2,mid
            NW_003566051.1_27612_28072 ,MW,Novel non-coding transcript,2,mid
            NW_003570251.1_1_205 ,MW,Novel non-coding transcript,2,mid
            NW_003566402.1_16205_16626,MW,Novel non-coding transcript,2,mid
            NW_003565578.1_15852_16245 ,MW,Novel non-coding transcript,2,mid
            NW_003566392.1_107329_108415 ,MW,Novel non-coding transcript,2,mid
            NW_003566392.1_117210_118602 ,MW,Novel non-coding transcript,2,mid
            NW_003568569.1_34_1912 ,MW,Novel non-coding transcript,2,mid
            NW_003566310.1_2148_2480 ,MW,Novel non-coding transcript,2,mid
            NW_003570241.1_212_817,MW,Unannotated not novel transcript,3,mid
            NW_003569590.1_181_1803,MW,Unannotated not novel transcript,3,mid
            NW_003566591.1_14033_16443,MW,Unannotated not novel transcript,3,mid
            NC_015773.1_9147329_9148959,MW,Unannotated not novel transcript,3,mid
            NW_003566126.1_931407_932585,MW,Unannotated not novel transcript,3,mid
            NC_015776.1_10821399_10824165,MW,Unannotated not novel transcript,3,mid
            NC_015770.1_14064922_14069970,MW,Unannotated not novel transcript,3,mid
            NC_015770.1_14055769_14060379,MW,Unannotated not novel transcript,3,mid
            NC_015775.1_1188284_1191363,MW,Unannotated not novel transcript,3,mid
            NC_015771.1_1524712_1526197,MW,Unannotated not novel transcript,3,mid
            NW_003566109.1_182735_182941,MW,Unannotated not novel transcript,3,mid
            NW_003566082.1_14813_15171,MW,Unannotated not novel transcript,3,mid
            NW_003566071.1_499_771,MW,Unannotated not novel transcript,3,mid
            NW_003567028.1_105_327 ,MW,Unannotated not novel transcript,3,mid
            NW_003567028.1_801_1078 ,MW,Unannotated not novel transcript,3,mid
            NW_003565683.1_17908_18179 ,MW,Unannotated not novel transcript,3,mid
            NW_003566559.1_125_356,MW,Unannotated not novel transcript,3,mid
            NW_003566428.1_1775_1985,MW,Unannotated not novel transcript,3,mid
            NW_003570107.1_1_628 ,MW,Unannotated not novel transcript,3,mid
            NC_015764.1_6688167_6688610 ,MW,Unannotated not novel transcript,3,mid
            NW_003566082.1_16824_21043 ,MW,Unannotated not novel transcript,3,mid
            NW_003565443.1_986853_987131,MW,Unannotated not novel transcript,3,mid
            NW_003566538.1_1206_1438,MW,Unannotated not novel transcript,3,mid
            NW_003566135.1_62598_62999,MW,Unannotated not novel transcript,3,mid
            NW_003566109.1_179772_180002 ,MW,Unannotated not novel transcript,3,mid
            NC_015770.1_2462213_2462450 ,MW,Unannotated not novel transcript,3,mid
            NW_003566533.1_1559_1963 ,MW,Unannotated not novel transcript,3,mid
            NW_003566430.1_254_471 ,MW,Unannotated not novel transcript,3,mid
            NW_003568322.1_527_834,MW,Unannotated not novel transcript,3,mid
            NW_003566068.1_954_1587,MW,Unannotated not novel transcript,3,mid
            NW_003566044.1_232183_232408 ,MW,Unannotated not novel transcript,3,mid
            NC_015778.1_728331_730698 ,MW,Unannotated not novel transcript,3,mid
            NW_003566548.1_68963_69363 ,MW,Unannotated not novel transcript,3,mid
            NW_003568913.1_2_604,MW,Unannotated not novel transcript,3,mid
            NW_003566548.1_70370_70612 ,MW,Unannotated not novel transcript,3,mid
            NC_015770.1_14734826_14735042,MW,Unannotated not novel transcript,3,mid
            NW_003569469.1_989_1280 ,MW,Unannotated not novel transcript,3,mid
            NW_003566421.1_5_273 ,MW,Unannotated not novel transcript,3,mid
            NW_003565429.1_998_1237 ,MW,Unannotated not novel transcript,3,mid
            NC_015775.1_4401317_4401726 ,MW,Unannotated not novel transcript,3,mid
            NW_003570404.1_512_781,MW,Unannotated not novel transcript,3,mid
            Combined_NW_003566430.1_1342_1550_NW_003566430.1_1560_1882,MW,Unannotated not novel transcript,3,mid
            Combined_ NW_003566484.1_500_835_NW_003566484.1_956_1360,MW,Unannotated not novel transcript,3,mid
            NC_015774.1_1138443_1155576,MW,Annotated transcript,4,mid
            NC_015773.1_8155580_8156578,MW,Annotated transcript,4,mid
            NC_015772.1_6823653_6833868,MW,Annotated transcript,4,mid
            NC_015771.1_1512923_1524641,MW,Annotated transcript,4,mid
            NW_003566082.1_28365_29984,MW,Annotated transcript,4,mid
            NW_003566521.1_318_1901,MW,Annotated transcript,4,mid
            NW_003566402.1_17423_17962,MW,Annotated transcript,4,mid
            NW_003566402.1_40433_42001,MW,Annotated transcript,4,mid
            NC_015772.1_3820761_3823080,MW,Annotated transcript,4,mid
            NW_003565979.1_101208_102105,MW,Annotated transcript,4,mid
            NC_015773.1_11570834_11583050,MW,Annotated transcript,4,mid
            NC_015777.1_3837406_3839960,MW,Annotated transcript,4,mid
            NC_015773.1_12803250_12806027,MW,Annotated transcript,4,mid
            NW_003566126.1_900832_957614,MW,Annotated transcript,4,mid
            NC_015770.1_3430376_3431260,MW,Annotated transcript,4,mid
            NW_003567140.1_2_1602,MW,Annotated transcript,4,mid
            NW_003566591.1_154845_159461,MW,Annotated transcript,4,mid
            NW_003566591.1_148357_151275,MW,Annotated transcript,4,mid
            NC_015778.1_743384_746560,MW,Annotated transcript,4,mid
            NW_003565426.1_918418_918958,MW,Annotated transcript,4,mid
            NC_015763.1_3838250_3840620,MW,Annotated transcript,4,mid
            NC_015773.1_75377_97543,MW,Annotated transcript,4,mid
            NC_015763.1_8068073_8072797,MW,Annotated transcript,4,mid
            NC_015773.1_10231842_10238223,MW,Annotated transcript,4,mid
            NC_015773.1_10241749_10246428,MW,Annotated transcript,4,mid
            NC_015773.1_10248242_10252788,MW,Annotated transcript,4,mid
            NC_015775.1_158184_158770,MW,Annotated transcript,4,mid
            NC_015770.1_14103731_14108442,MW,Annotated transcript,4,mid
            NC_015775.1_1623826_1626494,MW,Annotated transcript,4,mid
            NC_015774.1_8055053_8058493,MW,Annotated transcript,4,mid
            NC_015765.1_8124568_8133102,MW,Annotated transcript,4,mid
            NC_015764.1_14164685_14166014,MW,Annotated transcript,4,mid
            NC_015772.1_3712247_3716247,MW,Annotated transcript,4,mid
            NW_003566591.1_142411_146217,MW,Annotated transcript,4,mid
            NW_003566591.1_180278_183051,MW,Annotated transcript,4,mid", 
              col_types= (cols(
                Phenotype = col_factor(),
                Status = col_character(),
                Order = col_double(),
                Age = col_character()
              )))

View(transcriptstatus)

#Reorder the phenotype factors so they are in the correct order for the a axis.
transcriptstatus$Phenotype <- factor(transcriptstatus$Phenotype,levels = c("MQ", "MW", "LQ", "LW"))

#Proportional bar plot

ggplot(transcriptstatus, #the data you want to use
       aes(x = factor(Phenotype), #x axis
       fill = reorder(Status,Order)))+ #y axis ordered by the "order" column
      geom_bar(position = "fill") + #makes figure a proportional barplot
      scale_fill_grey(name="Status")+ #colour scale and name of the legend
      xlab("Phenotype")+ #x label
      ylab("Proportion") #y label

# Count bar plot

ggplot(transcriptstatus, #the data you want to use
       aes(x = factor (Phenotype), #x axis
           fill = reorder(Status,Order))) + #y axis ordered by the "order" column
        geom_bar() + #type of chart
        scale_fill_grey(name="Status") + #colour scale and name of the legend
  xlab("Phenotype")+ #x label      
  ylab("Count") #y label


###############################################################################################
###############################################################################################
######################### Make bar plot grid for presentation of mRNA-seq data ################

#To make the plots for mRNA-seq you need to make a dataframe from 'targetgenes'. This contains the normalised read number for each sample sent for after bioinformatic analysis. It contains the following variables:

#Genename:The full name of the gene of interest
#GeneCode:The appreviated name of the gene of interest
#GeneNumber: A number given to each gene to easily arrange it in the figure
#GeneID: The genes accession number on NCBI
#Sample: The sample that the normalised read count taken from for that gene
#Phenotype: The phenotype that the gene was expressed in (EQ,EW,MQ,MW,LQ,LW)
#Order: A number given to order the phenotypes
#Development: The developmental stage, i.e. whether it was in an early-, mid- , or late-instar larva
#ReadNumber: The normalised read count for the given gene in the given phenotype and sample number

Targetgenes <- read_csv("GeneName,GeneCode,GeneNumber,GeneID,Sample,Phenotype,Order,Development,ReadNumber
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mEQ3,EQ,1,Early,441.7
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mEQ1,EQ,1,Early,997.8
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mEQ2,EQ,1,Early,3900.4
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mEW3,EW,2,Early,93.9
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mEW2,EW,2,Early,245
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mEW1,EW,2,Early,1055.6
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mMQ1,MQ,3,Medium,12346.3
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mMQ4,MQ,3,Medium,18555.9
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mMQ2,MQ,3,Medium,22061.5
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mMW2,MW,4,Medium,27.8
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mMW1,MW,4,Medium,202.5
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mLQ1,LQ,5,Late,3633.1
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mLQ4,LQ,5,Late,4920
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mLQ2,LQ,5,Late,5359.9
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mLW2,LW,6,Late,65.3
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mLW3,LW,6,Late,109.1
pancreatic lipase-related protein 2-like,plrp2,1,XM_003398387.2,mLW1,LW,6,Late,263.7
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mEQ2,EQ,1,Early,210.5
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mEQ1,EQ,1,Early,558.1
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mEQ3,EQ,1,Early,777.2
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mEW2,EW,2,Early,35.9
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mEW1,EW,2,Early,181.1
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mEW3,EW,2,Early,317.1
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mMQ4,MQ,3,Medium,36.8
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mMQ2,MQ,3,Medium,36.9
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mMQ1,MQ,3,Medium,609.1
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mMW2,MW,4,Medium,50.7
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mMW1,MW,4,Medium,134.7
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mLQ4,LQ,5,Late,0.9
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mLQ1,LQ,5,Late,5.5
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mLQ2,LQ,5,Late,7
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mLW2,LW,6,Late,149.6
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mLW3,LW,6,Late,241.2
putative uncharacterized protein DDB_G0271606,unch.G0271606,2,XM_003399878.2,mLW1,LW,6,Late,700.3
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mEQ2,EQ,1,Early,10.9
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mEQ3,EQ,1,Early,16
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mEQ1,EQ,1,Early,47.2
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mEW1,EW,2,Early,17.1
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mEW2,EW,2,Early,28.5
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mEW3,EW,2,Early,38.9
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mMQ2,MQ,3,Medium,25.6
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mMQ1,MQ,3,Medium,32.6
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mMQ4,MQ,3,Medium,187.7
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mMW2,MW,4,Medium,3162.5
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mMW1,MW,4,Medium,3693.7
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mLQ1,LQ,5,Late,11294.6
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mLQ4,LQ,5,Late,20899.4
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mLQ2,LQ,5,Late,26511
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mLW1,LW,6,Late,9570.2
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mLW3,LW,6,Late,29399.3
P17/29C-like protein DDB_G0287399,p17.29C,3,XM_003398794.2,mLW2,LW,6,Late,42070.2
Hexamerin,Hex,4,XM_003401733.2,mEQ1,EQ,1,Early,5751.2
Hexamerin,Hex,4,XM_003401733.2,mEQ2,EQ,1,Early,6032.1
Hexamerin,Hex,4,XM_003401733.2,mEQ3,EQ,1,Early,6521.5
Hexamerin,Hex,4,XM_003401733.2,mEW2,EW,2,Early,3194.4
Hexamerin,Hex,4,XM_003401733.2,mEW3,EW,2,Early,4505
Hexamerin,Hex,4,XM_003401733.2,mEW1,EW,2,Early,4840.7
Hexamerin,Hex,4,XM_003401733.2,mMQ2,MQ,3,Medium,6385.3
Hexamerin,Hex,4,XM_003401733.2,mMQ1,MQ,3,Medium,7762
Hexamerin,Hex,4,XM_003401733.2,mMQ4,MQ,3,Medium,11397.8
Hexamerin,Hex,4,XM_003401733.2,mMW2,MW,4,Medium,1114.4
Hexamerin,Hex,4,XM_003401733.2,mMW1,MW,4,Medium,1200.5
Hexamerin,Hex,4,XM_003401733.2,mLQ1,LQ,5,Late,2143.1
Hexamerin,Hex,4,XM_003401733.2,mLQ4,LQ,5,Late,3833
Hexamerin,Hex,4,XM_003401733.2,mLQ2,LQ,5,Late,4938.4
Hexamerin,Hex,4,XM_003401733.2,mLW1,LW,6,Late,3383.7
Hexamerin,Hex,4,XM_003401733.2,mLW2,LW,6,Late,4024.7
Hexamerin,Hex,4,XM_003401733.2,mLW3,LW,6,Late,4626.5
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mEQ2,EQ,1,Early,421
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mEQ1,EQ,1,Early,425.7
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mEQ3,EQ,1,Early,448
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mEW3,EW,2,Early,299.2
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mEW1,EW,2,Early,390.1
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mEW2,EW,2,Early,409.1
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mMQ2,MQ,3,Medium,595.7
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mMQ4,MQ,3,Medium,629.5
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mMQ1,MQ,3,Medium,683.7
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mMW2,MW,4,Medium,199.7
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mMW1,MW,4,Medium,203.3
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mLQ2,LQ,5,Late,232.4
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mLQ1,LQ,5,Late,296.1
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mLQ4,LQ,5,Late,296.1
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mLW3,LW,6,Late,156.4
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mLW2,LW,6,Late,165.7
Kruppel homolog-1,Kruppel,5,NM_001280921.1,mLW1,LW,6,Late,232.7
Takeout,Takeout,6,XM_003397243.2,mEQ1,EQ,1,Early,893.8
Takeout,Takeout,6,XM_003397243.2,mEQ2,EQ,1,Early,1372.1
Takeout,Takeout,6,XM_003397243.2,mEQ3,EQ,1,Early,1379.1
Takeout,Takeout,6,XM_003397243.2,mEW3,EW,2,Early,818.5
Takeout,Takeout,6,XM_003397243.2,mEW1,EW,2,Early,1461.9
Takeout,Takeout,6,XM_003397243.2,mEW2,EW,2,Early,1514.8
Takeout,Takeout,6,XM_003397243.2,mMQ2,MQ,3,Medium,389.2
Takeout,Takeout,6,XM_003397243.2,mMQ4,MQ,3,Medium,558.7
Takeout,Takeout,6,XM_003397243.2,mMQ1,MQ,3,Medium,695.2
Takeout,Takeout,6,XM_003397243.2,mMW2,MW,4,Medium,141.1
Takeout,Takeout,6,XM_003397243.2,mMW1,MW,4,Medium,153.8
Takeout,Takeout,6,XM_003397243.2,mLQ4,LQ,5,Late,41.5
Takeout,Takeout,6,XM_003397243.2,mLQ1,LQ,5,Late,57.6
Takeout,Takeout,6,XM_003397243.2,mLQ2,LQ,5,Late,62.1
Takeout,Takeout,6,XM_003397243.2,mLW3,LW,6,Late,30.8
Takeout,Takeout,6,XM_003397243.2,mLW2,LW,6,Late,31.7
Takeout,Takeout,6,XM_003397243.2,mLW1,LW,6,Late,186.9
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mEQ3,EQ,1,Early,2002.5
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mEQ1,EQ,1,Early,3331.8
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mEQ2,EQ,1,Early,3331.8
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mEW3,EW,2,Early,2290.4
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mEW1,EW,2,Early,2507.9
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mEW2,EW,2,Early,3033.6
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mMQ1,MQ,3,Medium,3162.5
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mMQ4,MQ,3,Medium,4058.2
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mMQ2,MQ,3,Medium,4901.3
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mMW1,MW,4,Medium,1925
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mMW2,MW,4,Medium,2674.5
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mLQ4,LQ,5,Late,5359.9
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mLQ1,LQ,5,Late,5827.1
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mLQ2,LQ,5,Late,8744.4
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mLW3,LW,6,Late,2144.9
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mLW2,LW,6,Late,2513.6
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,7,XM_012315471.1,mLW1,LW,6,Late,2532.4
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mEQ2,EQ,1,Early,6.9
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mEQ3,EQ,1,Early,7.7
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mEQ1,EQ,1,Early,8.5
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mEW3,EW,2,Early,12.5
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mEW1,EW,2,Early,26.5
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mEW2,EW,2,Early,37.2
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mMQ1,MQ,3,Medium,47.9
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mMQ2,MQ,3,Medium,70.5
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mMQ4,MQ,3,Medium,83.5
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mMW2,MW,4,Medium,269.8
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mMW1,MW,4,Medium,355.3
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mLQ2,LQ,5,Late,432.2
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mLQ1,LQ,5,Late,587.6
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mLQ4,LQ,5,Late,648.7
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mLW1,LW,6,Late,340.5
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mLW3,LW,6,Late,400.7
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,8,XM_012314843.1,mLW2,LW,6,Late,597.9
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mEQ1,EQ,1,Early,0
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mEQ2,EQ,1,Early,0.8
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mEQ3,EQ,1,Early,1.2
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mEW3,EW,2,Early,0.4
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mEW2,EW,2,Early,1.9
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mEW1,EW,2,Early,2.5
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mMQ2,MQ,3,Medium,5
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mMQ1,MQ,3,Medium,9.4
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mMQ4,MQ,3,Medium,14.9
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mMW2,MW,4,Medium,23.2
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mMW1,MW,4,Medium,28.5
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mLQ1,LQ,5,Late,248.2
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mLQ4,LQ,5,Late,256
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mLQ2,LQ,5,Late,399.2
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mLW1,LW,6,Late,55.7
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mLW3,LW,6,Late,62.4
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,9,XM_012314831.1,mLW2,LW,6,Late,69.1
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mEQ1,EQ,1,Early,146.8
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mEQ3,EQ,1,Early,162.8
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mEQ2,EQ,1,Early,167
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mEW3,EW,2,Early,109
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mEW2,EW,2,Early,150.1
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mEW1,EW,2,Early,167.6
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mMQ2,MQ,3,Medium,352.9
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mMQ4,MQ,3,Medium,967.1
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mMQ1,MQ,3,Medium,1272.8
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mMW2,MW,4,Medium,68.8
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mMW1,MW,4,Medium,113.7
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mLQ2,LQ,5,Late,64.5
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mLQ4,LQ,5,Late,111.5
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mLQ1,LQ,5,Late,112.9
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mLW1,LW,6,Late,128.9
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mLW2,LW,6,Late,133.2
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,10,XM_003396620.2,mLW3,LW,6,Late,178.8
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.640031,11,XM_012319905.1,mEQ1,EQ,1,Early,24.7
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mEQ1,EQ,1,Early,32.5
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mEQ2,EQ,1,Early,36.2
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mEQ3,EQ,1,Early,46.9
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mEW3,EW,2,Early,36.5
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mEW2,EW,2,Early,72.8
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mEW1,EW,2,Early,89.7
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mMQ1,MQ,3,Medium,238.2
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mMQ2,MQ,3,Medium,273.9
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mMQ4,MQ,3,Medium,332.2
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mMW2,MW,4,Medium,35.2
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mMW1,MW,4,Medium,55.3
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mLQ2,LQ,5,Late,183
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mLQ4,LQ,5,Late,217.3
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mLQ1,LQ,5,Late,244.4
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mLW3,LW,6,Late,4.9
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mLW2,LW,6,Late,7.9
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,11,XM_012321236.1,mLW1,LW,6,Late,46.2
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mEQ2,EQ,1,Early,29.1
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mEQ3,EQ,1,Early,31.7
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mEW3,EW,2,Early,19.9
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mEW1,EW,2,Early,78.2
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mEW2,EW,2,Early,79.5
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mMQ2,MQ,3,Medium,195.3
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mMQ1,MQ,3,Medium,257.6
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mMQ4,MQ,3,Medium,324.1
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mMW2,MW,4,Medium,12.8
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mMW1,MW,4,Medium,39.2
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mLQ2,LQ,5,Late,419.3
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mLQ1,LQ,5,Late,441.2
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mLQ4,LQ,5,Late,481
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mLW2,LW,6,Late,13.6
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mLW3,LW,6,Late,15.4
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,12,XM_012319905.1,mLW1,LW,6,Late,40
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mEQ1,EQ,1,Early,406.1
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mEQ3,EQ,1,Early,487.2
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mEQ2,EQ,1,Early,492.1
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mEW1,EW,2,Early,349.6
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mEW2,EW,2,Early,381.3
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mEW3,EW,2,Early,402.3
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mMQ1,MQ,3,Medium,301.9
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mMQ4,MQ,3,Medium,559.1
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mMQ2,MQ,3,Medium,567.1
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mMW1,MW,4,Medium,348.6
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mMW2,MW,4,Medium,696.4
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mLQ1,LQ,5,Late,6843.3
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mLQ2,LQ,5,Late,8269.4
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mLQ4,LQ,5,Late,8626
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mLW1,LW,6,Late,618.2
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mLW2,LW,6,Late,826.3
chymotrypsin 2,Chymotry2,13,XM_012309549.1,mLW3,LW,6,Late,1083.3
Yellow protein,Yellow,14,XM_003399634.2,mEQ2,EQ,1,Early,50.4
Yellow protein,Yellow,14,XM_003399634.2,mEQ1,EQ,1,Early,71.9
Yellow protein,Yellow,14,XM_003399634.2,mEQ3,EQ,1,Early,114.8
Yellow protein,Yellow,14,XM_003399634.2,mEW3,EW,2,Early,72.1
Yellow protein,Yellow,14,XM_003399634.2,mEW2,EW,2,Early,82.6
Yellow protein,Yellow,14,XM_003399634.2,mEW1,EW,2,Early,88.8
Yellow protein,Yellow,14,XM_003399634.2,mMQ4,MQ,3,Medium,16
Yellow protein,Yellow,14,XM_003399634.2,mMQ2,MQ,3,Medium,18.9
Yellow protein,Yellow,14,XM_003399634.2,mMQ1,MQ,3,Medium,22.4
Yellow protein,Yellow,14,XM_003399634.2,mMW1,MW,4,Medium,62.7
Yellow protein,Yellow,14,XM_003399634.2,mMW2,MW,4,Medium,67.1
Yellow protein,Yellow,14,XM_003399634.2,mLQ4,LQ,5,Late,24.2
Yellow protein,Yellow,14,XM_003399634.2,mLQ2,LQ,5,Late,24.3
Yellow protein,Yellow,14,XM_003399634.2,mLQ1,LQ,5,Late,27.4
Yellow protein,Yellow,14,XM_003399634.2,mLW3,LW,6,Late,19.4
Yellow protein,Yellow,14,XM_003399634.2,mLW1,LW,6,Late,28.9
Yellow protein,Yellow,14,XM_003399634.2,mLW2,LW,6,Late,30.2
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mEQ3,EQ,1,Early,69.8
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mEQ1,EQ,1,Early,86.7
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mEQ2,EQ,1,Early,91.3
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mEW2,EW,2,Early,73.4
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mEW3,EW,2,Early,78.3
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mEW1,EW,2,Early,83.4
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mMQ1,MQ,3,Medium,68.5
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mMQ4,MQ,3,Medium,73.3
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mMQ2,MQ,3,Medium,77.9
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mMW2,MW,4,Medium,69.8
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mMW1,MW,4,Medium,71.5
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mLQ2,LQ,5,Late,68.8
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mLQ1,LQ,5,Late,73.5
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mLQ4,LQ,5,Late,74
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mLW1,LW,6,Late,72.6
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mLW2,LW,6,Late,86.2
DNA methyltransferase 3,dnmt3,15,XM_012317511.1,mLW3,LW,6,Late,87.3
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mEQ2,EQ,1,Early,18.2
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mEQ3,EQ,1,Early,20.2
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mEQ1,EQ,1,Early,99.6
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mEW1,EW,2,Early,14.6
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mEW2,EW,2,Early,15.3
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mEW3,EW,2,Early,20
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mMQ2,MQ,3,Medium,31.4
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mMQ1,MQ,3,Medium,68.7
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mMQ4,MQ,3,Medium,287.5
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mMW2,MW,4,Medium,21.2
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mMW1,MW,4,Medium,45.3
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mLQ1,LQ,5,Late,21412.4
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mLQ4,LQ,5,Late,27414.8
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mLQ2,LQ,5,Late,28666
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mLW1,LW,6,Late,9904.6
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mLW2,LW,6,Late,26958.7
Very High Density Lipoprotein,VHDL,16,NM_001331111.1,mLW3,LW,6,Late,43340.6")

View(Targetgenes)
glimpse(Targetgenes)

#filter out the rows for a gene of interest, pancreatic lipase in this case, and assign it to its own object

plrp2.df<- filter(Targetgenes, GeneCode=="plrp2")

#check it's correct

plrp2.df

#Do the same for the rest

unch.G0271606.df<- filter(Targetgenes, GeneCode=="unch.G0271606")

p17.29C.df<- filter(Targetgenes, GeneCode=="p17.29C")

Hex.df<- filter(Targetgenes, GeneCode=="Hex")

Kruppel.df<- filter(Targetgenes, GeneCode=="Kruppel")

Takeout.df<- filter(Targetgenes, GeneCode=="Takeout")

Cyt6A1.649469.df<- filter(Targetgenes, GeneCode=="Cyt6A1.649469")

Cyt6k1.642936.df<- filter(Targetgenes, GeneCode=="Cyt6k1.642936")

Cyt6k1.648995.df<- filter(Targetgenes, GeneCode=="Cyt6k1.648995")

Cyt305A1.647578.df<- filter(Targetgenes, GeneCode=="Cyt305A1.647578")

Nos.res.640031.df<- filter(Targetgenes, GeneCode=="Nos.res.640031")

Chymotry2.df<- filter(Targetgenes, GeneCode=="Chymotry2")

Yellow.df<- filter(Targetgenes, GeneCode=="Yellow")

dnmt3.df<- filter(Targetgenes, GeneCode=="dnmt3")

VHDL.df<- filter(Targetgenes, GeneCode=="VHDL")



#Make univariate scatter plots for each gene

plrp2_us <- ggplot(plrp2.df, aes(x=Phenotype,y=ReadNumber)) +
  theme_classic() +
  geom_point() +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Pancreatic lipase related protein 2", x = "Phenotype", y = "Normalised abundance")

plrp2_us


#########################Median and Ranges Univariate Line Plots#################################################
#################################################################################################################

#Manual addition of medians and ranges

#Median + Ranges plotted out when ranges are manually entered

plrp2m=data.frame(Phenotype=c("EQ","EW","MQ","MW","LQ","LW"), median=c(997.8,245,18555.9,115.15,4920,65.3),
                  lower=c(441.7,93.9,12346.3,27.8,3633.1,65.3), upper=c(3900.4,1055.6,22061.5,202.5,5359.9,263.7))
ggplot() + 
  geom_pointrange(data=plrp2m, mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  width=0.2, size=1, color="blue", fill="white", shape=22) +
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "Pancreatic lipase related protein 2", x = "Phenotype", y = "Normalised abundance")


#Better to use a function to make a separate dataframe for each gene

#Make a function that automates the calculation of the medians, mins and max stats and put them
#into a data frame

Gene.Function <- function(x,f.gene,stat = "med"){
  plrp2.df <- x[x$GeneCode==f.gene,]
  treatments <- c("EQ","EW","MQ","MW","LQ","LW")
  meds <- NULL
  mins <- NULL
  maxs <- NULL
  for(i in 1:length(treatments)){
    meds [i] <- median(plrp2.df$ReadNumber[plrp2.df$Phenotype == treatments[i]])
    mins [i] <- min(plrp2.df$ReadNumber[plrp2.df$Phenotype == treatments[i]])
    maxs [i] <- max(plrp2.df$ReadNumber[plrp2.df$Phenotype == treatments[i]])
  }
  #if(stat == "med"){
  # out <- meds
  #}
  #if(stat == "min"){
  # out <- mins
  #  }
  #if(stat == "max"){
  #   out <- maxs
  #}
  out <- data.frame(Phenotype=treatments, median=meds,
                    lower=mins, upper=maxs)
  out
}


plrp2m <- data.frame(Phenotype=treatments, median=meds,
                     lower=mins, upper=maxs)


#Using this function on each individual gene:


###plrp2###

#First make the dataframe to check you have the right data (not totally necessary) 


Gene.Function(Targetgenes, f.gene = "plrp2")

#Then add the function into the plot in place of the dataframe and it will plot out your desired figure for 
#your current gene of interest

plrp2.lp<-ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "plrp2"), 
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  width=0.2, size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "A", x = "Phenotype", y = "Normalised abundance")



####p17.29C####


Gene.Function(Targetgenes, f.gene = "p17.29C")

p17.29C.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "p17.29C"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "B", x = "Phenotype", y = "Normalised abundance")



####putative uncharacterized protein DDB_G0271606####

Gene.Function(Targetgenes, f.gene = "unch.G0271606")

unch.G0271606.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "unch.G0271606"), 
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  width=0.2, size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "C", x = "Phenotype", y = "Normalised abundance")



####hexamerin####


Gene.Function(Targetgenes, f.gene = "Hex")

Hex.lp<- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Hex"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "D", x = "Phenotype", y = "Normalised abundance")


####Kruppel homolog-1####


Gene.Function(Targetgenes, f.gene = "Kruppel")

Kruppel.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Kruppel"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "E", x = "Phenotype", y = "Normalised abundance")

####Takeout####


Gene.Function(Targetgenes, f.gene = "Takeout")

Takeout.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Takeout"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "F", x = "Phenotype", y = "Normalised abundance")



####Chymotry2####


Gene.Function(Targetgenes, f.gene = "Chymotry2")

Chymotry2.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Chymotry2"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "G", x = "Phenotype", y = "Normalised abundance")

####Yellow####


Gene.Function(Targetgenes, f.gene = "Yellow")

Yellow.lp  <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Yellow"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "H", x = "Phenotype", y = "Normalised abundance")


####Cyt6A1.649469####


Gene.Function(Targetgenes, f.gene = "Cyt6A1.649469")

Cyt6A1.649469.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Cyt6A1.649469"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "I", x = "Phenotype", y = "Normalised abundance")


####Cyt6k1.648995####


Gene.Function(Targetgenes, f.gene = "Cyt6k1.648995")

Cyt6k1.648995.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Cyt6k1.648995"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "J", x = "Phenotype", y = "Normalised abundance")

####Cyt305A1.647578####


Gene.Function(Targetgenes, f.gene = "Cyt305A1.647578")

Cyt305A1.647578.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Cyt305A1.647578"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "K", x = "Phenotype", y = "Normalised abundance")

####Cyt6k1.642936####


Gene.Function(Targetgenes, f.gene = "Cyt6k1.642936")

Cyt6k1.642936.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Cyt6k1.642936"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "L", x = "Phenotype", y = "Normalised abundance")


####Nos.res.640031####


Gene.Function(Targetgenes, f.gene = "Nos.res.640031")

Nos.res.640031.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Nos.res.640031"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "M", x = "Phenotype", y = "Normalised abundance")

####Nos.res.645614####


Gene.Function(Targetgenes, f.gene = "Nos.res.645614")

Nos.res.645614.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "Nos.res.645614"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "N", x = "Phenotype", y = "Normalised abundance")


####dnmt3####


Gene.Function(Targetgenes, f.gene = "dnmt3")

dnmt3.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "dnmt3"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "O", x = "Phenotype", y = "Normalised abundance")

####VHDL####


Gene.Function(Targetgenes, f.gene = "VHDL")

VHDL.lp <- ggplot() + 
  geom_pointrange(data=Gene.Function(Targetgenes, f.gene = "VHDL"),
                  mapping=aes(x=Phenotype, y=median, ymin=upper, ymax=lower), 
                  size=1, color="black", fill="white", shape=175) +
  theme(legend.position="none", #get rid of the legend 
        axis.title.x=element_blank(), #remove x axis title, 
        axis.title.y=element_blank())+
  scale_x_discrete(limits=c("EQ","EW","MQ","MW","LQ","LW")) +
  labs(title= "P", x = "Phenotype", y = "Normalised abundance")

#####Plot them all as a Grid################


lp.panel <- gridExtra::arrangeGrob(plrp2.lp,p17.29C.lp,unch.G0271606.lp,Hex.lp,
                                   Kruppel.lp,Takeout.lp,Chymotry2.lp,Yellow.lp,
                                   Cyt6A1.649469.lp,Cyt6k1.648995.lp,Cyt305A1.647578.lp,Cyt6k1.642936.lp,
                                   Nos.res.640031.lp, Nos.res.645614.lp,dnmt3.lp,VHDL.lp,
                                   ncol = 4,
                                   bottom=grid::textGrob("Phenotype"), #place a label at the bottom, and make the font size 20 points
                                   left=grid::textGrob("Relative Quantification", rot=90)) #place a font 20 label at side and rotate it
grid::grid.newpage() #clear current graphic and start a new page (not sure why this is necessary)
grid::grid.draw(lp.panel) #draw your plot on the newpage



###############################################################################################
###############################################################################################
######################### qRT-PCR reference genes ##########################################

##This script is to assess some B terrestris housekeeping genes for suitability as reference genes. *Note to author - This version includes the version of my qRTPCR data that used the October-November 2017 plates which better accounted for interplate variation and are summarised in datafile 'BB-MM001482-1_obj1_qRTPCR_rawdata_v1_20181201.xlsx'
#BterRawData is the file containing the raw Cq values for each reference gene after qRT-PCR analysis. Each row representss the housekeeping gene of interest, and the sample replicate, where each sample represents a replicate from each of the six phenotypes (EQ, EW, MQ, MW, LQ, LW) that were sequenced in the study. It contains the following variables
#Sample -The sample ID (phenotype and replicate number)
#Detector - The abbreviated name of each of the 8 housekeeping genes that was tested
#Cq- The normalised Cq value of each gene in each sample - at this stage the genes have been normalised between technical replicates (see "qpcrdatafifthtry171201"), and the purpose of these methods is to normalise them between samples from the same phenotype (representing biological replicates)

#Important Note 22/05/2020: the following script up until 1277 no longer appears to work, the following is the script that was used to generate the reference genes is supplied for posterity purposes, however it will no longer generate the same results as when the analysis was first run

#Find file

BterRaw <- read_tsv("Sample	Detector	Cq
EQ3	ATub	25.77
EQ4	ATub	25.55
EQ5	ATub	25.64
EQ6	ATub	24.91
EQ7	ATub	25.51
EQ8	ATub	25.23
EW3	ATub	24.83
EW4	ATub	24.89
EW5	ATub	25.04
EW6	ATub	23.99
EW7	ATub	24.71
EW8	ATub	24.44
MQ3	ATub	24.19
MQ4	ATub	24.28
MQ5	ATub	23.01
MQ6	ATub	23.28
MQ7	ATub	23.04
MQ8	ATub	23.95
MW3	ATub	24.98
MW4	ATub	22.99
MW5	ATub	23.93
MW6	ATub	22.84
MW7	ATub	25.08
MW8	ATub	23.51
LQ3	ATub	24.96
LQ4	ATub	25.29
LQ5	ATub	27.36
LQ6	ATub	24.81
LQ7	ATub	25.23
LQ8	ATub	25.40
LW3	ATub	26.24
LW4	ATub	25.70
LW5	ATub	25.28
LW6	ATub	25.21
LW7	ATub	25.74
LW8	ATub	26.37
EQ3	AK	26.34
EQ4	AK	27.03
EQ5	AK	26.18
EQ6	AK	24.75
EQ7	AK	26.98
EQ8	AK	25.45
EW3	AK	25.26
EW4	AK	28.99
EW5	AK	26.67
EW6	AK	24.33
EW7	AK	25.29
EW8	AK	24.43
MQ3	AK	25.55
MQ4	AK	25.54
MQ5	AK	25.39
MQ6	AK	24.04
MQ7	AK	24.11
MQ8	AK	24.70
MW3	AK	25.87
MW4	AK	24.53
MW5	AK	25.41
MW6	AK	24.83
MW7	AK	26.32
MW8	AK	25.48
LQ3	AK	25.33
LQ4	AK	25.77
LQ5	AK	28.44
LQ6	AK	25.13
LQ7	AK	25.59
LQ8	AK	26.00
LW3	AK	28.42
LW4	AK	26.97
LW5	AK	28.07
LW6	AK	24.40
LW7	AK	25.27
LW8	AK	25.91
EQ3	GAPDH	24.60
EQ4	GAPDH	24.82
EQ5	GAPDH	24.09
EQ6	GAPDH	24.29
EQ7	GAPDH	25.78
EQ8	GAPDH	25.58
EW3	GAPDH	23.12
EW4	GAPDH	27.40
EW5	GAPDH	24.91
EW6	GAPDH	24.19
EW7	GAPDH	24.82
EW8	GAPDH	24.05
MQ3	GAPDH	23.73
MQ4	GAPDH	23.51
MQ5	GAPDH	23.16
MQ6	GAPDH	23.49
MQ7	GAPDH	23.34
MQ8	GAPDH	24.18
MW3	GAPDH	23.32
MW4	GAPDH	22.67
MW5	GAPDH	24.31
MW6	GAPDH	24.00
MW7	GAPDH	25.91
MW8	GAPDH	24.62
LQ3	GAPDH	24.23
LQ4	GAPDH	24.60
LQ5	GAPDH	26.78
LQ6	GAPDH	25.36
LQ7	GAPDH	25.46
LQ8	GAPDH	25.93
LW3	GAPDH	26.62
LW4	GAPDH	24.82
LW5	GAPDH	26.96
LW6	GAPDH	25.67
LW7	GAPDH	26.04
LW8	GAPDH	26.47
EQ3	PPIA	24.90
EQ4	PPIA	24.95
EQ5	PPIA	24.92
EQ6	PPIA	23.94
EQ7	PPIA	25.42
EQ8	PPIA	24.98
EW3	PPIA	24.07
EW4	PPIA	23.74
EW5	PPIA	24.99
EW6	PPIA	23.98
EW7	PPIA	24.50
EW8	PPIA	24.35
MQ3	PPIA	24.38
MQ4	PPIA	24.55
MQ5	PPIA	23.81
MQ6	PPIA	23.78
MQ7	PPIA	24.20
MQ8	PPIA	24.52
MW3	PPIA	23.90
MW4	PPIA	23.53
MW5	PPIA	24.40
MW6	PPIA	23.67
MW7	PPIA	25.33
MW8	PPIA	24.09
LQ3	PPIA	24.99
LQ4	PPIA	25.00
LQ5	PPIA	26.33
LQ6	PPIA	25.37
LQ7	PPIA	25.05
LQ8	PPIA	25.41
LW3	PPIA	26.37
LW4	PPIA	25.06
LW5	PPIA	25.87
LW6	PPIA	25.04
LW7	PPIA	25.54
LW8	PPIA	25.83
EQ3	PA2	31.17
EQ4	PA2	31.45
EQ5	PA2	31.36
EQ6	PA2	30.75
EQ7	PA2	31.74
EQ8	PA2	31.70
EW3	PA2	31.10
EW4	PA2	31.20
EW5	PA2	31.22
EW6	PA2	31.36
EW7	PA2	30.73
EW8	PA2	31.12
MQ3	PA2	31.72
MQ4	PA2	30.91
MQ5	PA2	31.12
MQ6	PA2	30.59
MQ7	PA2	31.89
MQ8	PA2	31.90
MW3	PA2	31.49
MW4	PA2	30.52
MW5	PA2	31.09
MW6	PA2	30.70
MW7	PA2	32.79
MW8	PA2	30.76
LQ3	PA2	31.88
LQ4	PA2	31.43
LQ5	PA2	32.73
LQ6	PA2	32.91
LQ7	PA2	31.48
LQ8	PA2	33.02
LW3	PA2	33.13
LW4	PA2	32.61
LW5	PA2	32.09
LW6	PA2	32.96
LW7	PA2	32.31
LW8	PA2	32.73
EQ3	RP18S	23.16
EQ4	RP18S	23.39
EQ5	RP18S	23.18
EQ6	RP18S	22.26
EQ7	RP18S	22.60
EQ8	RP18S	22.58
EW3	RP18S	22.98
EW4	RP18S	22.62
EW5	RP18S	23.38
EW6	RP18S	22.35
EW7	RP18S	22.28
EW8	RP18S	21.83
MQ3	RP18S	23.80
MQ4	RP18S	23.51
MQ5	RP18S	22.93
MQ6	RP18S	23.12
MQ7	RP18S	22.54
MQ8	RP18S	22.21
MW3	RP18S	23.02
MW4	RP18S	22.81
MW5	RP18S	23.07
MW6	RP18S	22.96
MW7	RP18S	23.77
MW8	RP18S	22.47
LQ3	RP18S	23.83
LQ4	RP18S	23.14
LQ5	RP18S	23.74
LQ6	RP18S	24.03
LQ7	RP18S	22.36
LQ8	RP18S	23.11
LW3	RP18S	24.03
LW4	RP18S	23.48
LW5	RP18S	23.87
LW6	RP18S	23.41
LW7	RP18S	23.19
LW8	RP18S	22.64
EQ3	RPS5a	22.98
EQ4	RPS5a	23.94
EQ5	RPS5a	23.27
EQ6	RPS5a	22.46
EQ7	RPS5a	23.33
EQ8	RPS5a	22.99
EW3	RPS5a	22.78
EW4	RPS5a	22.24
EW5	RPS5a	24.16
EW6	RPS5a	22.40
EW7	RPS5a	22.91
EW8	RPS5a	22.45
MQ3	RPS5a	23.91
MQ4	RPS5a	23.21
MQ5	RPS5a	22.97
MQ6	RPS5a	22.49
MQ7	RPS5a	22.80
MQ8	RPS5a	22.96
MW3	RPS5a	23.45
MW4	RPS5a	22.33
MW5	RPS5a	23.35
MW6	RPS5a	22.74
MW7	RPS5a	23.27
MW8	RPS5a	22.60
LQ3	RPS5a	24.17
LQ4	RPS5a	23.52
LQ5	RPS5a	25.28
LQ6	RPS5a	23.69
LQ7	RPS5a	22.93
LQ8	RPS5a	23.55
LW3	RPS5a	24.46
LW4	RPS5a	23.96
LW5	RPS5a	23.92
LW6	RPS5a	23.33
LW7	RPS5a	23.61
LW8	RPS5a	23.39
EQ3	TBP	30.08
EQ4	TBP	30.35
EQ5	TBP	30.42
EQ6	TBP	29.18
EQ7	TBP	30.14
EQ8	TBP	29.83
EW3	TBP	29.56
EW4	TBP	30.01
EW5	TBP	30.27
EW6	TBP	29.28
EW7	TBP	29.59
EW8	TBP	29.09
MQ3	TBP	30.52
MQ4	TBP	30.48
MQ5	TBP	30.09
MQ6	TBP	29.16
MQ7	TBP	29.55
MQ8	TBP	29.72
MW3	TBP	29.98
MW4	TBP	29.34
MW5	TBP	30.25
MW6	TBP	29.78
MW7	TBP	30.70
MW8	TBP	29.30
LQ3	TBP	31.31
LQ4	TBP	30.76
LQ5	TBP	31.98
LQ6	TBP	31.04
LQ7	TBP	30.27
LQ8	TBP	30.21
LW3	TBP	30.89
LW4	TBP	31.12
LW5	TBP	30.74
LW6	TBP	30.76
LW7	TBP	30.28
LW8	TBP	30.88")

BterRaw
browseVignettes("ReadqPCR")


##############################################################
# Block A: The SDS output file, as obtained from the RTqPCR
# system, were loaded into R using the ReadqPCR packages,
# populating a qPCRBatch object
# These first two block measure the best reference gene using genenorm on raw CT values

BterRawData <- read.csv("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/RAnalysis/qRTPCRdata/normQBterdata/BterRawData.txt")
input.file.BterRaw <-  BterRawData
qPCRBatch.Cq <- read.qPCR(filename = input.file.BterRaw) # this command reads in the file
?read.qPCR()
##############################################################
# Block B: NormqPCR was then used to find the optimal reference
# genes, using selectHKs

# First potential housekeeping genes are set:
hkgs <- c("AK","ATub",
          "GAPDH","PA2",
          "PPIA","RP18S","TBP","RPS5a")

# Then selectHKs is run with log = TRUE since we are using Cq
# values, source("/home/jperkins/NormqPCR/R/deltaDeltaCq.R")

BterRaw.geNorm <-selectHKs (qPCRBatch.Cq [hkgs,], method = "geNorm", Symbols = hkgs, 
                            minNrHK = 2, log = TRUE)

#This concludes that TBP, RPS5a, and RP18S are the best reference genes

############################################################


###############################################################################################
###############################################################################################
######################### qRT-PCR stats and figures ##########################################

#Use qpcrdatafifthtry171201 to produce comparisons between samples in terms of their Cq values. Each row represents the Cq value from a qRTPCR experiment where we assessed gene expression for 8 reference genes, and 17 target genes. Only 16 of these are described in the manuscript, Novel1 was dropped from further analysis.

#Genename:The full name of the gene of interest
#GeneCode:The appreviated name of the gene of interest
#Sample: The sample that the normalised read count taken from for that gene
#Phenotype: The phenotype that the gene was expressed in (EQ,EW,MQ,MW,LQ,LW)
#Development: The developmental stage, i.e. whether it was in an early-, mid- , or late-instar larva
#Order: A number given to arrange the phenotypes in a coherent order
#SamOrder: A number given to arrange the samples within each phenotype (not relevant for this analysis)
#ConstantOrder: A number given to arrange the samples within each experiment (not relevant for this analysis)
#GeneType: Whether the gene was of interest as a reference gene (and therefore expected to be stable between samples), or as a target gene (where the aim was to assess for evidence of differential expression between phenotypes)
#ExperimentNo: A number given for the plate that the sample was run on - not relevant for this analysis
#Rep1: The raw Cq value for technical replicate 1
#Rep2: The raw Cq value for technical replicate 2
#Rep3: The raw Cq value for technical replicate 3
#Mean: The mean of the three technical replicates
#Min: The minimum value of the three technical replicates
#Max: The maximum value of the three technical replicates
#Range: The difference between the maximum and minimum value of the the three technical replicates
#InterPlateCal: The biological replicate that was used as an interplate calibrator (IPC) for multiple plates run with the same gene - therefore allowing us to normalise the Cq values between plates
#Cal1: The raw Cq value for technical replicate 1 of the IPC
#Cal2: The raw Cq value for technical replicate 2 of the IPC
#Cal3: The raw Cq value for technical replicate 3 of the IPC
#CalMean: The mean Cq value for the technical replicates of the IPC
#CalMin: The minimum Cq value for the technical replicates of the IPC
#CalMax: The maximum Cq value for the technical replicates of the IPC	
#CalRange: The difference between the maximum and minimum value of the the three technical replicates of the IPC
#CalHigh: The Cq value of the lowest IPC value for two plates that contained samples to be normalised
#PrimerEfficiency: The calculated primer efficiency for each gene, ideally this will be between 1.9 and 2.1, but is often lower than this if the primers are inefficient.
#NormMean	:  The IPC normalised mean Cq value of the three technical replicates for each biological replicate within in each phenotype. The value was calculated by multiplying the mean Cq value of the target gene for the focal sample by the mean Cq value of the IPC that was run on the same plate as the focal sample, and dividing by the mean Cq value of the lowest of the two mean Cq values for the two IPCs run on two plates that we wished to compare. With respect to the columns in the file this would be: (Mean*CalMean/CalHigh)
#BKI: The calculated Best Keeper Index which was calculated using best keeper software to produce a normalised Cq value using the three most stable reference genes across all replicates and phenotypes. The BKI depends on the Cq values for the three reference genes for the focal RNA sample.
#EquationCal: The normalised Cq value for the comparison sample, in this dataset LQ3 was used as the comparison sample to quantify every other sample against and to produce a mean quantification score.
#BKEquationCal:	The Best Keeper Index of the comparison sample (LQ3).
#BKEfficiencyMean: The average efficiency of the three genes that made up the Best Keeper Index
#pfafflnum: The calculated numerator for equation 1 in Pfaffl. 2001. A New Mathematical Model for Relative Quantification in Real-Time RT-PCR. Nucleic Acids Res 29(9):e45 using the qRT-PCR data from the present study. In this sheet this is calculated as (PrimerEfficiency^(EquationCal-NormMean))
#pfafflden:	The calculated denominator for equation 1 in Pfaffl. 2001. A New Mathematical Model for Relative Quantification in Real-Time RT-PCR. Nucleic Acids Res 29(9):e45 using the qRT-PCR data from the present study. In this sheet this is calculated as (BKEfficiencyMean^(BKEquationCal-BKI))
#numoverden: The value of pfafflnum/pfafflden
#RelQuan: The Relative quantification for each of the target genes in the focal samples of each of the focal phenotypes. This calculated by log2(numoverden), which completes the equation from Pfaffl. 2001. A New Mathematical Model for Relative Quantification in Real-Time RT-PCR. Nucleic Acids Res 29(9):e45. For the purposes of the following statistical calculations and figures, RelQuan is representative of the normalised gene expression for each gene and sample.


qpcrdatafifthtry171201 <- read_csv("GeneName,GeneCode,Sample,Phenotype,Development,Order,SamOrder,ConstantOrder,GeneType,ExperimentNo,Rep1,Rep2,Rep3,Mean,Min,Max,Range,InterPlateCal,Cal1,Cal2,Cal3,CalMean,CalMin,CalMax,CalRange,CalHigh,PrimerEfficiency,NormMean,BKI,EquationCal,BKEquationCal,BKEfficiencyMean,pfafflnum,pfafflden,numoverden,RelQuan
Alpha tubulin (Bter03085; XM_003396121),ATub,EQ3,EQ,Early,1,1,1,Control,120,25.74954033,25.92703438,25.64194489,25.77283986,25.64194489,25.92703438,0.285089493,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.77284,25.20347779,24.95721,26.22237753,1.993,0.63859907,2.019146866,0.316271729,-1.66076
Alpha tubulin (Bter03085; XM_003396121),ATub,EQ4,EQ,Early,1,1,1,Control,120,25.41412544,25.50836754,25.74171448,25.55473582,25.41412544,25.74171448,0.327589035,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.55474,25.71323328,24.95721,26.22237753,1.993,0.71996472,1.420667419,0.506779215,-0.98057
Alpha tubulin (Bter03085; XM_003396121),ATub,EQ5,EQ,Early,1,1,1,Control,120,25.94663811,25.53452492,25.42642403,25.63586235,25.42642403,25.94663811,0.520214081,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.63586,25.41256952,24.95721,26.22237753,1.993,0.68855449,1.748007975,0.393908095,-1.34407
Alpha tubulin (Bter03085; XM_003396121),ATub,EQ6,EQ,Early,1,2,1,Control,17,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,24.9137,24.43307804,24.95721,26.22237753,1.993,1.02421323,3.434853817,0.29818248,-1.74573
Alpha tubulin (Bter03085; XM_003396121),ATub,EQ7,EQ,Early,1,2,1,Control,17,25.59929657,25.44042397,25.437006,25.49224218,25.437006,25.59929657,0.162290573,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,25.51312,25.14049381,24.95721,26.22237753,1.993,0.73663078,2.108783874,0.349315446,-1.5174
Alpha tubulin (Bter03085; XM_003396121),ATub,EQ8,EQ,Early,1,2,1,Control,17,25.1687851,25.25496483,25.20417023,25.20930672,25.1687851,25.25496483,0.086179733,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,25.22995,24.92400707,24.95721,26.22237753,1.993,0.86073561,2.448338325,0.351559095,-1.50816
Alpha tubulin (Bter03085; XM_003396121),ATub,EW3,EW,Early,2,1,1,Control,120,24.71879578,24.85164261,24.92774391,24.83272743,24.71879578,24.92774391,0.208948135,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,24.83273,24.92159537,24.95721,26.22237753,1.993,1.07084295,2.452413802,0.436648556,-1.19546
Alpha tubulin (Bter03085; XM_003396121),ATub,EW4,EW,Early,2,1,1,Control,120,25.00065613,24.83971024,24.84148026,24.89394887,24.83971024,25.00065613,0.160945892,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,24.89395,24.71233996,24.95721,26.22237753,1.993,1.0353953,2.833135032,0.365459212,-1.45222
Alpha tubulin (Bter03085; XM_003396121),ATub,EW5,EW,Early,2,1,1,Control,120,24.91464043,25.24183464,24.96658897,25.04102135,24.91464043,25.24183464,0.327194214,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.04102,25.76160569,24.95721,26.22237753,1.993,0.95496069,1.374056303,0.694993854,-0.52493
Alpha tubulin (Bter03085; XM_003396121),ATub,EW6,EW,Early,2,2,1,Control,17,23.97255325,23.97978401,23.96603966,23.97279231,23.96603966,23.97978401,0.013744354,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,23.99242,24.47244454,24.95721,26.22237753,1.993,1.6997671,3.342856272,0.508477471,-0.97574
Alpha tubulin (Bter03085; XM_003396121),ATub,EW7,EW,Early,2,2,1,Control,17,24.68728256,24.76081657,24.62542152,24.69117355,24.62542152,24.76081657,0.13539505,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,24.71139,24.71929254,24.95721,26.22237753,1.993,1.14472334,2.819583273,0.405990259,-1.30048
Alpha tubulin (Bter03085; XM_003396121),ATub,EW8,EW,Early,2,2,1,Control,17,24.45856285,24.42285728,24.3664341,24.41595141,24.3664341,24.45856285,0.092128754,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,24.43594,24.24981237,24.95721,26.22237753,1.993,1.33191653,3.897603549,0.341727041,-1.54908
Alpha tubulin (Bter03085; XM_003396121),ATub,MQ3,MQ,Medium,3,1,1,Control,120,24.74272728,23.87503815,23.95921135,24.19232559,23.87503815,24.74272728,0.867689133,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,24.19233,25.89835631,24.95721,26.22237753,1.993,1.52283489,1.25039354,1.217884484,0.28438
Alpha tubulin (Bter03085; XM_003396121),ATub,MQ4,MQ,Medium,3,1,1,Control,120,24.75932121,24.0994091,23.97704315,24.27859116,23.97704315,24.75932121,0.782278061,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,24.27859,25.52299076,24.95721,26.22237753,1.993,1.45228792,1.619837474,0.896563975,-0.15752
Alpha tubulin (Bter03085; XM_003396121),ATub,MQ5,MQ,Medium,3,1,1,Control,120,22.88893318,23.05644035,23.07144356,23.0056057,22.88893318,23.07144356,0.182510376,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,23.00561,25.12014958,24.95721,26.22237753,1.993,2.92442053,2.138579087,1.367459613,0.4515
Alpha tubulin (Bter03085; XM_003396121),ATub,MQ6,MQ,Medium,3,2,1,Control,17,23.27304459,23.1979332,23.31581497,23.26226425,23.1979332,23.31581497,0.117881775,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,23.28131,24.7511079,24.95721,26.22237753,1.993,2.51305177,2.75839205,0.911056775,-0.13439
Alpha tubulin (Bter03085; XM_003396121),ATub,MQ7,MQ,Medium,3,2,1,Control,17,23.051754,23.04985809,22.97179985,23.02447065,22.97179985,23.051754,0.079954147,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,23.04332,24.76237619,24.95721,26.22237753,1.993,2.86439366,2.737039456,1.046529912,0.06561
Alpha tubulin (Bter03085; XM_003396121),ATub,MQ8,MQ,Medium,3,2,1,Control,17,23.95144081,23.93533516,23.89963341,23.92880313,23.89963341,23.95144081,0.051807404,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,23.9484,24.74641872,24.95721,26.22237753,1.993,1.74141618,2.767326721,0.629277405,-0.66823
Alpha tubulin (Bter03085; XM_003396121),ATub,MW3,MW,Medium,4,1,1,Control,120,25.24067116,25.23079681,24.46321678,24.97822825,24.46321678,25.24067116,0.777454376,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,24.97823,25.29041589,24.95721,26.22237753,1.993,0.98850847,1.901644416,0.519817721,-0.94392
Alpha tubulin (Bter03085; XM_003396121),ATub,MW4,MW,Medium,4,1,1,Control,120,23.3206234,23.07178307,22.58325577,22.99188741,22.58325577,23.3206234,0.73736763,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,22.99189,24.63093495,24.95721,26.22237753,1.993,2.946563,2.996737273,0.983257034,-0.02436
Alpha tubulin (Bter03085; XM_003396121),ATub,MW5,MW,Medium,4,1,1,Control,120,23.97172356,23.98334312,23.84080887,23.93195852,23.84080887,23.98334312,0.142534256,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,23.93196,25.35247425,24.95721,26.22237753,1.993,1.75722795,1.821974874,0.964463328,-0.0522
Alpha tubulin (Bter03085; XM_003396121),ATub,MW6,MW,Medium,4,2,1,Control,17,22.84619904,22.73298454,22.89413261,22.82443873,22.73298454,22.89413261,0.161148071,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,22.84313,24.95817119,24.95721,26.22237753,1.993,3.19771057,2.391327326,1.337211569,0.41923
Alpha tubulin (Bter03085; XM_003396121),ATub,MW7,MW,Medium,4,2,1,Control,17,25.07860565,25.08334923,25.0057373,25.05589739,25.0057373,25.08334923,0.077611923,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,25.07641,25.70248086,24.95721,26.22237753,1.993,0.93655587,1.431241263,0.654366177,-0.61183
Alpha tubulin (Bter03085; XM_003396121),ATub,MW8,MW,Medium,4,2,1,Control,17,23.59167099,23.47355652,23.40968895,23.49163882,23.40968895,23.59167099,0.18198204,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,23.51087,24.59640071,24.95721,26.22237753,1.993,2.21504353,3.068964925,0.721755897,-0.47042
Alpha tubulin (Bter03085; XM_003396121),ATub,LQ3,LQ,Late,5,1,1,Control,120,24.93774796,25.00972176,24.92415428,24.957208,24.92415428,25.00972176,0.085567474,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,24.95721,26.22237753,24.95721,26.22237753,1.993,1,1,1,0
Alpha tubulin (Bter03085; XM_003396121),ATub,LQ4,LQ,Late,5,1,1,Control,120,25.13065147,25.19545174,25.55693245,25.29434522,25.13065147,25.55693245,0.426280975,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.29435,25.57874041,24.95721,26.22237753,1.993,0.83079164,1.558741203,0.532988821,-0.90782
Alpha tubulin (Bter03085; XM_003396121),ATub,LQ5,LQ,Late,5,1,1,Control,120,27.49603844,27.22458076,27.35603523,27.35888481,27.22458076,27.49603844,0.271457672,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,27.35888,26.77102083,24.95721,26.22237753,1.993,0.26698264,0.684979114,0.389767568,-1.35931
Alpha tubulin (Bter03085; XM_003396121),ATub,LQ6,LQ,Late,5,2,1,Control,17,24.63429451,24.79469299,24.94633675,24.79177475,24.63429451,24.94633675,0.312042236,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,24.81208,26.04640784,24.95721,26.22237753,1.993,1.08307229,1.129026689,0.959297335,-0.05995
Alpha tubulin (Bter03085; XM_003396121),ATub,LQ7,LQ,Late,5,2,1,Control,17,25.27015114,25.20082474,25.16306877,25.21134822,25.16306877,25.27015114,0.107082367,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,25.23199,24.94259361,24.95721,26.22237753,1.993,0.85976916,2.417155701,0.35569457,-1.49129
Alpha tubulin (Bter03085; XM_003396121),ATub,LQ8,LQ,Late,5,2,1,Control,17,25.3941288,25.4125843,25.34061623,25.38244311,25.34061623,25.4125843,0.071968079,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,25.40323,25.43112605,24.95721,26.22237753,1.993,0.78251233,1.725780648,0.453425141,-1.14106
Alpha tubulin (Bter03085; XM_003396121),ATub,LW3,LW,Late,6,1,1,Control,120,26.16055107,26.11499405,26.43378067,26.23644193,26.11499405,26.43378067,0.318786621,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,26.23644,26.28165682,24.95721,26.22237753,1.993,0.49490379,0.959942941,0.515555419,-0.9558
Alpha tubulin (Bter03085; XM_003396121),ATub,LW4,LW,Late,6,1,1,Control,120,25.50890541,26.10432816,25.4732933,25.69550896,25.4732933,26.10432816,0.631034851,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.69551,25.96422745,24.95721,26.22237753,1.993,0.66633832,1.194862228,0.557669586,-0.84252
Alpha tubulin (Bter03085; XM_003396121),ATub,LW5,LW,Late,6,1,1,Control,120,25.4530735,25.26184654,25.13415337,25.28302447,25.13415337,25.4530735,0.318920135,eq6,24.9181633,24.84104347,24.981884,24.91369692,24.84104347,24.981884,0.14084053,24.91369692,1.733,25.28302,25.98699537,24.95721,26.22237753,1.993,0.83597924,1.176247401,0.710717187,-0.49265
Alpha tubulin (Bter03085; XM_003396121),ATub,LW6,LW,Late,6,2,1,Control,17,25.11269951,24.83271027,24.91984367,24.95508448,24.83271027,25.11269951,0.279989243,eq6,24.62824821,24.69111824,24.67347717,24.66428121,24.62824821,24.69111824,0.062870026,24.91369692,1.733,25.20744,25.6146532,24.95721,26.22237753,1.993,0.87145454,1.520610107,0.573095322,-0.80315
Alpha tubulin (Bter03085; XM_003396121),ATub,LW7,LW,Late,6,2,1,Control,17,25.85536957,25.67505836,25.61236191,25.71426328,25.61236191,25.85536957,0.24300766,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,25.73532,25.49879946,24.95721,26.22237753,1.993,0.65191073,1.647088425,0.395795828,-1.33717
Alpha tubulin (Bter03085; XM_003396121),ATub,LW8,LW,Late,6,2,1,Control,17,26.36575508,26.28004456,26.40882301,26.35154088,26.28004456,26.40882301,0.128778458,eq6,24.9203167,24.93157196,24.82805061,24.89331309,24.82805061,24.93157196,0.103521347,24.91369692,1.733,26.37312,25.38152329,24.95721,26.22237753,1.993,0.45907383,1.78583768,0.257063583,-1.9598
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EQ3,EQ,Early,1,1,1,Control,130,26.41402435,26.20033264,26.41024208,26.34153303,26.20033264,26.41402435,0.213691711,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,26.34153,25.20347779,25.33085,26.22237753,1.993,0.56966872,2.019146866,0.282133377,-1.82555
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EQ4,EQ,Early,1,1,1,Control,130,26.91649055,27.25425148,26.92582703,27.03218969,26.91649055,27.25425148,0.337760925,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,27.03219,25.71323328,25.33085,26.22237753,1.993,0.38781513,1.420667419,0.272980943,-1.87313
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EQ5,EQ,Early,1,1,1,Control,130,26.11267281,26.23409081,26.19474411,26.18050257,26.11267281,26.23409081,0.121417999,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,26.1805,25.41256952,25.33085,26.22237753,1.993,0.62310151,1.748007975,0.356463768,-1.48817
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EQ6,EQ,Early,1,2,1,Control,18,24.48452759,24.47648048,24.52757454,24.4961942,24.47648048,24.52757454,0.051094055,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.7476,24.43307804,25.33085,26.22237753,1.993,1.38365719,3.434853817,0.402828552,-1.31176
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EQ7,EQ,Early,1,2,1,Control,18,26.68611717,26.82317162,26.60955048,26.70627975,26.60955048,26.82317162,0.21362114,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,26.98037,25.14049381,25.33085,26.22237753,1.993,0.39916779,2.108783874,0.189288147,-2.40134
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EQ8,EQ,Early,1,2,1,Control,18,25.16293716,25.29655647,25.11070824,25.19006729,25.11070824,25.29655647,0.185848236,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.44859,24.92400707,25.33085,26.22237753,1.993,0.93655068,2.448338325,0.382525027,-1.38637
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EW3,EW,Early,2,1,1,Control,130,25.45951843,25.23666,25.0717144,25.25596428,25.0717144,25.45951843,0.387804031,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.25596,24.92159537,25.33085,26.22237753,1.993,1.04257652,2.452413802,0.425122595,-1.23405
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EW4,EW,Early,2,1,1,Control,130,29.10202217,28.8743515,NA,28.98818684,28.8743515,29.10202217,0.22767067,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,28.98819,24.71233996,25.33085,26.22237753,1.993,0.13051895,2.833135032,0.046068737,-4.44007
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EW5,EW,Early,2,1,1,Control,130,26.78035736,26.59047318,26.62930489,26.66671181,26.59047318,26.78035736,0.189884186,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,26.66671,25.76160569,25.33085,26.22237753,1.993,0.47533045,1.374056303,0.345932291,-1.53144
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EW6,EW,Early,2,2,1,Control,18,24.12599754,24.05374336,24.07068634,24.08347575,24.05374336,24.12599754,0.072254181,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.33064,24.47244454,25.33085,26.22237753,1.993,1.74520368,3.342856272,0.522069613,-0.93769
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EW7,EW,Early,2,2,1,Control,18,24.95027542,25.09556961,25.06005096,25.03529867,24.95027542,25.09556961,0.145294189,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.29224,24.71929254,25.33085,26.22237753,1.993,1.02173365,2.819583273,0.362370447,-1.46446
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,EW8,EW,Early,2,2,1,Control,18,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.43258,24.24981237,25.33085,26.22237753,1.993,1.64891791,3.897603549,0.423059424,-1.24107
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MQ3,MQ,Medium,3,1,1,Control,130,25.59692955,25.52114868,25.54653549,25.55487124,25.52114868,25.59692955,0.075780869,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.55487,25.89835631,25.33085,26.22237753,1.993,0.88274162,1.25039354,0.70597103,-0.50232
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MQ4,MQ,Medium,3,1,1,Control,130,25.63899422,25.50071335,25.48693657,25.54221471,25.48693657,25.63899422,0.152057648,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.54221,25.52299076,25.33085,26.22237753,1.993,0.88898389,1.619837474,0.548810549,-0.86562
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MQ5,MQ,Medium,3,1,1,Control,130,25.27511406,25.68108177,25.200737,25.38564428,25.200737,25.68108177,0.480344772,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.38564,25.12014958,25.33085,26.22237753,1.993,0.96995572,2.138579087,0.453551484,-1.14066
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MQ6,MQ,Medium,3,2,1,Control,18,23.83075523,23.78974533,23.75878334,23.79309464,23.75878334,23.83075523,0.071971893,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.03728,24.7511079,25.33085,26.22237753,1.993,2.05484723,2.75839205,0.744943864,-0.4248
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MQ7,MQ,Medium,3,2,1,Control,18,23.89137077,23.852808,23.85343361,23.86587079,23.852808,23.89137077,0.038562775,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.11081,24.76237619,25.33085,26.22237753,1.993,1.97243182,2.737039456,0.720644277,-0.47264
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MQ8,MQ,Medium,3,2,1,Control,18,24.4739933,24.50804138,24.37812042,24.45338504,24.37812042,24.50804138,0.129920959,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.70435,24.74641872,25.33085,26.22237753,1.993,1.41737836,2.767326721,0.512183238,-0.96527
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MW3,MW,Medium,4,1,1,Control,130,25.71489143,26.28716469,25.62214279,25.87473297,25.62214279,26.28716469,0.665021896,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.87473,25.29041589,25.33085,26.22237753,1.993,0.7387415,1.901644416,0.388475097,-1.36411
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MW4,MW,Medium,4,1,1,Control,130,24.55984688,24.41773987,24.62398148,24.53385607,24.41773987,24.62398148,0.206241608,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,24.53386,24.63093495,25.33085,26.22237753,1.993,1.55851291,2.996737273,0.520069919,-0.94322
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MW5,MW,Medium,4,1,1,Control,130,25.37984657,25.39986801,25.45455933,25.41142464,25.37984657,25.45455933,0.074712753,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.41142,25.35247425,25.33085,26.22237753,1.993,0.95613306,1.821974874,0.5247784,-0.93022
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MW6,MW,Medium,4,2,1,Control,18,24.75507355,24.42736053,24.5582943,24.58024279,24.42736053,24.75507355,0.327713013,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.83251,24.95817119,25.33085,26.22237753,1.993,1.31976732,2.391327326,0.551897393,-0.85753
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MW7,MW,Medium,4,2,1,Control,18,26.1132946,26.09984207,25.95087814,26.05467161,25.95087814,26.1132946,0.162416458,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,26.32207,25.70248086,25.33085,26.22237753,1.993,0.57587518,1.431241263,0.402360659,-1.31344
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,MW8,MW,Medium,4,2,1,Control,18,25.44401169,25.12947273,25.08096313,25.21814919,25.08096313,25.44401169,0.363048553,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.47696,24.59640071,25.33085,26.22237753,1.993,0.92187391,3.068964925,0.30038594,-1.73511
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LQ3,LQ,Late,5,1,1,Control,130,25.45887184,25.21414757,25.31954193,25.33085378,25.21414757,25.45887184,0.244724274,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.33085,26.22237753,25.33085,26.22237753,1.993,1,1,1,0
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LQ4,LQ,Late,5,1,1,Control,130,26.02511406,25.65532112,25.61689568,25.76577695,25.61689568,26.02511406,0.408218384,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,25.76578,25.57874041,25.33085,26.22237753,1.993,0.78494199,1.558741203,0.503574289,-0.98972
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LQ5,LQ,Late,5,1,1,Control,130,28.46827316,28.29141998,28.56369019,28.44112778,28.29141998,28.56369019,0.272270203,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,28.44113,26.77102083,25.33085,26.22237753,1.993,0.1769905,0.684979114,0.258388171,-1.95239
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LQ6,LQ,Late,5,2,1,Control,18,24.75657272,24.76680183,25.08817291,24.87051582,24.75657272,25.08817291,0.331600189,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.12576,26.04640784,25.33085,26.22237753,1.993,1.12096061,1.129026689,0.992855723,-0.01034
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LQ7,LQ,Late,5,2,1,Control,18,25.45812798,25.2067337,25.33536148,25.33340772,25.2067337,25.45812798,0.251394272,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.5934,24.94259361,25.33085,26.22237753,1.993,0.86400552,2.417155701,0.357447195,-1.4842
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LQ8,LQ,Late,5,2,1,Control,18,25.8076725,25.78689384,25.61870384,25.73775673,25.61870384,25.8076725,0.188968658,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,26.0019,25.43112605,25.33085,26.22237753,1.993,0.68824548,1.725780648,0.398802408,-1.32625
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LW3,LW,Late,6,1,1,Control,130,28.69249153,28.27706146,28.29927826,28.42294375,28.27706146,28.69249153,0.415430069,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,28.42294,26.28165682,25.33085,26.22237753,1.993,0.17879146,0.959942941,0.186252175,-2.42467
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LW4,LW,Late,6,1,1,Control,130,27.11413193,26.68571663,27.12473488,26.97486115,26.68571663,27.12473488,0.43901825,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,26.97486,25.96422745,25.33085,26.22237753,1.993,0.40039305,1.194862228,0.335095579,-1.57736
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LW5,LW,Late,6,1,1,Control,130,27.99466705,27.79139519,28.41835594,28.06813939,27.79139519,28.41835594,0.626960754,ew8,24.50312614,24.2505722,24.54403496,24.43257777,24.2505722,24.54403496,0.293462753,24.43257777,1.745,28.06814,25.98699537,25.33085,26.22237753,1.993,0.21783961,1.176247401,0.1851988,-2.43285
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LW6,LW,Late,6,2,1,Control,18,24.1865406,24.11751366,24.15719223,24.15374883,24.11751366,24.1865406,0.069026947,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,24.40164,25.6146532,25.33085,26.22237753,1.993,1.67756758,1.520610107,1.103220063,0.14172
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LW7,LW,Late,6,2,1,Control,18,25.12590218,24.96860123,24.93891907,25.01114082,24.93891907,25.12590218,0.186983109,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.26783,25.49879946,25.33085,26.22237753,1.993,1.03571175,1.647088425,0.628813688,-0.6693
Arginine Kinase (Bter08073; XM_012317366; XM_003401454),AK,LW8,LW,Late,6,2,1,Control,18,25.67605782,25.60389519,25.6529007,25.64428457,25.60389519,25.67605782,0.072162628,ew8,24.27244186,24.10513878,24.17554092,24.18437386,24.10513878,24.27244186,0.167303085,24.43257777,1.745,25.90747,25.38152329,25.33085,26.22237753,1.993,0.72539818,1.78583768,0.406194915,-1.29976
GAPDH (Bter05491; XM_003398087),GAPDH,EQ3,EQ,Early,1,1,1,Control,129,24.91206741,24.43627739,24.44542885,24.59792455,24.43627739,24.91206741,0.475790024,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.59792,25.20347779,24.23425,26.22237753,1.993,0.78146426,2.019146866,0.387026951,-1.36949
GAPDH (Bter05491; XM_003398087),GAPDH,EQ4,EQ,Early,1,1,1,Control,129,24.75634766,24.84217262,24.87073135,24.82308388,24.75634766,24.87073135,0.114383698,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.82308,25.71323328,24.23425,26.22237753,1.993,0.67082201,1.420667419,0.472187929,-1.08257
GAPDH (Bter05491; XM_003398087),GAPDH,EQ5,EQ,Early,1,1,1,Control,129,24.06806564,24.08815765,24.10023308,24.08548546,24.06806564,24.10023308,0.032167435,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.08549,25.41256952,24.23425,26.22237753,1.993,1.10612733,1.748007975,0.632793069,-0.66019
GAPDH (Bter05491; XM_003398087),GAPDH,EQ6,EQ,Early,1,2,1,Control,24,22.49689102,22.56336212,22.6314621,22.56390508,22.49689102,22.6314621,0.134571075,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.28588,24.43307804,24.23425,26.22237753,1.993,0.96559329,3.434853817,0.28111627,-1.83076
GAPDH (Bter05491; XM_003398087),GAPDH,EQ7,EQ,Early,1,2,1,Control,24,23.91828728,23.93341446,23.99798965,23.94989713,23.91828728,23.99798965,0.079702377,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.77765,25.14049381,24.23425,26.22237753,1.993,0.35117196,2.108783874,0.166528191,-2.58616
GAPDH (Bter05491; XM_003398087),GAPDH,EQ8,EQ,Early,1,2,1,Control,24,23.7656517,23.84302139,23.68109894,23.76325734,23.68109894,23.84302139,0.161922455,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.57677,24.92400707,24.23425,26.22237753,1.993,0.4024141,2.448338325,0.16436213,-2.60505
GAPDH (Bter05491; XM_003398087),GAPDH,EW3,EW,Early,2,1,1,Control,129,23.21421242,23.1210289,23.01558113,23.11694082,23.01558113,23.21421242,0.198631287,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,23.11694,24.92159537,24.23425,26.22237753,1.993,2.13308818,2.452413802,0.869791297,-0.20126
GAPDH (Bter05491; XM_003398087),GAPDH,EW4,EW,Early,2,1,1,Control,129,27.38426399,27.48329163,27.32852554,27.39869372,27.32852554,27.48329163,0.154766083,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,27.39869,24.71233996,24.23425,26.22237753,1.993,0.11699766,2.833135032,0.041296183,-4.59785
GAPDH (Bter05491; XM_003398087),GAPDH,EW5,EW,Early,2,1,1,Control,129,24.91310883,24.89247322,24.93061829,24.91206678,24.89247322,24.93061829,0.038145065,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.91207,25.76160569,24.23425,26.22237753,1.993,0.63154579,1.374056303,0.459621476,-1.12148
GAPDH (Bter05491; XM_003398087),GAPDH,EW6,EW,Early,2,2,1,Control,24,22.43984222,22.50960922,22.48015404,22.47653516,22.43984222,22.50960922,0.069766998,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.19185,24.47244454,24.23425,26.22237753,1.993,1.02916531,3.342856272,0.307870046,-1.69961
GAPDH (Bter05491; XM_003398087),GAPDH,EW7,EW,Early,2,2,1,Control,24,22.99027061,23.1172657,23.08207512,23.06320381,22.99027061,23.1172657,0.126995087,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.82329,24.71929254,24.23425,26.22237753,1.993,0.67072916,2.819583273,0.237882374,-2.07168
GAPDH (Bter05491; XM_003398087),GAPDH,EW8,EW,Early,2,2,1,Control,24,22.46686935,22.30665207,22.26391411,22.34581184,22.26391411,22.46686935,0.202955246,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.05115,24.24981237,24.23425,26.22237753,1.993,1.13218242,3.897603549,0.290481679,-1.78348
GAPDH (Bter05491; XM_003398087),GAPDH,MQ3,MQ,Medium,3,1,1,Control,129,23.84801292,23.68585014,23.6594696,23.73111089,23.6594696,23.84801292,0.18854332,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,23.73111,25.89835631,24.23425,26.22237753,1.993,1.40655418,1.25039354,1.124889191,0.16978
GAPDH (Bter05491; XM_003398087),GAPDH,MQ4,MQ,Medium,3,1,1,Control,129,23.39458466,23.57845116,23.54954529,23.50752703,23.39458466,23.57845116,0.183866501,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,23.50753,25.52299076,24.23425,26.22237753,1.993,1.63679526,1.619837474,1.010468822,0.01502
GAPDH (Bter05491; XM_003398087),GAPDH,MQ5,MQ,Medium,3,1,1,Control,129,23.31988144,23.01498985,23.14555931,23.16014353,23.01498985,23.31988144,0.304891586,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,23.16014,25.12014958,24.23425,26.22237753,1.993,2.07151016,2.138579087,0.968638557,-0.04597
GAPDH (Bter05491; XM_003398087),GAPDH,MQ6,MQ,Medium,3,2,1,Control,24,21.85686493,21.903862,21.70406914,21.82159869,21.70406914,21.903862,0.199792862,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,23.48693,24.7511079,24.23425,26.22237753,1.993,1.65981545,2.75839205,0.601732972,-0.7328
GAPDH (Bter05491; XM_003398087),GAPDH,MQ7,MQ,Medium,3,2,1,Control,24,21.67775345,21.62182617,21.75875854,21.68611272,21.62182617,21.75875854,0.136932373,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,23.3411,24.76237619,24.23425,26.22237753,1.993,1.83231669,2.737039456,0.669452054,-0.57895
GAPDH (Bter05491; XM_003398087),GAPDH,MQ8,MQ,Medium,3,2,1,Control,24,22.55595207,22.42200089,22.4317627,22.46990522,22.42200089,22.55595207,0.133951187,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.18471,24.74641872,24.23425,26.22237753,1.993,1.03415688,2.767326721,0.373702488,-1.42004
GAPDH (Bter05491; XM_003398087),GAPDH,MW3,MW,Medium,4,1,1,Control,129,23.28741264,23.32790375,23.34094238,23.31875292,23.28741264,23.34094238,0.053529739,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,23.31875,25.29041589,24.23425,26.22237753,1.993,1.8602956,1.901644416,0.978256285,-0.03172
GAPDH (Bter05491; XM_003398087),GAPDH,MW4,MW,Medium,4,1,1,Control,129,22.71736908,22.69108772,22.60532761,22.67126147,22.60532761,22.71736908,0.112041473,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,22.67126,24.63093495,24.23425,26.22237753,1.993,2.88566774,2.996737273,0.962936514,-0.05449
GAPDH (Bter05491; XM_003398087),GAPDH,MW5,MW,Medium,4,1,1,Control,129,24.27483368,24.34562111,24.31109047,24.31051509,24.27483368,24.34562111,0.07078743,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.31052,25.35247425,24.23425,26.22237753,1.993,0.94960173,1.821974874,0.521193649,-0.94011
GAPDH (Bter05491; XM_003398087),GAPDH,MW6,MW,Medium,4,2,1,Control,24,22.25027275,22.33125496,22.32084846,22.30079206,22.25027275,22.33125496,0.080982208,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.00269,24.95817119,24.23425,26.22237753,1.993,1.16999746,2.391327326,0.489266962,-1.03131
GAPDH (Bter05491; XM_003398087),GAPDH,MW7,MW,Medium,4,2,1,Control,24,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.91395,25.70248086,24.23425,26.22237753,1.993,0.32017295,1.431241263,0.223702992,-2.16034
GAPDH (Bter05491; XM_003398087),GAPDH,MW8,MW,Medium,4,2,1,Control,24,22.88829803,22.82891083,22.90127563,22.87282817,22.82891083,22.90127563,0.072364807,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,24.61838,24.59640071,24.23425,26.22237753,1.993,0.77069862,3.068964925,0.251126567,-1.99351
GAPDH (Bter05491; XM_003398087),GAPDH,LQ3,LQ,Late,5,1,1,Control,129,24.19168854,24.23388481,24.27716637,24.23424657,24.19168854,24.27716637,0.085477829,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.23425,26.22237753,24.23425,26.22237753,1.993,1,1,1,0
GAPDH (Bter05491; XM_003398087),GAPDH,LQ4,LQ,Late,5,1,1,Control,129,24.63498306,24.5620842,24.596035,24.59770075,24.5620842,24.63498306,0.072898865,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.5977,25.57874041,24.23425,26.22237753,1.993,0.78158284,1.558741203,0.50141925,-0.99591
GAPDH (Bter05491; XM_003398087),GAPDH,LQ5,LQ,Late,5,1,1,Control,129,26.83956528,26.79712105,26.69919777,26.77862803,26.69919777,26.83956528,0.140367508,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,26.77863,26.77102083,24.23425,26.22237753,1.993,0.1781417,0.684979114,0.260068806,-1.94303
GAPDH (Bter05491; XM_003398087),GAPDH,LQ6,LQ,Late,5,2,1,Control,24,23.60577202,23.4861145,23.59090424,23.56093025,23.4861145,23.60577202,0.119657516,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.359,26.04640784,24.23425,26.22237753,1.993,0.46644287,1.129026689,0.413137152,-1.27531
GAPDH (Bter05491; XM_003398087),GAPDH,LQ7,LQ,Late,5,2,1,Control,24,23.6021328,23.6303463,23.71878815,23.65042241,23.6021328,23.71878815,0.11665535,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.45532,24.94259361,24.23425,26.22237753,1.993,0.43695321,2.417155701,0.180771645,-2.46776
GAPDH (Bter05491; XM_003398087),GAPDH,LQ8,LQ,Late,5,2,1,Control,24,24.34797859,23.96710587,23.96765327,24.09424591,23.96710587,24.34797859,0.380872726,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.93301,25.43112605,24.23425,26.22237753,1.993,0.31606043,1.725780648,0.18314056,-2.44898
GAPDH (Bter05491; XM_003398087),GAPDH,LW3,LW,Late,6,1,1,Control,129,26.82068634,26.56377411,26.48084831,26.62176959,26.48084831,26.82068634,0.339838028,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,26.62177,26.28165682,24.23425,26.22237753,1.993,0.19813222,0.959942941,0.206399998,-2.27649
GAPDH (Bter05491; XM_003398087),GAPDH,LW4,LW,Late,6,1,1,Control,129,24.81895828,24.87557983,24.76982117,24.82145309,24.76982117,24.87557983,0.105758667,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,24.82145,25.96422745,24.23425,26.22237753,1.993,0.67156416,1.194862228,0.562043176,-0.83125
GAPDH (Bter05491; XM_003398087),GAPDH,LW5,LW,Late,6,1,1,Control,129,26.90208054,26.93702507,27.03534698,26.95815086,26.90208054,27.03534698,0.133266449,MW7,25.88874245,25.91699982,25.93610191,25.91394806,25.88874245,25.93610191,0.047359467,25.91394806,1.97,26.95815,25.98699537,24.23425,26.22237753,1.993,0.1577256,1.176247401,0.134092197,-2.8987
GAPDH (Bter05491; XM_003398087),GAPDH,LW6,LW,Late,6,2,1,Control,24,23.79363251,23.87203598,23.88958549,23.85175133,23.79363251,23.88958549,0.095952988,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,25.67201,25.6146532,24.23425,26.22237753,1.993,0.3772472,1.520610107,0.24808937,-2.01107
GAPDH (Bter05491; XM_003398087),GAPDH,LW7,LW,Late,6,2,1,Control,24,24.21774292,24.14757156,24.20331001,24.1895415,24.14757156,24.21774292,0.070171356,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,26.03558,25.49879946,24.23425,26.22237753,1.993,0.29482702,1.647088425,0.178998901,-2.48198
GAPDH (Bter05491; XM_003398087),GAPDH,LW8,LW,Late,6,2,1,Control,24,24.57852554,24.44684601,24.74696922,24.59078026,24.44684601,24.74696922,0.300123215,MW7,24.04982376,24.13998985,24.03977966,24.07653109,24.03977966,24.13998985,0.10021019,25.91394806,1.97,26.46744,25.38152329,24.23425,26.22237753,1.993,0.21998812,1.78583768,0.123184833,-3.0211
Peptidylprolyl isomerase A,PPIA,EQ3,EQ,Early,1,1,1,Control,119,24.87050629,24.83188438,24.99140358,24.89793142,24.83188438,24.99140358,0.159519196,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.89793,25.20347779,24.99129,26.22237753,1.993,1.06017733,2.019146866,0.525062019,-0.92944
Peptidylprolyl isomerase A,PPIA,EQ4,EQ,Early,1,1,1,Control,119,24.91679001,25.00720024,24.93945694,24.9544824,24.91679001,25.00720024,0.090410233,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.95448,25.71323328,24.99129,26.22237753,1.993,1.0233062,1.420667419,0.720299616,-0.47333
Peptidylprolyl isomerase A,PPIA,EQ5,EQ,Early,1,1,1,Control,119,24.94062805,24.94618988,24.86946106,24.91875966,24.86946106,24.94618988,0.076728821,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.91876,25.41256952,24.99129,26.22237753,1.993,1.04644529,1.748007975,0.598650183,-0.74021
Peptidylprolyl isomerase A,PPIA,EQ6,EQ,Early,1,2,1,Control,16,23.51783562,23.76676369,23.51238823,23.59899584,23.51238823,23.76676369,0.254375458,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,23.94432,24.43307804,24.99129,26.22237753,1.993,1.92578869,3.434853817,0.560661034,-0.8348
Peptidylprolyl isomerase A,PPIA,EQ7,EQ,Early,1,2,1,Control,16,25.19945145,24.93836594,25.03079796,25.05620511,24.93836594,25.19945145,0.26108551,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.42286,25.14049381,24.99129,26.22237753,1.993,0.76327648,2.108783874,0.361951023,-1.46613
Peptidylprolyl isomerase A,PPIA,EQ8,EQ,Early,1,2,1,Control,16,24.52829933,24.68946457,24.63518143,24.61764844,24.52829933,24.68946457,0.161165237,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,24.97788,24.92400707,24.99129,26.22237753,1.993,1.00842677,2.448338325,0.411882118,-1.2797
Peptidylprolyl isomerase A,PPIA,EW3,EW,Early,2,1,1,Control,119,24.20983124,24.03105545,23.98186874,24.07425181,23.98186874,24.20983124,0.227962494,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.07425,24.92159537,24.99129,26.22237753,1.993,1.7753701,2.452413802,0.723927624,-0.46608
Peptidylprolyl isomerase A,PPIA,EW4,EW,Early,2,1,1,Control,119,23.73291206,23.84271622,23.64976883,23.74179904,23.64976883,23.84271622,0.192947388,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,23.7418,24.71233996,24.99129,26.22237753,1.993,2.18606664,2.833135032,0.771606935,-0.37406
Peptidylprolyl isomerase A,PPIA,EW5,EW,Early,2,1,1,Control,119,25.1430645,24.91292381,24.9096508,24.98854637,24.9096508,25.1430645,0.233413696,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.98855,25.76160569,24.99129,26.22237753,1.993,1.00171828,1.374056303,0.72902273,-0.45596
Peptidylprolyl isomerase A,PPIA,EW6,EW,Early,2,2,1,Control,16,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,23.98309,24.47244454,24.99129,26.22237753,1.993,1.87961647,3.342856272,0.562278577,-0.83064
Peptidylprolyl isomerase A,PPIA,EW7,EW,Early,2,2,1,Control,16,24.05797386,24.20338631,24.17971802,24.14702606,24.05797386,24.20338631,0.145412445,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,24.50037,24.71929254,24.99129,26.22237753,1.993,1.35972538,2.819583273,0.482243384,-1.05217
Peptidylprolyl isomerase A,PPIA,EW8,EW,Early,2,2,1,Control,16,23.95874786,23.98536873,24.03924561,23.99445407,23.95874786,24.03924561,0.080497742,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,24.34557,24.24981237,24.99129,26.22237753,1.993,1.49807494,3.897603549,0.384357957,-1.37948
Peptidylprolyl isomerase A,PPIA,MQ3,MQ,Medium,3,1,1,Control,119,24.45910835,24.30418968,24.37823296,24.38051033,24.30418968,24.45910835,0.154918671,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.38051,25.89835631,24.99129,26.22237753,1.993,1.46566627,1.25039354,1.172163979,0.22917
Peptidylprolyl isomerase A,PPIA,MQ4,MQ,Medium,3,1,1,Control,119,24.56461525,24.59337807,24.50148582,24.55315971,24.50148582,24.59337807,0.091892242,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.55316,25.52299076,24.99129,26.22237753,1.993,1.31553318,1.619837474,0.812138997,-0.3002
Peptidylprolyl isomerase A,PPIA,MQ5,MQ,Medium,3,1,1,Control,119,23.76585007,23.92354584,23.73823929,23.80921173,23.73823929,23.92354584,0.185306549,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,23.80921,25.12014958,24.99129,26.22237753,1.993,2.09574202,2.138579087,0.97996938,-0.02919
Peptidylprolyl isomerase A,PPIA,MQ6,MQ,Medium,3,2,1,Control,16,23.52456474,23.31193542,23.47681236,23.43777084,23.31193542,23.52456474,0.212629318,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,23.78074,24.7511079,24.99129,26.22237753,1.993,2.13342625,2.75839205,0.773431118,-0.37066
Peptidylprolyl isomerase A,PPIA,MQ7,MQ,Medium,3,2,1,Control,16,23.8692627,23.95159531,23.7193203,23.8467261,23.7193203,23.95159531,0.232275009,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,24.19568,24.76237619,24.99129,26.22237753,1.993,1.64543145,2.737039456,0.601171987,-0.73415
Peptidylprolyl isomerase A,PPIA,MQ8,MQ,Medium,3,2,1,Control,16,24.08543587,24.15802002,24.24641418,24.16329002,24.08543587,24.24641418,0.160978317,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,24.51688,24.74641872,24.99129,26.22237753,1.993,1.34575278,2.767326721,0.486300649,-1.04008
Peptidylprolyl isomerase A,PPIA,MW3,MW,Medium,4,1,1,Control,119,23.96127892,23.94730949,23.78920364,23.89926402,23.78920364,23.96127892,0.172075272,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,23.89926,25.29041589,24.99129,26.22237753,1.993,1.98087861,1.901644416,1.041666148,0.05889
Peptidylprolyl isomerase A,PPIA,MW4,MW,Medium,4,1,1,Control,119,23.59807777,23.48451996,23.52219963,23.53493245,23.48451996,23.59807777,0.113557816,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,23.53493,24.63093495,24.99129,26.22237753,1.993,2.4882748,2.996737273,0.830327977,-0.26825
Peptidylprolyl isomerase A,PPIA,MW5,MW,Medium,4,1,1,Control,119,24.45255089,24.36309433,24.37566566,24.39710363,24.36309433,24.45255089,0.089456558,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.3971,25.35247425,24.99129,26.22237753,1.993,1.45052208,1.821974874,0.79612628,-0.32893
Peptidylprolyl isomerase A,PPIA,MW6,MW,Medium,4,2,1,Control,16,23.52868843,23.23213005,23.21191788,23.32424545,23.21191788,23.52868843,0.316770554,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,23.66555,24.95817119,24.99129,26.22237753,1.993,2.2929266,2.391327326,0.958850999,-0.06062
Peptidylprolyl isomerase A,PPIA,MW7,MW,Medium,4,2,1,Control,16,25.02876282,25.05074501,24.82114029,24.96688271,24.82114029,25.05074501,0.229604721,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.33223,25.70248086,24.99129,26.22237753,1.993,0.8078277,1.431241263,0.564424549,-0.82515
Peptidylprolyl isomerase A,PPIA,MW8,MW,Medium,4,2,1,Control,16,23.69641113,23.74757576,23.78208542,23.7420241,23.69641113,23.78208542,0.085674286,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,24.08945,24.59640071,24.99129,26.22237753,1.993,1.7585659,3.068964925,0.573015965,-0.80335
Peptidylprolyl isomerase A,PPIA,LQ3,LQ,Late,5,1,1,Control,119,24.9840641,24.99740791,24.9923954,24.99128914,24.9840641,24.99740791,0.013343811,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.99129,26.22237753,24.99129,26.22237753,1.993,1,1,1,0
Peptidylprolyl isomerase A,PPIA,LQ4,LQ,Late,5,1,1,Control,119,25.07180786,24.99045944,24.92518044,24.99581591,24.92518044,25.07180786,0.146627426,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,24.99582,25.57874041,24.99129,26.22237753,1.993,0.99717053,1.558741203,0.639728087,-0.64447
Peptidylprolyl isomerase A,PPIA,LQ5,LQ,Late,5,1,1,Control,119,26.15399742,26.3902092,26.43406487,26.32609049,26.15399742,26.43406487,0.280067444,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,26.32609,26.77102083,24.99129,26.22237753,1.993,0.4336561,0.684979114,0.633093895,-0.65951
Peptidylprolyl isomerase A,PPIA,LQ6,LQ,Late,5,2,1,Control,16,24.95391655,24.96801567,25.09404182,25.00532468,24.95391655,25.09404182,0.140125275,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.37123,26.04640784,24.99129,26.22237753,1.993,0.78834387,1.129026689,0.698250874,-0.51818
Peptidylprolyl isomerase A,PPIA,LQ7,LQ,Late,5,2,1,Control,16,24.80206871,24.67831039,24.58798599,24.68945503,24.58798599,24.80206871,0.214082718,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.05074,24.94259361,24.99129,26.22237753,1.993,0.96347114,2.417155701,0.398597053,-1.327
Peptidylprolyl isomerase A,PPIA,LQ8,LQ,Late,5,2,1,Control,16,24.96518707,25.22072601,24.95692825,25.04761378,24.95692825,25.22072601,0.26379776,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.41414,25.43112605,24.99129,26.22237753,1.993,0.76745256,1.725780648,0.444698785,-1.1691
Peptidylprolyl isomerase A,PPIA,LW3,LW,Late,6,1,1,Control,119,26.18377686,26.60946846,26.30664825,26.36663119,26.18377686,26.60946846,0.425691605,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,26.36663,26.28165682,24.99129,26.22237753,1.993,0.4227901,0.959942941,0.44043253,-1.18301
Peptidylprolyl isomerase A,PPIA,LW4,LW,Late,6,1,1,Control,119,25.18603325,24.97498512,25.00643539,25.05581792,24.97498512,25.18603325,0.211048126,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,25.05582,25.96422745,24.99129,26.22237753,1.993,0.9604138,1.194862228,0.803786226,-0.31512
Peptidylprolyl isomerase A,PPIA,LW5,LW,Late,6,1,1,Control,119,25.94371796,25.81283188,25.86149597,25.87268194,25.81283188,25.94371796,0.130886078,ew6,23.91893387,23.94215584,24.08819389,23.98309453,23.91893387,24.08819389,0.169260025,23.98309453,1.87,25.87268,25.98699537,24.99129,26.22237753,1.993,0.5759712,1.176247401,0.489668412,-1.03012
Peptidylprolyl isomerase A,PPIA,LW6,LW,Late,6,2,1,Control,16,24.44255829,24.77818871,24.48774147,24.56949615,24.44255829,24.77818871,0.335630417,ew6,23.24247169,23.9883976,23.37168694,23.53418541,23.24247169,23.9883976,0.745925903,23.98309453,1.87,25.03815,25.6146532,24.99129,26.22237753,1.993,0.97109181,1.520610107,0.638619859,-0.64697
Peptidylprolyl isomerase A,PPIA,LW7,LW,Late,6,2,1,Control,16,25.30870628,25.15282059,25.05127144,25.17093277,25.05127144,25.30870628,0.257434845,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.53926,25.49879946,24.99129,26.22237753,1.993,0.70963946,1.647088425,0.430844783,-1.21476
Peptidylprolyl isomerase A,PPIA,LW8,LW,Late,6,2,1,Control,16,25.46641541,25.37747002,25.5359745,25.45995331,25.37747002,25.5359745,0.158504486,ew6,23.59986496,23.74976349,23.56199265,23.63720703,23.56199265,23.74976349,0.187770844,23.98309453,1.87,25.83251,25.38152329,24.99129,26.22237753,1.993,0.59063646,1.78583768,0.330733562,-1.59626
Phospholipase A2 (Bter07778; XM_003400908),PA2,EQ3,EQ,Early,1,1,1,Control,121,29.5983181,29.71265221,29.86693001,29.72596677,29.5983181,29.86693001,0.268611908,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.17342,25.20347779,31.87633,26.22237753,1.993,1.59388489,2.019146866,0.789385319,-0.3412
Phospholipase A2 (Bter07778; XM_003400908),PA2,EQ4,EQ,Early,1,1,1,Control,121,29.58731461,30.19406128,30.20191956,29.99443181,29.58731461,30.20191956,0.61460495,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.45496,25.71323328,31.87633,26.22237753,1.993,1.32241259,1.420667419,0.930838969,-0.1034
Phospholipase A2 (Bter07778; XM_003400908),PA2,EQ5,EQ,Early,1,1,1,Control,121,29.9835453,29.91746712,29.80116844,29.90072695,29.80116844,29.9835453,0.182376862,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.35669,25.41256952,31.87633,26.22237753,1.993,1.41146648,1.748007975,0.807471419,-0.30852
Phospholipase A2 (Bter07778; XM_003400908),PA2,EQ6,EQ,Early,1,2,1,Control,65,30.76279259,30.85488319,30.61744881,30.74504153,30.61744881,30.85488319,0.237434387,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,30.74504,24.43307804,31.87633,26.22237753,1.993,2.1175835,3.434853817,0.61649887,-0.69783
Phospholipase A2 (Bter07778; XM_003400908),PA2,EQ7,EQ,Early,1,2,1,Control,65,31.91845322,31.89804649,31.41417503,31.74355825,31.41417503,31.91845322,0.504278183,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.74356,25.14049381,31.87633,26.22237753,1.993,1.09204927,2.108783874,0.517857368,-0.94937
Phospholipase A2 (Bter07778; XM_003400908),PA2,EQ8,EQ,Early,1,2,1,Control,65,31.8838253,31.83666039,31.37043953,31.69697507,31.37043953,31.8838253,0.513385773,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.69698,24.92400707,31.87633,26.22237753,1.993,1.12631372,2.448338325,0.460031896,-1.12019
Phospholipase A2 (Bter07778; XM_003400908),PA2,EW3,EW,Early,2,1,1,Control,121,29.80817223,29.72428322,29.44768906,29.66004817,29.44768906,29.80817223,0.36048317,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.10429,24.92159537,31.87633,26.22237753,1.993,1.66865939,2.452413802,0.680415103,-0.55551
Phospholipase A2 (Bter07778; XM_003400908),PA2,EW4,EW,Early,2,1,1,Control,121,30.07168388,29.45988655,29.72184944,29.75113996,29.45988655,30.07168388,0.611797333,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.19982,24.71233996,31.87633,26.22237753,1.993,1.56622223,2.833135032,0.552823007,-0.85511
Phospholipase A2 (Bter07778; XM_003400908),PA2,EW5,EW,Early,2,1,1,Control,121,29.68628693,29.75805283,29.85533142,29.76655706,29.68628693,29.85533142,0.169044495,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.21599,25.76160569,31.87633,26.22237753,1.993,1.54951806,1.374056303,1.12769619,0.17338
Phospholipase A2 (Bter07778; XM_003400908),PA2,EW6,EW,Early,2,2,1,Control,65,31.33247757,31.25357246,31.48913574,31.35839526,31.25357246,31.48913574,0.235563278,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.3584,24.47244454,31.87633,26.22237753,1.993,1.4098698,3.342856272,0.421756033,-1.24552
Phospholipase A2 (Bter07778; XM_003400908),PA2,EW7,EW,Early,2,2,1,Control,65,30.63671112,30.74251175,30.79714203,30.72545497,30.63671112,30.79714203,0.160430908,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,30.72545,24.71929254,31.87633,26.22237753,1.993,2.14527007,2.819583273,0.760846503,-0.39432
Phospholipase A2 (Bter07778; XM_003400908),PA2,EW8,EW,Early,2,2,1,Control,65,30.91497803,30.81124687,31.62770271,31.11797587,30.81124687,31.62770271,0.816455841,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.11798,24.24981237,31.87633,26.22237753,1.993,1.65358274,3.897603549,0.424256269,-1.23699
Phospholipase A2 (Bter07778; XM_003400908),PA2,MQ3,MQ,Medium,3,1,1,Control,121,30.74809074,29.86618233,30.12872696,30.24766668,29.86618233,30.74809074,0.881908417,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.72052,25.89835631,31.87633,26.22237753,1.993,1.10886144,1.25039354,0.886809955,-0.1733
Phospholipase A2 (Bter07778; XM_003400908),PA2,MQ4,MQ,Medium,3,1,1,Control,121,29.6409893,29.07511902,29.71085548,29.4756546,29.07511902,29.71085548,0.635736465,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,30.91092,25.52299076,31.87633,26.22237753,1.993,1.89698479,1.619837474,1.171095757,0.22786
Phospholipase A2 (Bter07778; XM_003400908),PA2,MQ5,MQ,Medium,3,1,1,Control,121,29.6576767,29.4767952,29.8968029,29.6770916,29.4767952,29.8968029,0.420007706,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.12216,25.12014958,31.87633,26.22237753,1.993,1.64899645,2.138579087,0.771071064,-0.37506
Phospholipase A2 (Bter07778; XM_003400908),PA2,MQ6,MQ,Medium,3,2,1,Control,65,30.53046608,30.68361282,30.55613708,30.590072,30.53046608,30.68361282,0.153146744,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,30.59007,24.7511079,31.87633,26.22237753,1.993,2.34679809,2.75839205,0.850784823,-0.23313
Phospholipase A2 (Bter07778; XM_003400908),PA2,MQ7,MQ,Medium,3,2,1,Control,65,31.89284325,32.08221054,31.70589256,31.89364878,31.70589256,32.08221054,0.376317978,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.89365,24.76237619,31.87633,26.22237753,1.993,0.98858115,2.737039456,0.361186298,-1.46918
Phospholipase A2 (Bter07778; XM_003400908),PA2,MQ8,MQ,Medium,3,2,1,Control,65,31.61898804,32.02055359,32.05805588,31.89919917,31.61898804,32.05805588,0.439067841,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.8992,24.74641872,31.87633,26.22237753,1.993,0.98494884,2.767326721,0.355920691,-1.49037
Phospholipase A2 (Bter07778; XM_003400908),PA2,MW3,MW,Medium,4,1,1,Control,121,30.30455208,29.8400631,29.93300629,30.02587382,29.8400631,30.30455208,0.464488983,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.48793,25.29041589,31.87633,26.22237753,1.993,1.29380822,1.901644416,0.680362852,-0.55562
Phospholipase A2 (Bter07778; XM_003400908),PA2,MW4,MW,Medium,4,1,1,Control,121,29.04710388,29.28736115,28.98045349,29.10497284,28.98045349,29.28736115,0.306907654,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,30.52219,24.63093495,31.87633,26.22237753,1.993,2.45486961,2.996737273,0.819180792,-0.28775
Phospholipase A2 (Bter07778; XM_003400908),PA2,MW5,MW,Medium,4,1,1,Control,121,29.53559685,29.64141846,29.74843597,29.64181709,29.53559685,29.74843597,0.212839127,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.08517,25.35247425,31.87633,26.22237753,1.993,1.68995212,1.821974874,0.927538652,-0.10852
Phospholipase A2 (Bter07778; XM_003400908),PA2,MW6,MW,Medium,4,2,1,Control,65,31.08344078,30.44115257,30.56969643,30.69809659,30.44115257,31.08344078,0.642288208,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,30.6981,24.95817119,31.87633,26.22237753,1.993,2.18454947,2.391327326,0.913530092,-0.13048
Phospholipase A2 (Bter07778; XM_003400908),PA2,MW7,MW,Medium,4,2,1,Control,65,NA,32.54095078,33.04457474,32.79276276,32.54095078,33.04457474,0.503623962,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,32.79276,25.70248086,31.87633,26.22237753,1.993,0.5445585,1.431241263,0.380479877,-1.39411
Phospholipase A2 (Bter07778; XM_003400908),PA2,MW8,MW,Medium,4,2,1,Control,65,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,30.75634,24.59640071,31.87633,26.22237753,1.993,2.10177975,3.068964925,0.684849714,-0.54614
Phospholipase A2 (Bter07778; XM_003400908),PA2,LQ3,LQ,Late,5,1,1,Control,121,30.82754707,30.06310081,30.29807854,30.39624214,30.06310081,30.82754707,0.764446259,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.87633,26.22237753,31.87633,26.22237753,1.993,1,1,1,0
Phospholipase A2 (Bter07778; XM_003400908),PA2,LQ4,LQ,Late,5,1,1,Control,121,30.15324974,29.9346447,29.83311653,29.97367032,29.83311653,30.15324974,0.320133209,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,31.43318,25.57874041,31.87633,26.22237753,1.993,1.34164616,1.558741203,0.860724129,-0.21638
Phospholipase A2 (Bter07778; XM_003400908),PA2,LQ5,LQ,Late,5,1,1,Control,121,31.18336105,31.12186241,31.3341217,31.21311506,31.12186241,31.3341217,0.212259293,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,32.73298,26.77102083,31.87633,26.22237753,1.993,0.56658255,0.684979114,0.827153034,-0.27377
Phospholipase A2 (Bter07778; XM_003400908),PA2,LQ6,LQ,Late,5,2,1,Control,65,32.8296814,NA,32.98488617,32.90728378,32.8296814,32.98488617,0.155204773,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,32.90728,26.04640784,31.87633,26.22237753,1.993,0.50473051,1.129026689,0.447049232,-1.16149
Phospholipase A2 (Bter07778; XM_003400908),PA2,LQ7,LQ,Late,5,2,1,Control,65,31.3485527,31.22663879,31.85759735,31.47759628,31.22663879,31.85759735,0.630958557,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,31.4776,24.94259361,31.87633,26.22237753,1.993,1.30270497,2.417155701,0.538941271,-0.8918
Phospholipase A2 (Bter07778; XM_003400908),PA2,LQ8,LQ,Late,5,2,1,Control,65,32.57296371,33.45378113,33.03846741,33.02173742,32.57296371,33.45378113,0.880817413,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,33.02174,25.43112605,31.87633,26.22237753,1.993,0.46783638,1.725780648,0.271086816,-1.88317
Phospholipase A2 (Bter07778; XM_003400908),PA2,LW3,LW,Late,6,1,1,Control,121,31.88539314,31.38862801,31.48800659,31.58734258,31.38862801,31.88539314,0.496765137,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,33.12543,26.28165682,31.87633,26.22237753,1.993,0.43674467,0.959942941,0.45496941,-1.13616
Phospholipase A2 (Bter07778; XM_003400908),PA2,LW4,LW,Late,6,1,1,Control,121,31.06599426,31.31542397,30.90855598,31.09665807,30.90855598,31.31542397,0.406867981,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,32.61085,25.96422745,31.87633,26.22237753,1.993,0.61438284,1.194862228,0.51418718,-0.95963
Phospholipase A2 (Bter07778; XM_003400908),PA2,LW5,LW,Late,6,1,1,Control,121,30.73991966,30.34951019,30.70002937,30.59648641,30.34951019,30.73991966,0.39040947,mw8,29.40106773,29.28207588,29.30160904,29.32825089,29.28207588,29.40106773,0.118991852,30.75633685,1.941,32.08633,25.98699537,31.87633,26.22237753,1.993,0.86999376,1.176247401,0.739635008,-0.43511
Phospholipase A2 (Bter07778; XM_003400908),PA2,LW6,LW,Late,6,2,1,Control,65,32.90286255,32.6778717,33.30805969,32.96293132,32.6778717,33.30805969,0.630187988,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,32.96293,25.6146532,31.87633,26.22237753,1.993,0.48644266,1.520610107,0.31989966,-1.64431
Phospholipase A2 (Bter07778; XM_003400908),PA2,LW7,LW,Late,6,2,1,Control,65,32.52821732,32.13199615,32.28403854,32.31475067,32.13199615,32.52821732,0.396221161,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,32.31475,25.49879946,31.87633,26.22237753,1.993,0.74769457,1.647088425,0.45394926,-1.1394
Phospholipase A2 (Bter07778; XM_003400908),PA2,LW8,LW,Late,6,2,1,Control,65,32.66399384,NA,32.80010986,32.73205185,32.66399384,32.80010986,0.136116028,mw8,30.46079254,30.95024109,30.85797691,30.75633685,30.46079254,30.95024109,0.489448547,30.75633685,1.941,32.73205,25.38152329,31.87633,26.22237753,1.993,0.56693182,1.78583768,0.317459882,-1.65535
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EQ3,EQ,Early,1,1,1,Control,116,22.10266113,22.08180237,22.16851044,22.11765798,22.08180237,22.16851044,0.086708069,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.16116,25.20347779,23.82684,26.22237753,1.993,1.59685021,2.019146866,0.790853919,-0.33852
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EQ4,EQ,Early,1,1,1,Control,116,22.3410244,22.32194901,22.3554039,22.3394591,22.32194901,22.3554039,0.033454895,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.39343,25.71323328,23.82684,26.22237753,1.993,1.35625562,1.420667419,0.954660891,-0.06694
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EQ5,EQ,Early,1,1,1,Control,116,22.15655136,22.16333389,22.09861755,22.13950094,22.09861755,22.16333389,0.064716339,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.18404,25.41256952,23.82684,26.22237753,1.993,1.57137459,1.748007975,0.898951613,-0.15368
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EQ6,EQ,Early,1,2,1,Control,71,22.27620125,22.27774429,22.23395157,22.26263237,22.23395157,22.27774429,0.043792725,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.26263,24.43307804,23.82684,26.22237753,1.993,3.00352685,3.434853817,0.874426398,-0.19359
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EQ7,EQ,Early,1,2,1,Control,71,22.64824867,22.64657974,22.49610329,22.59697723,22.49610329,22.64824867,0.152145386,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.59698,25.14049381,23.82684,26.22237753,1.993,2.37431785,2.108783874,1.125918059,0.1711
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EQ8,EQ,Early,1,2,1,Control,71,22.7389431,22.59047699,22.40195847,22.57712619,22.40195847,22.7389431,0.336984634,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.57713,24.92400707,23.82684,26.22237753,1.993,2.40768907,2.448338325,0.983397208,-0.02415
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EW3,EW,Early,2,1,1,Control,116,21.96050453,21.94771194,21.93632317,21.94817988,21.93632317,21.96050453,0.024181366,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,22.98369,24.92159537,23.82684,26.22237753,1.993,1.8090731,2.452413802,0.737670413,-0.43895
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EW4,EW,Early,2,1,1,Control,116,21.58567429,21.60601997,21.60273361,21.59814262,21.58567429,21.60601997,0.020345688,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,22.61714,24.71233996,23.82684,26.22237753,1.993,2.34090206,2.833135032,0.826258554,-0.27533
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EW5,EW,Early,2,1,1,Control,116,22.36401939,22.28722382,22.31834221,22.32319514,22.28722382,22.36401939,0.076795578,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.3764,25.76160569,23.82684,26.22237753,1.993,1.37259395,1.374056303,0.998935739,-0.00154
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EW6,EW,Early,2,2,1,Control,71,22.40208626,22.4005909,22.23493004,22.34586906,22.23493004,22.40208626,0.167156219,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.34587,24.47244454,23.82684,26.22237753,1.993,2.83279456,3.342856272,0.847417398,-0.23886
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EW7,EW,Early,2,2,1,Control,71,22.32283401,22.29423904,22.22128296,22.27945201,22.22128296,22.32283401,0.101551056,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.27945,24.71929254,23.82684,26.22237753,1.993,2.96821681,2.819583273,1.052714717,0.07411
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,EW8,EW,Early,2,2,1,Control,71,21.84645081,21.73234367,21.92106438,21.83328629,21.73234367,21.92106438,0.188720703,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,21.83329,24.24981237,23.82684,26.22237753,1.993,4.06193467,3.897603549,1.042162092,0.05958
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MQ3,MQ,Medium,3,1,1,Control,116,22.68690109,22.69065666,22.81877136,22.73210971,22.68690109,22.81877136,0.13187027,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.8046,25.89835631,23.82684,26.22237753,1.993,1.01575385,1.25039354,0.812347326,-0.29983
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MQ4,MQ,Medium,3,1,1,Control,116,22.38986397,22.42315292,22.52582932,22.44628207,22.38986397,22.52582932,0.135965347,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.50529,25.52299076,23.82684,26.22237753,1.993,1.25367244,1.619837474,0.773949525,-0.36969
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MQ5,MQ,Medium,3,1,1,Control,116,21.92627716,21.90995216,21.86778641,21.90133858,21.86778641,21.92627716,0.058490753,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,22.93464,25.12014958,23.82684,26.22237753,1.993,1.87255243,2.138579087,0.875605883,-0.19165
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MQ6,MQ,Medium,3,2,1,Control,71,23.07622528,23.11198807,23.17412567,23.12077967,23.07622528,23.17412567,0.097900391,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,23.12078,24.7511079,23.82684,26.22237753,1.993,1.64283923,2.75839205,0.595578583,-0.74764
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MQ7,MQ,Medium,3,2,1,Control,71,22.57452393,22.56281853,22.4812355,22.53952599,22.4812355,22.57452393,0.093288422,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.53953,24.76237619,23.82684,26.22237753,1.993,2.47218903,2.737039456,0.903234705,-0.14683
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MQ8,MQ,Medium,3,2,1,Control,71,22.28075409,22.15603256,22.18028641,22.20569102,22.15603256,22.28075409,0.124721527,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.20569,24.74641872,23.82684,26.22237753,1.993,3.12621353,2.767326721,1.129687183,0.17592
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MW3,MW,Medium,4,1,1,Control,116,21.99693298,21.94152451,21.99670601,21.97838783,21.94152451,21.99693298,0.055408478,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.01532,25.29041589,23.82684,26.22237753,1.993,1.76928131,1.901644416,0.930395448,-0.10408
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MW4,MW,Medium,4,1,1,Control,116,21.74370956,21.78490257,21.80696487,21.77852567,21.74370956,21.80696487,0.06325531,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,22.80603,24.63093495,23.82684,26.22237753,1.993,2.04976638,2.996737273,0.683999361,-0.54793
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MW5,MW,Medium,4,1,1,Control,116,22.06335258,22.03977585,21.99704933,22.03339259,21.99704933,22.06335258,0.066303253,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.07292,25.35247425,23.82684,26.22237753,1.993,1.69905991,1.821974874,0.932537508,-0.10077
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MW6,MW,Medium,4,2,1,Control,71,22.97146988,22.97136879,22.93859673,22.96047846,22.93859673,22.97146988,0.032873154,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.96048,24.95817119,23.82684,26.22237753,1.993,1.8388371,2.391327326,0.768960854,-0.37902
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MW7,MW,Medium,4,2,1,Control,71,23.84467888,23.73852348,23.72694588,23.77004941,23.72694588,23.84467888,0.117733002,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,23.77005,25.70248086,23.82684,26.22237753,1.993,1.0407338,1.431241263,0.727154694,-0.45967
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,MW8,MW,Medium,4,2,1,Control,71,22.58205986,22.38028526,22.45775795,22.47336769,22.38028526,22.58205986,0.201774597,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.47337,24.59640071,23.82684,26.22237753,1.993,2.58990122,3.068964925,0.84390056,-0.24486
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LQ3,LQ,Late,5,1,1,Control,116,22.80093765,22.68119621,22.77788544,22.75333977,22.68119621,22.80093765,0.11974144,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.82684,26.22237753,23.82684,26.22237753,1.993,1,1,1,0
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LQ4,LQ,Late,5,1,1,Control,116,22.21260643,21.99259949,22.07901764,22.09474119,21.99259949,22.21260643,0.220006943,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.13716,25.57874041,23.82684,26.22237753,1.993,1.62402235,1.558741203,1.041880686,0.05919
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LQ5,LQ,Late,5,1,1,Control,116,22.66474724,22.63156891,22.70362282,22.66664632,22.63156891,22.70362282,0.072053909,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.73605,26.77102083,23.82684,26.22237753,1.993,1.0659109,0.684979114,1.556121748,0.63795
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LQ6,LQ,Late,5,2,1,Control,71,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,24.02754,26.04640784,23.82684,26.22237753,1.993,0.86839169,1.129026689,0.76915072,-0.37866
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LQ7,LQ,Late,5,2,1,Control,71,22.43199921,22.33936119,22.3145752,22.36197853,22.3145752,22.43199921,0.117424011,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.36198,24.94259361,23.82684,26.22237753,1.993,2.80088987,2.417155701,1.158754427,0.21257
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LQ8,LQ,Late,5,2,1,Control,71,23.07364845,23.02852631,23.22756386,23.10991287,23.02852631,23.22756386,0.199037552,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,23.10991,25.43112605,23.82684,26.22237753,1.993,1.65543929,1.725780648,0.959240844,-0.06004
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LW3,LW,Late,6,1,1,Control,116,22.93077278,22.94536209,22.96056366,22.94556618,22.93077278,22.96056366,0.029790878,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,24.02813,26.28165682,23.82684,26.22237753,1.993,0.86802876,0.959942941,0.904250371,-0.14521
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LW4,LW,Late,6,1,1,Control,116,22.42598152,22.40384865,22.43489647,22.42157555,22.40384865,22.43489647,0.031047821,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.47942,25.96422745,23.82684,26.22237753,1.993,1.27668625,1.194862228,1.068479882,0.09556
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LW5,LW,Late,6,1,1,Control,116,22.83018494,22.76849747,22.7936554,22.79744593,22.76849747,22.83018494,0.061687469,lq6,22.96159554,22.93584251,22.93755722,22.94499842,22.93584251,22.96159554,0.025753021,24.02753639,2.02,23.87302,25.98699537,23.82684,26.22237753,1.993,0.9680476,1.176247401,0.822996588,-0.28104
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LW6,LW,Late,6,2,1,Control,71,23.27482986,23.33304214,23.62517929,23.4110171,23.27482986,23.62517929,0.350349426,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,23.41102,25.6146532,23.82684,26.22237753,1.993,1.33958621,1.520610107,0.880953116,-0.18286
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LW7,LW,Late,6,2,1,Control,71,23.22806549,23.18474197,23.16762733,23.19347827,23.16762733,23.22806549,0.060438156,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,23.19348,25.49879946,23.82684,26.22237753,1.993,1.5609768,1.647088425,0.947718884,-0.07747
Ribosomal protein 18S (Bter07703; XM_003400778),RP18S,LW8,LW,Late,6,2,1,Control,71,22.61252403,22.77902794,22.53116798,22.64090665,22.53116798,22.77902794,0.247859955,lq6,24.06656647,23.9568615,24.05918121,24.02753639,23.9568615,24.06656647,0.109704971,24.02753639,2.02,22.64091,25.38152329,23.82684,26.22237753,1.993,2.30210405,1.78583768,1.289089192,0.36635
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EQ3,EQ,Early,1,1,1,Control,131,23.00222015,22.95056725,22.98786926,22.98021889,22.95056725,23.00222015,0.051652908,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,22.98022,25.20347779,24.16923,26.22237753,1.993,2.2298931,2.019146866,1.1043739,0.14323
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EQ4,EQ,Early,1,1,1,Control,131,23.94218826,23.96899414,23.91613579,23.9424394,23.91613579,23.96899414,0.052858353,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.94244,25.71323328,24.16923,26.22237753,1.993,1.16527955,1.420667419,0.820233878,-0.28589
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EQ5,EQ,Early,1,1,1,Control,131,23.15739822,23.11906433,23.54394913,23.27347056,23.11906433,23.54394913,0.424884796,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.27347,25.41256952,24.16923,26.22237753,1.993,1.82972017,1.748007975,1.046745893,0.06591
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EQ6,EQ,Early,1,2,1,Control,23,22.11611176,22.1179657,22.04029465,22.09145737,22.04029465,22.1179657,0.077671051,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.45665,24.43307804,24.16923,26.22237753,1.993,3.1743026,3.434853817,0.924144888,-0.11381
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EQ7,EQ,Early,1,2,1,Control,23,22.89861298,22.98456383,22.97130394,22.95149358,22.89861298,22.98456383,0.085950851,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.3309,25.14049381,24.16923,26.22237753,1.993,1.76019897,2.108783874,0.83469861,-0.26067
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EQ8,EQ,Early,1,2,1,Control,23,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.99048,24.92400707,24.16923,26.22237753,1.993,2.21451014,2.448338325,0.904495149,-0.14482
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EW3,EW,Early,2,1,1,Control,131,22.70775986,22.82052231,22.81703186,22.78177134,22.70775986,22.82052231,0.112762451,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,22.78177,24.92159537,24.16923,26.22237753,1.993,2.54925546,2.452413802,1.039488302,0.05587
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EW4,EW,Early,2,1,1,Control,131,22.4332684,22.32772097,22.25867125,22.33988687,22.25867125,22.4332684,0.174597155,EQ8,23.00677872,23.10577202,23.17715645,23.09656906,23.00677872,23.17715645,0.170377731,22.99048233,1.963,22.23728,24.71233996,24.16923,26.22237753,1.993,3.68050302,2.833135032,1.299091989,0.3775
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EW5,EW,Early,2,1,1,Control,131,24.16993904,24.21772003,24.1051178,24.16425896,24.1051178,24.21772003,0.112602234,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,24.16426,25.76160569,24.16923,26.22237753,1.993,1.00335512,1.374056303,0.730213979,-0.45361
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EW6,EW,Early,2,2,1,Control,23,22.03810692,22.01754951,22.06185341,22.03916995,22.01754951,22.06185341,0.044303894,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.4035,24.47244454,24.16923,26.22237753,1.993,3.2901641,3.342856272,0.984237379,-0.02292
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EW7,EW,Early,2,2,1,Control,23,22.51380348,22.56316566,22.54711342,22.54136086,22.51380348,22.56316566,0.049362183,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.91399,24.71929254,24.16923,26.22237753,1.993,2.3317604,2.819583273,0.826987597,-0.27406
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,EW8,EW,Early,2,2,1,Control,23,22.07732391,22.09625816,22.08643532,22.08667247,22.07732391,22.09625816,0.01893425,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.45178,24.24981237,24.16923,26.22237753,1.993,3.18473345,3.897603549,0.817100407,-0.29141
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MQ3,MQ,Medium,3,1,1,Control,131,23.88947105,23.98023796,23.85728264,23.90899722,23.85728264,23.98023796,0.122955322,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.909,25.89835631,24.16923,26.22237753,1.993,1.19186212,1.25039354,0.9531896,-0.06916
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MQ4,MQ,Medium,3,1,1,Control,131,23.06689644,23.13906288,23.42529869,23.21041934,23.06689644,23.42529869,0.358402252,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.21042,25.52299076,24.16923,26.22237753,1.993,1.90920982,1.619837474,1.178642825,0.23713
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MQ5,MQ,Medium,3,1,1,Control,131,23.04243469,22.954216,22.91259575,22.96974881,22.91259575,23.04243469,0.129838943,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,22.96975,25.12014958,24.16923,26.22237753,1.993,2.24569587,2.138579087,1.050087829,0.07051
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MQ6,MQ,Medium,3,2,1,Control,23,22.13017654,22.14640427,22.09604836,22.12420972,22.09604836,22.14640427,0.050355911,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.48994,24.7511079,24.16923,26.22237753,1.993,3.10381552,2.75839205,1.125226386,0.17022
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MQ7,MQ,Medium,3,2,1,Control,23,22.40862656,22.40924072,22.45717239,22.42501322,22.40862656,22.45717239,0.048545837,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.79572,24.76237619,24.16923,26.22237753,1.993,2.52538669,2.737039456,0.922670911,-0.11611
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MQ8,MQ,Medium,3,2,1,Control,23,22.71056557,22.49227905,22.55447197,22.5857722,22.49227905,22.71056557,0.218286514,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.95914,24.74641872,24.16923,26.22237753,1.993,2.26182963,2.767326721,0.817333786,-0.291
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MW3,MW,Medium,4,1,1,Control,131,23.40284729,23.38587761,23.55114937,23.44662476,23.38587761,23.55114937,0.165271759,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.44662,25.29041589,24.16923,26.22237753,1.993,1.6280371,1.901644416,0.856120673,-0.22411
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MW4,MW,Medium,4,1,1,Control,131,22.34172249,22.35059357,22.30164909,22.33132172,22.30164909,22.35059357,0.048944473,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,22.33132,24.63093495,24.16923,26.22237753,1.993,3.4542932,2.996737273,1.152684697,0.205
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MW5,MW,Medium,4,1,1,Control,131,23.36494064,23.40412712,23.26970291,23.34625689,23.26970291,23.40412712,0.13442421,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.34626,25.35247425,24.16923,26.22237753,1.993,1.7420639,1.821974874,0.956140464,-0.06471
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MW6,MW,Medium,4,2,1,Control,23,22.34125328,22.32448006,22.43471718,22.36681684,22.32448006,22.43471718,0.110237122,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.73656,24.95817119,24.16923,26.22237753,1.993,2.62818898,2.391327326,1.099050286,0.13626
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MW7,MW,Medium,4,2,1,Control,23,22.75924301,22.94736481,22.95880318,22.88847033,22.75924301,22.95880318,0.199560165,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.26684,25.70248086,24.16923,26.22237753,1.993,1.83792474,1.431241263,1.284147397,0.36081
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,MW8,MW,Medium,4,2,1,Control,23,22.23554802,22.23601723,22.21175194,22.22777239,22.21175194,22.23601723,0.024265289,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.59522,24.59640071,24.16923,26.22237753,1.993,2.89107139,3.068964925,0.94203468,-0.08615
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LQ3,LQ,Late,5,1,1,Control,131,24.17894363,24.13760567,24.19112587,24.16922506,24.13760567,24.19112587,0.053520203,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,24.16923,26.22237753,24.16923,26.22237753,1.993,1,1,1,0
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LQ4,LQ,Late,5,1,1,Control,131,23.58798981,23.42746162,23.53443527,23.5166289,23.42746162,23.58798981,0.160528183,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.51663,25.57874041,24.16923,26.22237753,1.993,1.55295425,1.558741203,0.996287419,-0.00537
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LQ5,LQ,Late,5,1,1,Control,131,25.44577599,25.3928318,24.9912796,25.27662913,24.9912796,25.44577599,0.454496384,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,25.27663,26.77102083,24.16923,26.22237753,1.993,0.47382598,0.684979114,0.691737847,-0.5317
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LQ6,LQ,Late,5,2,1,Control,23,23.33759689,23.26947784,23.30964661,23.30557378,23.26947784,23.33759689,0.068119049,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.69084,26.04640784,24.16923,26.22237753,1.993,1.38079735,1.129026689,1.22299797,0.29042
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LQ7,LQ,Late,5,2,1,Control,23,22.5136776,22.54583168,22.60235214,22.55395381,22.5136776,22.60235214,0.088674545,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,22.92679,24.94259361,24.16923,26.22237753,1.993,2.31171458,2.417155701,0.95637802,-0.06435
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LQ8,LQ,Late,5,2,1,Control,23,23.17152214,23.17590523,23.16789436,23.17177391,23.16789436,23.17590523,0.008010864,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.55482,25.43112605,24.16923,26.22237753,1.993,1.51345856,1.725780648,0.876970406,-0.1894
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LW3,LW,Late,6,1,1,Control,131,24.44840813,24.44467926,24.49205017,24.46171252,24.44467926,24.49205017,0.047370911,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,24.46171,26.28165682,24.16923,26.22237753,1.993,0.8209647,0.959942941,0.855222396,-0.22563
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LW4,LW,Late,6,1,1,Control,131,23.92955399,23.93728256,24.00510216,23.9573129,23.92955399,24.00510216,0.075548172,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.95731,25.96422745,24.16923,26.22237753,1.993,1.15364814,1.194862228,0.965507251,-0.05064
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LW5,LW,Late,6,1,1,Control,131,23.97618675,23.9911232,23.78606224,23.91779073,23.78606224,23.9911232,0.205060959,EQ8,23.01318359,22.9841423,22.97412109,22.99048233,22.97412109,23.01318359,0.0390625,22.99048233,1.963,23.91779,25.98699537,24.16923,26.22237753,1.993,1.18481411,1.176247401,1.007283085,0.01047
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LW6,LW,Late,6,2,1,Control,23,22.95176125,22.93391228,22.98019791,22.95529048,22.93391228,22.98019791,0.046285629,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.33476,25.6146532,24.16923,26.22237753,1.993,1.7556227,1.520610107,1.154551516,0.20733
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LW7,LW,Late,6,2,1,Control,23,23.22722435,23.23097038,23.21470642,23.22430038,23.21470642,23.23097038,0.016263962,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.60822,25.49879946,24.16923,26.22237753,1.993,1.45992356,1.647088425,0.886366231,-0.17403
Ribosomal protein S5a (Bter10521; XM_012308702; XM_012308701),RPS5a,LW8,LW,Late,6,2,1,Control,23,23.0330658,22.99742126,22.98984909,23.00677872,22.98984909,23.0330658,0.043216705,EQ8,22.60146904,22.65632629,22.59203339,22.61660957,22.59203339,22.65632629,0.064292908,22.99048233,1.963,23.3871,25.38152329,24.16923,26.22237753,1.993,1.69472767,1.78583768,0.948981919,-0.07555
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EQ3,EQ,Early,1,1,1,Control,122,29.93016624,30.00736237,30.30020142,30.07924334,29.93016624,30.30020142,0.370035172,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.07924,25.20347779,31.31029,26.22237753,1.993,2.34304053,2.019146866,1.160411145,0.21464
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EQ4,EQ,Early,1,1,1,Control,122,30.66065598,30.09192467,30.30771637,30.35343234,30.09192467,30.66065598,0.568731308,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.35343,25.71323328,31.31029,26.22237753,1.993,1.93829306,1.420667419,1.364353848,0.44822
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EQ5,EQ,Early,1,1,1,Control,122,30.32364464,30.91099167,30.01198769,30.41554133,30.01198769,30.91099167,0.899003983,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.41554,25.41256952,31.31029,26.22237753,1.993,1.85679202,1.748007975,1.062233152,0.0871
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EQ6,EQ,Early,1,2,1,Control,25,28.89063072,28.65345573,28.78640556,28.77683067,28.65345573,28.89063072,0.237174988,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.17515,24.43307804,31.31029,26.22237753,1.993,4.37874758,3.434853817,1.274798817,0.35027
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EQ7,EQ,Early,1,2,1,Control,25,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.13974,25.14049381,31.31029,26.22237753,1.993,2.24702219,2.108783874,1.065553573,0.0916
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EQ8,EQ,Early,1,2,1,Control,25,29.60889816,29.40029144,29.25570107,29.42163022,29.25570107,29.60889816,0.353197098,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.82887,24.92400707,31.31029,26.22237753,1.993,2.78602733,2.448338325,1.137925791,0.18641
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EW3,EW,Early,2,1,1,Control,122,29.32170296,29.66352654,29.69801331,29.56108093,29.32170296,29.69801331,0.376310349,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,29.56108,24.92159537,31.31029,26.22237753,1.993,3.3529295,2.452413802,1.367195659,0.45122
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EW4,EW,Early,2,1,1,Control,122,29.76888657,30.10099983,30.15107918,30.00698853,29.76888657,30.15107918,0.382192612,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.00699,24.71233996,31.31029,26.22237753,1.993,2.46310861,2.833135032,0.869393299,-0.20192
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EW5,EW,Early,2,1,1,Control,122,30.48864746,30.19138145,30.12062836,30.26688576,30.12062836,30.48864746,0.368019104,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.26689,25.76160569,31.31029,26.22237753,1.993,2.05786147,1.374056303,1.497654398,0.5827
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EW6,EW,Early,2,2,1,Control,25,28.91651154,28.78268242,28.9311676,28.87678719,28.78268242,28.9311676,0.148485184,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.27649,24.47244454,31.31029,26.22237753,1.993,4.08234379,3.342856272,1.221214272,0.28832
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EW7,EW,Early,2,2,1,Control,25,29.24736023,29.21671677,29.085495,29.18319066,29.085495,29.24736023,0.161865234,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.58713,24.71929254,31.31029,26.22237753,1.993,3.29305323,2.819583273,1.167921962,0.22394
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,EW8,EW,Early,2,2,1,Control,25,28.68439293,28.75141907,28.64498901,28.69360034,28.64498901,28.75141907,0.106430054,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.09077,24.24981237,31.31029,26.22237753,1.993,4.64190798,3.897603549,1.190964633,0.25213
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MQ3,MQ,Medium,3,1,1,Control,122,30.41025543,30.30241585,30.84938622,30.52068583,30.30241585,30.84938622,0.546970367,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.52069,25.89835631,31.31029,26.22237753,1.993,1.72655397,1.25039354,1.380808457,0.46551
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MQ4,MQ,Medium,3,1,1,Control,122,30.62309265,30.38425636,30.41814804,30.47516569,30.38425636,30.62309265,0.238836288,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.47517,25.52299076,31.31029,26.22237753,1.993,1.78177727,1.619837474,1.099972864,0.13747
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MQ5,MQ,Medium,3,1,1,Control,122,30.18767738,30.04660606,30.03495979,30.08974775,30.03495979,30.18767738,0.15271759,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.08975,25.12014958,31.31029,26.22237753,1.993,2.32607925,2.138579087,1.087675115,0.12125
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MQ6,MQ,Medium,3,2,1,Control,25,28.78410721,28.74051666,28.76215172,28.76225853,28.74051666,28.78410721,0.043590546,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.16038,24.7511079,31.31029,26.22237753,1.993,4.4237202,2.75839205,1.603731492,0.68143
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MQ7,MQ,Medium,3,2,1,Control,25,29.11377144,29.14326286,29.18696594,29.14800008,29.11377144,29.18696594,0.073194504,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.55146,24.76237619,31.31029,26.22237753,1.993,3.37532458,2.737039456,1.23320275,0.30241
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MQ8,MQ,Medium,3,2,1,Control,25,29.26716042,29.17239571,29.51703453,29.31886355,29.17239571,29.51703453,0.344638824,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.72468,24.74641872,31.31029,26.22237753,1.993,2.99420479,2.767326721,1.081984563,0.11368
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MW3,MW,Medium,4,1,1,Control,122,29.88480377,29.94456863,30.09797668,29.97578303,29.88480377,30.09797668,0.213172913,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,29.97578,25.29041589,31.31029,26.22237753,1.993,2.51684812,1.901644416,1.32351143,0.40437
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MW4,MW,Medium,4,1,1,Control,122,29.21837807,29.69626045,29.10910416,29.34124756,29.10910416,29.69626045,0.587156296,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,29.34125,24.63093495,31.31029,26.22237753,1.993,3.90353022,2.996737273,1.302593408,0.38139
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MW5,MW,Medium,4,1,1,Control,122,30.5180645,30.35146141,29.88378143,30.25110245,29.88378143,30.5180645,0.634283066,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.2511,25.35247425,31.31029,26.22237753,1.993,2.0804491,1.821974874,1.141864869,0.19139
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MW6,MW,Medium,4,2,1,Control,25,29.17254639,29.92666817,29.02265739,29.37395732,29.02265739,29.92666817,0.904010773,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.78054,24.95817119,31.31029,26.22237753,1.993,2.88073612,2.391327326,1.204659893,0.26863
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MW7,MW,Medium,4,2,1,Control,25,30.64100075,30.20707703,29.99855423,30.28221067,29.99855423,30.64100075,0.642446518,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.70137,25.70248086,31.31029,26.22237753,1.993,1.52372933,1.431241263,1.06462088,0.09034
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,MW8,MW,Medium,4,2,1,Control,25,29.0267868,28.94169235,28.74390602,28.90412839,28.74390602,29.0267868,0.282880783,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,29.30421,24.59640071,31.31029,26.22237753,1.993,4.00482183,3.068964925,1.304942196,0.38399
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LQ3,LQ,Late,5,1,1,Control,122,31.41435814,31.26503563,31.2514801,31.31029129,31.2514801,31.41435814,0.162878036,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,31.31029,26.22237753,31.31029,26.22237753,1.993,1,1,1,0
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LQ4,LQ,Late,5,1,1,Control,122,30.58612251,30.92669106,30.75994873,30.75758743,30.58612251,30.92669106,0.340568542,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.75759,25.57874041,31.31029,26.22237753,1.993,1.46561572,1.558741203,0.94025597,-0.08887
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LQ5,LQ,Late,5,1,1,Control,122,31.86176682,32.19898605,31.87675095,31.97916794,31.86176682,32.19898605,0.337219238,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,31.97917,26.77102083,31.31029,26.22237753,1.993,0.62962814,0.684979114,0.919193189,-0.12156
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LQ6,LQ,Late,5,2,1,Control,25,30.25445557,30.66176033,30.93925285,30.61848958,30.25445557,30.93925285,0.684797287,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,31.0423,26.04640784,31.31029,26.22237753,1.993,1.20364624,1.129026689,1.066091929,0.09233
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LQ7,LQ,Late,5,2,1,Control,25,29.83623886,29.84575844,29.87973213,29.85390981,29.83623886,29.87973213,0.043493271,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.26714,24.94259361,31.31029,26.22237753,1.993,2.05750437,2.417155701,0.851208871,-0.23241
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LQ8,LQ,Late,5,2,1,Control,25,29.76493645,29.87771225,29.76398277,29.80221049,29.76398277,29.87771225,0.113729477,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.21472,25.43112605,31.31029,26.22237753,1.993,2.13346273,1.725780648,1.236230533,0.30595
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LW3,LW,Late,6,1,1,Control,122,30.69132042,31.00723648,30.95724678,30.88526789,30.69132042,31.00723648,0.315916061,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.88527,26.28165682,31.31029,26.22237753,1.993,1.34173796,0.959942941,1.397726779,0.48308
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LW4,LW,Late,6,1,1,Control,122,30.99901581,30.97287178,31.37976837,31.11721865,30.97287178,31.37976837,0.406896591,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,31.11722,25.96422745,31.31029,26.22237753,1.993,1.14286461,1.194862228,0.956482333,-0.06419
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LW5,LW,Late,6,1,1,Control,122,30.9402504,30.55802155,30.70811272,30.73546155,30.55802155,30.9402504,0.382228851,eq7,30.28770065,30.13851547,29.99300957,30.1397419,29.99300957,30.28770065,0.294691086,30.1397419,1.997,30.73546,25.98699537,31.31029,26.22237753,1.993,1.48821694,1.176247401,1.265224425,0.33939
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LW6,LW,Late,6,2,1,Control,25,30.41369057,30.54808426,30.06998062,30.34391848,30.06998062,30.54808426,0.478103638,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.76393,25.6146532,31.31029,26.22237753,1.993,1.45920258,1.520610107,0.959616523,-0.05947
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LW7,LW,Late,6,2,1,Control,25,29.89184189,30.00531006,29.69732285,29.86482493,29.69732285,30.00531006,0.307987213,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.2782,25.49879946,31.31029,26.22237753,1.993,2.04181656,1.647088425,1.239652062,0.30994
TATA-binding protein (Bter08799; XM_003402044; XM_012317984),TBP,LW8,LW,Late,6,2,1,Control,25,30.66467285,30.43214607,30.27956581,30.45879491,30.27956581,30.66467285,0.38510704,eq7,29.85289192,29.79227257,29.53959846,29.72825432,29.53959846,29.85289192,0.313293457,30.1397419,1.997,30.88039,25.38152329,31.31029,26.22237753,1.993,1.34626826,1.78583768,0.753858133,-0.40764
chymotrypsin 2,Chymotry2,EQ3,EQ,Early,1,1,1,Target,114,29.57045746,29.41567039,29.35970688,29.44861158,29.35970688,29.57045746,0.21075058,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,32.09033,25.20347779,28.10006,26.22237753,1.993,0.0845995,2.019146866,0.041898637,-4.57695
chymotrypsin 2,Chymotry2,EQ4,EQ,Early,1,1,1,Target,114,29.32081223,28.90813255,29.24750328,29.15881602,28.90813255,29.32081223,0.412679672,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,31.77454,25.71323328,28.10006,26.22237753,1.993,0.10286235,1.420667419,0.072404243,-3.78778
chymotrypsin 2,Chymotry2,EQ5,EQ,Early,1,1,1,Target,114,28.45073891,28.39951897,27.998703,28.28298696,27.998703,28.45073891,0.452035904,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.82014,25.41256952,28.10006,26.22237753,1.993,0.18569896,1.748007975,0.106234621,-3.23467
chymotrypsin 2,Chymotry2,EQ6,EQ,Early,1,2,1,Target,58,29.52674675,29.82597542,29.90149307,29.75140508,29.52674675,29.90149307,0.374746323,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,29.75141,24.43307804,28.10006,26.22237753,1.993,0.35983052,3.434853817,0.104758611,-3.25486
chymotrypsin 2,Chymotry2,EQ7,EQ,Early,1,2,1,Target,58,31.3710556,31.09569168,30.85744286,31.10806338,30.85744286,31.3710556,0.513612747,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,31.10806,25.14049381,28.10006,26.22237753,1.993,0.15538607,2.108783874,0.073685158,-3.76248
chymotrypsin 2,Chymotry2,EQ8,EQ,Early,1,2,1,Target,58,30.81302071,31.0744648,30.80662346,30.89803632,30.80662346,31.0744648,0.267841339,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,30.89804,24.92400707,28.10006,26.22237753,1.993,0.17695786,2.448338325,0.07227672,-3.79033
chymotrypsin 2,Chymotry2,EW3,EW,Early,2,1,1,Target,114,28.07791901,28.08514977,28.36127281,28.17478053,28.07791901,28.36127281,0.283353806,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.70223,24.92159537,28.10006,26.22237753,1.993,0.19975881,2.452413802,0.081453959,-3.61787
chymotrypsin 2,Chymotry2,EW4,EW,Early,2,1,1,Target,114,26.96387863,26.93682289,26.89285851,26.93118668,26.89285851,26.96387863,0.071020126,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,29.34707,24.71233996,28.10006,26.22237753,1.993,0.46215413,2.833135032,0.16312464,-2.61595
chymotrypsin 2,Chymotry2,EW5,EW,Early,2,1,1,Target,114,27.47232056,27.53792763,27.70078087,27.57034302,27.47232056,27.70078087,0.228460312,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.04357,25.76160569,28.10006,26.22237753,1.993,0.30030418,1.374056303,0.218553036,-2.19394
chymotrypsin 2,Chymotry2,EW6,EW,Early,2,2,1,Target,58,30.38126183,30.11675835,30.06489944,30.18763987,30.06489944,30.38126183,0.316362381,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,30.18764,24.47244454,28.10006,26.22237753,1.993,0.2746837,3.342856272,0.082170359,-3.60524
chymotrypsin 2,Chymotry2,EW7,EW,Early,2,2,1,Target,58,30.68241882,30.66372299,30.45872498,30.60162226,30.45872498,30.68241882,0.223693848,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,30.60162,24.71929254,28.10006,26.22237753,1.993,0.21259326,2.819583273,0.075398823,-3.72931
chymotrypsin 2,Chymotry2,EW8,EW,Early,2,2,1,Target,58,28.89706612,28.88331604,28.82901573,28.8697993,28.82901573,28.89706612,0.068050385,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,28.8698,24.24981237,28.10006,26.22237753,1.993,0.6209893,3.897603549,0.159325929,-2.64995
chymotrypsin 2,Chymotry2,MQ3,MQ,Medium,3,1,1,Target,114,29.53447533,29.78111649,29.76631927,29.69397036,29.53447533,29.78111649,0.246641159,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,32.3577,25.89835631,28.10006,26.22237753,1.993,0.07169615,1.25039354,0.05733887,-4.12434
chymotrypsin 2,Chymotry2,MQ4,MQ,Medium,3,1,1,Target,114,29.54281998,29.54288864,29.52297401,29.53622754,29.52297401,29.54288864,0.019914627,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,32.1858,25.52299076,28.10006,26.22237753,1.993,0.07974488,1.619837474,0.049230174,-4.34431
chymotrypsin 2,Chymotry2,MQ5,MQ,Medium,3,1,1,Target,114,28.25574303,28.17436218,27.93128586,28.12046369,27.93128586,28.25574303,0.324457169,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.64304,25.12014958,28.10006,26.22237753,1.993,0.20721289,2.138579087,0.096892787,-3.36747
chymotrypsin 2,Chymotry2,MQ6,MQ,Medium,3,2,1,Target,58,30.17556381,30.21697617,30.24695969,30.21316655,30.17556381,30.24695969,0.071395874,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,30.21317,24.7511079,28.10006,26.22237753,1.993,0.27037779,2.75839205,0.098020072,-3.35078
chymotrypsin 2,Chymotry2,MQ7,MQ,Medium,3,2,1,Target,58,28.53684807,28.47146034,28.61120033,28.53983625,28.47146034,28.61120033,0.13973999,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,28.53984,24.76237619,28.10006,26.22237753,1.993,0.76169677,2.737039456,0.278292215,-1.84533
chymotrypsin 2,Chymotry2,MQ8,MQ,Medium,3,2,1,Target,58,29.29718018,29.0976944,29.25282288,29.21589915,29.0976944,29.29718018,0.199485779,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,29.2159,24.74641872,28.10006,26.22237753,1.993,0.50124305,2.767326721,0.181128972,-2.46491
chymotrypsin 2,Chymotry2,MW3,MW,Medium,4,1,1,Target,114,27.32230186,27.19105339,27.26947021,27.26094182,27.19105339,27.32230186,0.131248474,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,29.70641,25.29041589,28.10006,26.22237753,1.993,0.36999257,1.901644416,0.194564541,-2.36168
chymotrypsin 2,Chymotry2,MW4,MW,Medium,4,1,1,Target,114,27.82502365,28.20936966,27.75691986,27.93043772,27.75691986,28.20936966,0.452449799,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.43596,24.63093495,28.10006,26.22237753,1.993,0.23554847,2.996737273,0.078601643,-3.6693
chymotrypsin 2,Chymotry2,MW5,MW,Medium,4,1,1,Target,114,28.29733276,28.46711922,28.14978981,28.30474726,28.14978981,28.46711922,0.317329407,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.84385,25.35247425,28.10006,26.22237753,1.993,0.18299336,1.821974874,0.100436817,-3.31564
chymotrypsin 2,Chymotry2,MW6,MW,Medium,4,2,1,Target,58,30.45195198,30.62497902,30.38091087,30.48594729,30.38091087,30.62497902,0.244068146,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,30.48595,24.95817119,28.10006,26.22237753,1.993,0.22837276,2.391327326,0.095500418,-3.38835
chymotrypsin 2,Chymotry2,MW7,MW,Medium,4,2,1,Target,58,31.80997276,31.67367744,31.69147491,31.72504171,31.67367744,31.80997276,0.136295319,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,31.72504,25.70248086,28.10006,26.22237753,1.993,0.10606228,1.431241263,0.074105104,-3.75428
chymotrypsin 2,Chymotry2,MW8,MW,Medium,4,2,1,Target,58,31.75255203,32.11532593,31.54037666,31.80275154,31.54037666,32.11532593,0.574949265,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,31.80275,24.59640071,28.10006,26.22237753,1.993,0.10108149,3.068964925,0.032936672,-4.92416
chymotrypsin 2,Chymotry2,LQ3,LQ,Late,5,1,1,Target,114,25.73594475,25.79620552,25.82832527,25.78682518,25.73594475,25.82832527,0.092380524,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,28.10006,26.22237753,28.10006,26.22237753,1.993,1,1,1,0
chymotrypsin 2,Chymotry2,LQ4,LQ,Late,5,1,1,Target,114,25.44696999,25.31250191,25.51926613,25.42624601,25.31250191,25.51926613,0.206764221,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,27.70713,25.57874041,28.10006,26.22237753,1.993,1.27533122,1.558741203,0.818180219,-0.28951
chymotrypsin 2,Chymotry2,LQ5,LQ,Late,5,1,1,Target,114,24.79040527,24.75824738,24.74090004,24.76318423,24.74090004,24.79040527,0.049505234,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,26.98459,26.77102083,28.10006,26.22237753,1.993,1.9945778,0.684979114,2.911881193,1.54195
chymotrypsin 2,Chymotry2,LQ6,LQ,Late,5,2,1,Target,58,28.20557022,29.06518364,29.14469337,28.80514908,28.20557022,29.14469337,0.939123154,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,28.80515,26.04640784,28.10006,26.22237753,1.993,0.64634273,1.129026689,0.572477814,-0.80471
chymotrypsin 2,Chymotry2,LQ7,LQ,Late,5,2,1,Target,58,27.37611198,26.92327309,26.91369057,27.07102521,26.91369057,27.37611198,0.462421417,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,27.07103,24.94259361,28.10006,26.22237753,1.993,1.89067119,2.417155701,0.782188416,-0.35441
chymotrypsin 2,Chymotry2,LQ8,LQ,Late,5,2,1,Target,58,28.54307175,28.28681183,28.4086647,28.41284943,28.28681183,28.54307175,0.256259918,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,28.41285,25.43112605,28.10006,26.22237753,1.993,0.82398192,1.725780648,0.477454609,-1.06656
chymotrypsin 2,Chymotry2,LW3,LW,Late,6,1,1,Target,114,30.42106438,30.37414742,30.78105354,30.52542178,30.37414742,30.78105354,0.406906128,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,33.26373,26.28165682,28.10006,26.22237753,1.993,0.04092062,0.959942941,0.042628176,-4.55205
chymotrypsin 2,Chymotry2,LW4,LW,Late,6,1,1,Target,114,26.86265945,27.1104641,26.88041496,26.9511795,26.86265945,27.1104641,0.247804642,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,29.36886,25.96422745,28.10006,26.22237753,1.993,0.45596386,1.194862228,0.381603712,-1.38985
chymotrypsin 2,Chymotry2,LW5,LW,Late,6,1,1,Target,114,28.36895752,28.30659294,28.26193619,28.31249555,28.26193619,28.36895752,0.107021332,lw8,26.20892525,26.16324806,26.32615089,26.23277473,26.16324806,26.32615089,0.162902832,28.58601093,1.857,30.8523,25.98699537,28.10006,26.22237753,1.993,0.18203951,1.176247401,0.154762941,-2.69187
chymotrypsin 2,Chymotry2,LW6,LW,Late,6,2,1,Target,58,27.49920845,27.21578979,27.35351181,27.35617002,27.21578979,27.49920845,0.283418655,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,27.35617,25.6146532,28.10006,26.22237753,1.993,1.58476798,1.520610107,1.042192192,0.05962
chymotrypsin 2,Chymotry2,LW7,LW,Late,6,2,1,Target,58,28.66552353,28.5333004,28.70299339,28.63393911,28.5333004,28.70299339,0.169692993,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,28.63394,25.49879946,28.10006,26.22237753,1.993,0.71859824,1.647088425,0.436283945,-1.19666
chymotrypsin 2,Chymotry2,LW8,LW,Late,6,2,1,Target,58,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,lw8,28.70400047,NA,28.46802139,28.58601093,28.46802139,28.70400047,0.23597908,28.58601093,1.857,28.58601,25.38152329,28.10006,26.22237753,1.993,0.74023533,1.78583768,0.414503144,-1.27055
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EQ3,EQ,Early,1,1,1,Target,117,32.5723114,32.92770004,NA,32.75000572,32.5723114,32.92770004,0.355388641,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.75001,25.20347779,31.99603,26.22237753,1.993,0.60509592,2.019146866,0.299679005,-1.73851
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EQ4,EQ,Early,1,1,1,Target,117,33.65901947,33.22420883,33.43947983,33.44090271,33.22420883,33.65901947,0.434810638,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,33.4409,25.71323328,31.99603,26.22237753,1.993,0.38185824,1.420667419,0.268787918,-1.89546
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EQ5,EQ,Early,1,1,1,Target,117,33.17416,33.4756813,33.39875412,33.34953181,33.17416,33.4756813,0.301521301,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,33.34953,25.41256952,31.99603,26.22237753,1.993,0.4058278,1.748007975,0.232165873,-2.10677
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EQ6,EQ,Early,1,2,1,Target,41,32.17938995,32.3589859,31.93613243,32.15816943,31.93613243,32.3589859,0.42285347,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,32.79837,24.43307804,31.99603,26.22237753,1.993,0.58590816,3.434853817,0.170577319,-2.5515
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EQ7,EQ,Early,1,2,1,Target,41,32.69641113,32.55103683,32.39569473,32.54771423,32.39569473,32.69641113,0.3007164,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.19567,25.14049381,31.99603,26.22237753,1.993,0.44963955,2.108783874,0.213222208,-2.22957
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EQ8,EQ,Early,1,2,1,Target,41,32.00284576,32.21360397,32.05162048,32.08935674,32.00284576,32.21360397,0.210758209,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,32.72819,24.92400707,31.99603,26.22237753,1.993,0.61395706,2.448338325,0.250764796,-1.99559
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EW3,EW,Early,2,1,1,Target,117,31.91641617,32.0904541,31.9654808,31.99078369,31.91641617,32.0904541,0.174037933,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,31.99078,24.92159537,31.99603,26.22237753,1.993,1.00349967,2.452413802,0.409188558,-1.28916
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EW4,EW,Early,2,1,1,Target,117,32.02261734,32.48761368,32.3212204,32.27715047,32.02261734,32.48761368,0.464996338,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.27715,24.71233996,31.99603,26.22237753,1.993,0.82918691,2.833135032,0.292674687,-1.77263
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EW5,EW,Early,2,1,1,Target,117,33.15146637,33.04536438,32.7242775,32.97370275,32.7242775,33.15146637,0.427188873,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.9737,25.76160569,31.99603,26.22237753,1.993,0.52130744,1.374056303,0.379393072,-1.39823
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EW6,EW,Early,2,2,1,Target,41,31.02464294,31.41745567,31.3861866,31.27609507,31.02464294,31.41745567,0.392812729,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,31.89873,24.47244454,31.99603,26.22237753,1.993,1.06697231,3.342856272,0.319179835,-1.64756
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EW7,EW,Early,2,2,1,Target,41,32.59467316,32.60866928,32.2456665,32.48300298,32.2456665,32.60866928,0.363002777,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.12967,24.71929254,31.99603,26.22237753,1.993,0.46985356,2.819583273,0.166639362,-2.5852
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,EW8,EW,Early,2,2,1,Target,41,31.02501297,31.24364853,31.31508446,31.19458199,31.02501297,31.31508446,0.290071487,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,31.8156,24.24981237,31.99603,26.22237753,1.993,1.1277422,3.897603549,0.289342459,-1.78915
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MQ3,MQ,Medium,3,1,1,Target,117,29.904562,30.10918045,30.32287979,30.11220741,29.904562,30.32287979,0.418317795,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,30.11221,25.89835631,31.99603,26.22237753,1.993,3.50843347,1.25039354,2.805863398,1.48844
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MQ4,MQ,Medium,3,1,1,Target,117,30.94125938,30.95086479,30.73322487,30.87511635,30.73322487,30.95086479,0.217639923,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,30.87512,25.52299076,31.99603,26.22237753,1.993,2.11034459,1.619837474,1.30281255,0.38163
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MQ5,MQ,Medium,3,1,1,Target,117,31.27940178,31.23721695,31.41853905,31.31171926,31.23721695,31.41853905,0.181322098,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,31.31172,25.12014958,31.99603,26.22237753,1.993,1.57766706,2.138579087,0.737717427,-0.43886
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MQ6,MQ,Medium,3,2,1,Target,41,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,29.29156,24.7511079,31.99603,26.22237753,1.993,6.06149036,2.75839205,2.197472386,1.13585
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MQ7,MQ,Medium,3,2,1,Target,41,29.70074081,29.95859337,30.15030098,29.93654505,29.70074081,30.15030098,0.449560165,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,30.53252,24.76237619,31.99603,26.22237753,1.993,2.65149,2.737039456,0.9687438,-0.04581
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MQ8,MQ,Medium,3,2,1,Target,41,29.92821884,29.76143456,29.51300812,29.7342205,29.51300812,29.92821884,0.415210724,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,30.32616,24.74641872,31.99603,26.22237753,1.993,3.04229519,2.767326721,1.099362488,0.13667
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MW3,MW,Medium,4,1,1,Target,117,31.79305077,31.79781532,31.74607849,31.77898153,31.74607849,31.79781532,0.051736832,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,31.77898,25.29041589,31.99603,26.22237753,1.993,1.15559477,1.901644416,0.607681836,-0.71861
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MW4,MW,Medium,4,1,1,Target,117,32.08586502,32.42887497,31.88261223,32.13245074,31.88261223,32.42887497,0.546262741,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.13245,24.63093495,31.99603,26.22237753,1.993,0.91311108,2.996737273,0.304701746,-1.71453
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MW5,MW,Medium,4,1,1,Target,117,33.22049332,32.95978546,33.12698746,33.10242208,32.95978546,33.22049332,0.260707855,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,33.10242,25.35247425,31.99603,26.22237753,1.993,0.47846143,1.821974874,0.262605943,-1.92903
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MW6,MW,Medium,4,2,1,Target,41,32.89741516,32.68775177,32.52897263,32.70471319,32.52897263,32.89741516,0.368442535,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.35579,24.95817119,31.99603,26.22237753,1.993,0.40413824,2.391327326,0.169001638,-2.56489
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MW7,MW,Medium,4,2,1,Target,41,33.70588684,34.00052643,33.7530899,33.81983439,33.70588684,34.00052643,0.294639587,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,34.49311,25.70248086,31.99603,26.22237753,1.993,0.18942095,1.431241263,0.132347323,-2.9176
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,MW8,MW,Medium,4,2,1,Target,41,32.76914215,32.9332962,32.28038788,32.66094208,32.28038788,32.9332962,0.652908325,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.31115,24.59640071,31.99603,26.22237753,1.993,0.41633983,3.068964925,0.135661319,-2.88192
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LQ3,LQ,Late,5,1,1,Target,117,32.14227676,31.77741432,32.06838989,31.99602699,31.77741432,32.14227676,0.364862442,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,31.99603,26.22237753,31.99603,26.22237753,1.993,1,1,1,0
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LQ4,LQ,Late,5,1,1,Target,117,33.70827484,34.02420425,33.26028824,33.66425578,33.26028824,34.02420425,0.763916016,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,33.66426,25.57874041,31.99603,26.22237753,1.993,0.32905719,1.558741203,0.211104443,-2.24397
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LQ5,LQ,Late,5,1,1,Target,117,33.69075775,33.29528809,33.03218079,33.33940887,33.03218079,33.69075775,0.658576965,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,33.33941,26.77102083,31.99603,26.22237753,1.993,0.40857428,0.684979114,0.596476989,-0.74546
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LQ6,LQ,Late,5,2,1,Target,41,32.81244278,32.70795822,32.78907394,32.76982498,32.70795822,32.81244278,0.104484558,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.4222,26.04640784,31.99603,26.22237753,1.993,0.38664618,1.129026689,0.34245973,-1.54599
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LQ7,LQ,Late,5,2,1,Target,41,33.07961273,33.07355499,33.42037582,33.19118118,33.07355499,33.42037582,0.346820831,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.85195,24.94259361,31.99603,26.22237753,1.993,0.29037572,2.417155701,0.120131159,-3.05732
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LQ8,LQ,Late,5,2,1,Target,41,33.15650558,33.1697464,33.88676834,33.40434011,33.15650558,33.88676834,0.730262756,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,34.06935,25.43112605,31.99603,26.22237753,1.993,0.25121839,1.725780648,0.145567972,-2.78024
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LW3,LW,Late,6,1,1,Target,117,31.94018936,32.04188919,32.25728607,32.07978821,31.94018936,32.25728607,0.31709671,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.07979,26.28165682,31.99603,26.22237753,1.993,0.94571953,0.959942941,0.985183063,-0.02154
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LW4,LW,Late,6,1,1,Target,117,NA,32.61062241,33.03787994,32.82425117,32.61062241,33.03787994,0.427257538,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.82425,25.96422745,31.99603,26.22237753,1.993,0.57589077,1.194862228,0.481972533,-1.05298
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LW5,LW,Late,6,1,1,Target,117,32.31230164,31.86345482,32.02320862,32.06632169,31.86345482,32.31230164,0.448846817,mq6,29.17470741,29.39232063,29.30766487,29.29156431,29.17470741,29.39232063,0.21761322,29.29156431,1.947,32.06632,25.98699537,31.99603,26.22237753,1.993,0.95424327,1.176247401,0.811260687,-0.30176
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LW6,LW,Late,6,2,1,Target,41,34.22609329,34.10928726,33.98849487,34.10795848,33.98849487,34.22609329,0.237598419,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,34.78697,25.6146532,31.99603,26.22237753,1.993,0.15573812,1.520610107,0.10241818,-3.28746
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LW7,LW,Late,6,2,1,Target,41,32.06354523,31.88721848,31.72697067,31.89257813,31.72697067,32.06354523,0.336574554,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,32.52749,25.49879946,31.99603,26.22237753,1.993,0.70179877,1.647088425,0.426084452,-1.23079
Cytochrome P450 305A1 (LOC100647578),Cyt305A1.647578,LW8,LW,Late,6,2,1,Target,41,32.97445679,32.73745346,32.95695877,32.88962301,32.73745346,32.97445679,0.237003326,mq6,28.67350388,28.82612991,28.65981102,28.71981494,28.65981102,28.82612991,0.166318893,29.29156431,1.947,33.54438,25.38152329,31.99603,26.22237753,1.993,0.35641677,1.78583768,0.199579602,-2.32496
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EQ3,EQ,Early,1,1,1,Target,118,25.46961594,25.38253212,25.40411949,25.41875585,25.38253212,25.46961594,0.087083817,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.87103,25.20347779,27.07046,26.22237753,1.993,0.58951283,2.019146866,0.291961345,-1.77615
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EQ4,EQ,Early,1,1,1,Target,118,25.33742332,25.35497665,25.37553406,25.35597801,25.33742332,25.37553406,0.038110733,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.8022,25.71323328,27.07046,26.22237753,1.993,0.61691703,1.420667419,0.434244509,-1.20342
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EQ5,EQ,Early,1,1,1,Target,118,25.29314041,25.44976616,25.74672127,25.49654261,25.29314041,25.74672127,0.453580856,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.95632,25.41256952,27.07046,26.22237753,1.993,0.55723949,1.748007975,0.318785442,-1.64934
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EQ6,EQ,Early,1,2,1,Target,66,27.36457634,27.47835159,27.45609474,27.43300756,27.36457634,27.47835159,0.113775253,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,27.43301,24.43307804,27.07046,26.22237753,1.993,0.78716501,3.434853817,0.229169873,-2.12551
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EQ7,EQ,Early,1,2,1,Target,66,27.7047863,27.74138069,27.56933212,27.67183304,27.56933212,27.74138069,0.172048569,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,27.67183,25.14049381,27.07046,26.22237753,1.993,0.67235566,2.108783874,0.318835735,-1.64911
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EQ8,EQ,Early,1,2,1,Target,66,27.69600868,27.7604332,27.56201935,27.67282041,27.56201935,27.7604332,0.198413849,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,27.67282,24.92400707,27.07046,26.22237753,1.993,0.67191758,2.448338325,0.274438206,-1.86545
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EW3,EW,Early,2,1,1,Target,118,24.12146378,24.12981796,24.10528183,24.11885452,24.10528183,24.12981796,0.024536133,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.44572,24.92159537,27.07046,26.22237753,1.993,1.51043592,2.452413802,0.615897659,-0.69924
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EW4,EW,Early,2,1,1,Target,118,23.9348774,23.87110901,23.89574432,23.90057691,23.87110901,23.9348774,0.063768387,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.20638,24.71233996,27.07046,26.22237753,1.993,1.76894921,2.833135032,0.624378714,-0.67951
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EW5,EW,Early,2,1,1,Target,118,24.0350647,24.05990028,24.38645363,24.16047287,24.0350647,24.38645363,0.351388931,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.49135,25.76160569,27.07046,26.22237753,1.993,1.46561552,1.374056303,1.066634255,0.09307
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EW6,EW,Early,2,2,1,Target,66,26.95765686,26.72694397,26.86885262,26.85115115,26.72694397,26.95765686,0.230712891,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.85115,24.47244454,27.07046,26.22237753,1.993,1.15577406,3.342856272,0.345744466,-1.53222
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EW7,EW,Early,2,2,1,Target,66,26.93841934,26.97268867,26.76980972,26.89363925,26.76980972,26.97268867,0.202878952,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.89364,24.71929254,27.07046,26.22237753,1.993,1.12380876,2.819583273,0.398572644,-1.32709
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,EW8,EW,Early,2,2,1,Target,66,26.42745018,26.45697975,26.38986206,26.424764,26.38986206,26.45697975,0.067117691,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.42476,24.24981237,27.07046,26.22237753,1.993,1.5314754,3.897603549,0.392927443,-1.34767
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MQ3,MQ,Medium,3,1,1,Target,118,24.0477581,24.0526638,24.07753181,24.05931791,24.0477581,24.07753181,0.029773712,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.38044,25.89835631,27.07046,26.22237753,1.993,1.57694649,1.25039354,1.261160134,0.33475
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MQ4,MQ,Medium,3,1,1,Target,118,24.85794067,25.06275558,25.04525375,24.98865,24.85794067,25.06275558,0.204814911,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.39943,25.52299076,27.07046,26.22237753,1.993,0.80480765,1.619837474,0.496844688,-1.00913
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MQ5,MQ,Medium,3,1,1,Target,118,23.5761795,23.44093132,23.89781189,23.63830757,23.44093132,23.89781189,0.456880569,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,25.91881,25.12014958,27.07046,26.22237753,1.993,2.13873371,2.138579087,1.000072303,0.0001
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MQ6,MQ,Medium,3,2,1,Target,66,25.76529121,25.8136692,25.73704338,25.77200127,25.73704338,25.8136692,0.076625824,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,25.772,24.7511079,27.07046,26.22237753,1.993,2.35637646,2.75839205,0.85425727,-0.22726
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MQ7,MQ,Medium,3,2,1,Target,66,25.89699554,26.09606171,25.96486855,25.98597527,25.89699554,26.09606171,0.199066162,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,25.98598,24.76237619,27.07046,26.22237753,1.993,2.0459848,2.737039456,0.747517466,-0.41982
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MQ8,MQ,Medium,3,2,1,Target,66,26.32817459,26.18644905,26.10565758,26.20676041,26.10565758,26.32817459,0.222517014,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.20676,24.74641872,27.07046,26.22237753,1.993,1.76850994,2.767326721,0.639067994,-0.64596
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MW3,MW,Medium,4,1,1,Target,118,24.07217216,24.16810799,23.98228455,24.07418823,23.98228455,24.16810799,0.185823441,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.39674,25.29041589,27.07046,26.22237753,1.993,1.56006481,1.901644416,0.820376718,-0.28564
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MW4,MW,Medium,4,1,1,Target,118,23.79571342,23.73759651,23.68937302,23.74089432,23.68937302,23.79571342,0.106340408,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.0313,24.63093495,27.07046,26.22237753,1.993,1.98568224,2.996737273,0.662614722,-0.59376
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MW5,MW,Medium,4,1,1,Target,118,23.34509277,23.38544273,23.70349884,23.47801145,23.34509277,23.70349884,0.358406067,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,25.74305,25.35247425,27.07046,26.22237753,1.993,2.40183936,1.821974874,1.318261519,0.39864
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MW6,MW,Medium,4,2,1,Target,66,25.56773567,25.13899422,25.4583931,25.38837433,25.13899422,25.56773567,0.428741455,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,25.38837,24.95817119,27.07046,26.22237753,1.993,3.03545094,2.391327326,1.269358197,0.3441
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MW7,MW,Medium,4,2,1,Target,66,27.27843094,27.02513695,27.09591103,27.13315964,27.02513695,27.27843094,0.253293991,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,27.13316,25.70248086,27.07046,26.22237753,1.993,0.95945912,1.431241263,0.670368542,-0.57697
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,MW8,MW,Medium,4,2,1,Target,66,28.33098221,28.21114922,28.11151886,28.21788343,28.11151886,28.33098221,0.219463348,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,28.21788,24.59640071,27.07046,26.22237753,1.993,0.46887471,3.068964925,0.152779429,-2.71048
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LQ3,LQ,Late,5,1,1,Target,118,24.47632408,24.77762032,24.81194496,24.68862979,24.47632408,24.81194496,0.33562088,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.07046,26.22237753,27.07046,26.22237753,1.993,1,1,1,0
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LQ4,LQ,Late,5,1,1,Target,118,24.58530235,24.54750443,24.66827965,24.60036214,24.54750443,24.66827965,0.120775223,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,26.97368,25.57874041,27.07046,26.22237753,1.993,1.0659723,1.558741203,0.683867404,-0.54821
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LQ5,LQ,Late,5,1,1,Target,118,23.46407318,23.57530403,23.84346008,23.62761243,23.46407318,23.84346008,0.379386902,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,25.90709,26.77102083,27.07046,26.22237753,1.993,2.155354,0.684979114,3.146598128,1.65379
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LQ6,LQ,Late,5,2,1,Target,66,28.35870361,28.18958855,28.35185051,28.30004756,28.18958855,28.35870361,0.169115067,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,28.30005,26.04640784,27.07046,26.22237753,1.993,0.44412162,1.129026689,0.393366804,-1.34605
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LQ7,LQ,Late,5,2,1,Target,66,26.67309189,26.3596096,26.3021431,26.4449482,26.3021431,26.67309189,0.370948792,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.44495,24.94259361,27.07046,26.22237753,1.993,1.51120575,2.417155701,0.625200003,-0.67761
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LQ8,LQ,Late,5,2,1,Target,66,27.95396042,28.11182785,27.62235641,27.89604823,27.62235641,28.11182785,0.489471436,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,27.89605,25.43112605,27.07046,26.22237753,1.993,0.57985691,1.725780648,0.335996876,-1.57348
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LW3,LW,Late,6,1,1,Target,118,26.80009842,26.79597855,26.84261131,26.81289609,26.79597855,26.84261131,0.046632767,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,29.39967,26.28165682,27.07046,26.22237753,1.993,0.21491232,0.959942941,0.223880303,-2.1592
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LW4,LW,Late,6,1,1,Target,118,25.155056,25.16475677,25.28145599,25.20042292,25.155056,25.28145599,0.126399994,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.63163,25.96422745,27.07046,26.22237753,1.993,0.69043645,1.194862228,0.577837707,-0.79126
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LW5,LW,Late,6,1,1,Target,118,24.56772423,24.88405609,25.17761612,24.87646548,24.56772423,25.17761612,0.609891891,lw6,24.58464432,24.56979561,24.61177063,24.58873685,24.56979561,24.61177063,0.041975021,26.96093432,1.935,27.27642,25.98699537,27.07046,26.22237753,1.993,0.87288295,1.176247401,0.742091286,-0.43033
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LW6,LW,Late,6,2,1,Target,66,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.96093,25.6146532,27.07046,26.22237753,1.993,1.07497953,1.520610107,0.706939621,-0.50034
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LW7,LW,Late,6,2,1,Target,66,27.75458717,27.63883018,27.66638184,27.68659973,27.63883018,27.75458717,0.115756989,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,27.6866,25.49879946,27.07046,26.22237753,1.993,0.66583364,1.647088425,0.404248872,-1.30668
Cytochrome P450 6A1-like (LOC100649469),Cyt6A1.649469,LW8,LW,Late,6,2,1,Target,66,27.21907043,26.94024277,26.83659554,26.99863625,26.83659554,27.21907043,0.382474899,lw6,26.78355598,27.14247704,26.95676994,26.96093432,26.78355598,27.14247704,0.358921051,26.96093432,1.935,26.99864,25.38152329,27.07046,26.22237753,1.993,1.04855634,1.78583768,0.587150979,-0.7682
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EQ3,EQ,Early,1,1,1,Target,127,30.06845665,30.18653488,30.10205269,30.11901474,30.06845665,30.18653488,0.118078232,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,30.11901,25.20347779,31.76822,26.22237753,1.993,2.98546401,2.019146866,1.478576946,0.56421
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EQ4,EQ,Early,1,1,1,Target,127,29.99512672,30.43705177,30.17686653,30.20301501,29.99512672,30.43705177,0.441925049,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,30.20302,25.71323328,31.76822,26.22237753,1.993,2.82369391,1.420667419,1.987582647,0.99101
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EQ5,EQ,Early,1,1,1,Target,127,32.554039,32.05929565,31.98326492,32.19886653,31.98326492,32.554039,0.570774078,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,32.19887,25.41256952,31.76822,26.22237753,1.993,0.75155612,1.748007975,0.429950052,-1.21776
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EQ6,EQ,Early,1,2,1,Target,40,29.02910423,28.84080887,28.81695366,28.89562225,28.81695366,29.02910423,0.212150574,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,28.91685,24.43307804,31.76822,26.22237753,1.993,6.62620598,3.434853817,1.929108582,0.94793
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EQ7,EQ,Early,1,2,1,Target,40,29.36522865,29.60648727,29.33177376,29.43449656,29.33177376,29.60648727,0.274713516,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,29.45612,25.14049381,31.76822,26.22237753,1.993,4.63384314,2.108783874,2.197400692,1.1358
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EQ8,EQ,Early,1,2,1,Target,40,28.89372635,28.87622261,28.80748177,28.85914358,28.80748177,28.89372635,0.086244583,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,28.88035,24.92400707,31.76822,26.22237753,1.993,6.78858785,2.448338325,2.772732749,1.47131
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EW3,EW,Early,2,1,1,Target,127,32.1450386,31.83924294,31.80376244,31.92934799,31.80376244,32.1450386,0.341276169,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.92935,24.92159537,31.76822,26.22237753,1.993,0.89864805,2.452413802,0.366434102,-1.44837
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EW4,EW,Early,2,1,1,Target,127,29.28791618,29.15317726,29.23928642,29.22679329,29.15317726,29.28791618,0.134738922,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,29.22679,24.71233996,31.76822,26.22237753,1.993,5.39503667,2.833135032,1.904263865,0.92923
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EW5,EW,Early,2,1,1,Target,127,30.02805138,29.93121529,29.96840096,29.97588921,29.93121529,30.02805138,0.09683609,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,29.97589,25.76160569,31.76822,26.22237753,1.993,3.28273366,1.374056303,2.389082347,1.25646
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EW6,EW,Early,2,2,1,Target,40,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,28.98778,24.47244454,31.76822,26.22237753,1.993,6.32174525,3.342856272,1.891120866,0.91924
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EW7,EW,Early,2,2,1,Target,40,28.45527649,28.42075157,28.33744431,28.40449079,28.33744431,28.45527649,0.117832184,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,28.42536,24.71929254,31.76822,26.22237753,1.993,9.17967432,2.819583273,3.255684768,1.70296
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,EW8,EW,Early,2,2,1,Target,40,28.14135361,28.11037445,28.13684273,28.1295236,28.11037445,28.14135361,0.030979156,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,28.15019,24.24981237,31.76822,26.22237753,1.993,11.01749946,3.897603549,2.826736818,1.49914
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MQ3,MQ,Medium,3,1,1,Target,127,28.92655373,28.35752678,28.72523689,28.66977247,28.35752678,28.92655373,0.569026947,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,28.66977,25.89835631,31.76822,26.22237753,1.993,7.80603729,1.25039354,6.242864374,2.64221
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MQ4,MQ,Medium,3,1,1,Target,127,31.06816292,30.91874313,30.95375252,30.98021952,30.91874313,31.06816292,0.149419785,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,30.98022,25.52299076,31.76822,26.22237753,1.993,1.68640884,1.619837474,1.041097561,0.05811
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MQ5,MQ,Medium,3,1,1,Target,127,31.90233994,31.94924545,31.76890564,31.87349701,31.76890564,31.94924545,0.180339813,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.8735,25.12014958,31.76822,26.22237753,1.993,0.93255863,2.138579087,0.436064597,-1.19739
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MQ6,MQ,Medium,3,2,1,Target,40,28.43721199,28.35013008,28.5477581,28.44503339,28.35013008,28.5477581,0.197628021,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,28.46593,24.7511079,31.76822,26.22237753,1.993,8.93596355,2.75839205,3.239555285,1.6958
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MQ7,MQ,Medium,3,2,1,Target,40,33.02201462,33.39398956,32.97689056,33.13096491,32.97689056,33.39398956,0.417098999,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,33.15531,24.76237619,31.76822,26.22237753,1.993,0.39854904,2.737039456,0.145613188,-2.77979
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MQ8,MQ,Medium,3,2,1,Target,40,32.90652466,32.61271286,32.89046097,32.80323283,32.61271286,32.90652466,0.293811798,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,32.82734,24.74641872,31.76822,26.22237753,1.993,0.49538897,2.767326721,0.179013548,-2.48186
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MW3,MW,Medium,4,1,1,Target,127,31.09654617,31.3149128,30.83062744,31.08069547,30.83062744,31.3149128,0.484285355,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.0807,25.29041589,31.76822,26.22237753,1.993,1.57769567,1.901644416,0.829648094,-0.26943
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MW4,MW,Medium,4,1,1,Target,127,31.75652122,31.53681946,31.61574364,31.63636144,31.53681946,31.75652122,0.219701767,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.63636,24.63093495,31.76822,26.22237753,1.993,1.09138315,2.996737273,0.364190468,-1.45723
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MW5,MW,Medium,4,1,1,Target,127,30.73868942,30.9980011,30.62554932,30.78741328,30.62554932,30.9980011,0.372451782,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,30.78741,25.35247425,31.76822,26.22237753,1.993,1.91644338,1.821974874,1.05184951,0.07293
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MW6,MW,Medium,4,2,1,Target,40,31.44709778,31.04755402,31.62453842,31.37306341,31.04755402,31.62453842,0.576984406,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,31.39612,24.95817119,31.76822,26.22237753,1.993,1.27989492,2.391327326,0.535223643,-0.90179
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MW7,MW,Medium,4,2,1,Target,40,31.41949654,31.30549812,31.26549911,31.33016459,31.26549911,31.41949654,0.153997421,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,31.35318,25.70248086,31.76822,26.22237753,1.993,1.31685922,1.431241263,0.920081928,-0.12017
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,MW8,MW,Medium,4,2,1,Target,40,29.23556137,29.25228119,29.09303856,29.19362704,29.09303856,29.25228119,0.15924263,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,29.21508,24.59640071,31.76822,26.22237753,1.993,5.43711949,3.068964925,1.771646017,0.82509
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LQ3,LQ,Late,5,1,1,Target,127,31.80318642,31.9601078,31.54135132,31.76821518,31.54135132,31.9601078,0.418756485,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.76822,26.22237753,31.76822,26.22237753,1.993,1,1,1,0
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LQ4,LQ,Late,5,1,1,Target,127,31.9319458,31.98293877,31.82653427,31.91380628,31.82653427,31.98293877,0.156404495,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.91381,25.57874041,31.76822,26.22237753,1.993,0.9079586,1.558741203,0.582494769,-0.77968
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LQ5,LQ,Late,5,1,1,Target,127,31.46861076,32.18412399,31.88259888,31.84511121,31.46861076,32.18412399,0.715513229,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,31.84511,26.77102083,31.76822,26.22237753,1.993,0.95028085,0.684979114,1.38731362,0.47229
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LQ6,LQ,Late,5,2,1,Target,40,30.74359894,30.9414978,30.80792236,30.83100637,30.74359894,30.9414978,0.197898865,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,30.85366,26.04640784,31.76822,26.22237753,1.993,1.83406744,1.129026689,1.624467745,0.69997
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LQ7,LQ,Late,5,2,1,Target,40,30.64970016,30.70199203,30.81061172,30.72076797,30.64970016,30.81061172,0.16091156,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,30.74334,24.94259361,31.76822,26.22237753,1.993,1.97328631,2.417155701,0.816367067,-0.29271
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LQ8,LQ,Late,5,2,1,Target,40,32.41772461,31.96828651,32.13552094,32.17384402,31.96828651,32.41772461,0.449438095,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,32.19748,25.43112605,31.76822,26.22237753,1.993,0.75224547,1.725780648,0.435887071,-1.19797
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LW3,LW,Late,6,1,1,Target,127,35.91869736,35.18061829,35.85519028,35.65150197,35.18061829,35.91869736,0.738079071,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,35.6515,26.28165682,31.76822,26.22237753,1.993,0.07612271,0.959942941,0.079299203,-3.65655
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LW4,LW,Late,6,1,1,Target,127,32.47270584,32.60563278,32.24192429,32.44008764,32.24192429,32.60563278,0.363708496,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,32.44009,25.96422745,31.76822,26.22237753,1.993,0.64044756,1.194862228,0.536001177,-0.89969
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LW5,LW,Late,6,1,1,Target,127,36.58141708,35.57318115,NA,36.07729912,35.57318115,36.58141708,1.008235931,ew6,29.00437927,28.9483242,29.01062965,28.98777771,28.9483242,29.01062965,0.06230545,28.98777771,1.941,36.0773,25.98699537,31.76822,26.22237753,1.993,0.05739496,1.176247401,0.048794977,-4.35712
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LW6,LW,Late,6,2,1,Target,40,33.07368088,33.03202057,32.76538849,32.95702998,32.76538849,33.07368088,0.308292389,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,32.98125,25.6146532,31.76822,26.22237753,1.993,0.4473179,1.520610107,0.294170019,-1.76528
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LW7,LW,Late,6,2,1,Target,40,33.99365616,33.7709465,33.69781113,33.8208046,33.69781113,33.99365616,0.295845032,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,33.84565,25.49879946,31.76822,26.22237753,1.993,0.2521415,1.647088425,0.15308316,-2.70761
Cytochrome P450 6k1-like (LOC100642936),Cyt6k1.642936,LW8,LW,Late,6,2,1,Target,40,32.45653534,32.54871368,32.61598969,32.5404129,32.45653534,32.61598969,0.159454346,ew6,28.9426384,28.85450935,29.10233498,28.96649424,28.85450935,29.10233498,0.247825623,28.98777771,1.941,32.56432,25.38152329,31.76822,26.22237753,1.993,0.58979462,1.78583768,0.330262166,-1.59832
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EQ3,EQ,Early,1,1,1,Target,112,40,40,40,40,40,40,0,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,42.12005,25.20347779,34.01324,26.22237753,1.993,0.00410031,2.019146866,0.002030715,-8.9438
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EQ4,EQ,Early,1,1,1,Target,112,NA,37.29777527,37.15759277,37.22768402,37.15759277,37.29777527,0.140182495,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,39.2008,25.71323328,34.01324,26.22237753,1.993,0.02967824,1.420667419,0.020890351,-5.58102
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EQ5,EQ,Early,1,1,1,Target,112,40,NA,40,40,40,40,0,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,42.12005,25.41256952,34.01324,26.22237753,1.993,0.00410031,1.748007975,0.002345706,-8.73576
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EQ6,EQ,Early,1,2,1,Target,60,39.2870636,40,40,39.76235453,39.2870636,40,0.712936401,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,39.76235,24.43307804,34.01324,26.22237753,1.993,0.0202805,3.434853817,0.005904328,-7.40401
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EQ7,EQ,Early,1,2,1,Target,60,37.02188492,37.48595428,NA,37.2539196,37.02188492,37.48595428,0.464069366,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,37.25392,25.14049381,34.01324,26.22237753,1.993,0.11110359,2.108783874,0.052686095,-4.24643
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EQ8,EQ,Early,1,2,1,Target,60,40,40,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,24.92400707,34.01324,26.22237753,1.993,0.01726236,2.448338325,0.007050643,-7.14803
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EW3,EW,Early,2,1,1,Target,112,36.46503448,NA,37.11994171,36.7924881,36.46503448,37.11994171,0.654907227,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,38.74254,24.92159537,34.01324,26.22237753,1.993,0.04049308,2.452413802,0.016511519,-5.92038
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EW4,EW,Early,2,1,1,Target,112,36.37365341,35.02533722,NA,35.69949532,35.02533722,36.37365341,1.348316193,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,37.59161,24.71233996,34.01324,26.22237753,1.993,0.08836671,2.833135032,0.031190433,-5.00275
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EW5,EW,Early,2,1,1,Target,112,NA,36.96751404,36.30041122,36.63396263,36.30041122,36.96751404,0.667102814,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,38.57561,25.76160569,34.01324,26.22237753,1.993,0.04534562,1.374056303,0.033001279,-4.92133
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EW6,EW,Early,2,2,1,Target,60,40,40,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,24.47244454,34.01324,26.22237753,1.993,0.01726236,3.342856272,0.005163955,-7.59731
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EW7,EW,Early,2,2,1,Target,60,38.10348129,NA,37.28739929,37.69544029,37.28739929,38.10348129,0.816082001,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,37.69544,24.71929254,34.01324,26.22237753,1.993,0.08235978,2.819583273,0.029209911,-5.0974
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,EW8,EW,Early,2,2,1,Target,60,40,38.64887238,NA,39.32443619,38.64887238,40,1.351127625,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,39.32444,24.24981237,34.01324,26.22237753,1.993,0.02729172,3.897603549,0.007002179,-7.15798
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MQ3,MQ,Medium,3,1,1,Target,112,37.10730362,NA,35.78305817,36.44518089,35.78305817,37.10730362,1.324245453,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,38.37682,25.89835631,34.01324,26.22237753,1.993,0.05188854,1.25039354,0.04149777,-4.59082
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MQ4,MQ,Medium,3,1,1,Target,112,40,NA,40,40,40,40,0,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,42.12005,25.52299076,34.01324,26.22237753,1.993,0.00410031,1.619837474,0.002531311,-8.6259
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MQ5,MQ,Medium,3,1,1,Target,112,40,40,NA,40,40,40,0,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,42.12005,25.12014958,34.01324,26.22237753,1.993,0.00410031,2.138579087,0.001917307,-9.0267
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MQ6,MQ,Medium,3,2,1,Target,60,40,NA,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,24.7511079,34.01324,26.22237753,1.993,0.01726236,2.75839205,0.006258124,-7.32005
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MQ7,MQ,Medium,3,2,1,Target,60,NA,40,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,24.76237619,34.01324,26.22237753,1.993,0.01726236,2.737039456,0.006306946,-7.30884
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MQ8,MQ,Medium,3,2,1,Target,60,40,40,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,24.74641872,34.01324,26.22237753,1.993,0.01726236,2.767326721,0.006237919,-7.32472
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MW3,MW,Medium,4,1,1,Target,112,36.152668,36.2616272,NA,36.2071476,36.152668,36.2616272,0.108959198,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,38.12617,25.29041589,34.01324,26.22237753,1.993,0.06150059,1.901644416,0.032340739,-4.9505
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MW4,MW,Medium,4,1,1,Target,112,40,40,40,40,40,40,0,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,42.12005,24.63093495,34.01324,26.22237753,1.993,0.00410031,2.996737273,0.001368259,-9.51344
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MW5,MW,Medium,4,1,1,Target,112,34.46023941,34.62921906,33.92749405,34.33898417,33.92749405,34.62921906,0.701725006,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,36.15899,25.35247425,34.01324,26.22237753,1.993,0.23342473,1.821974874,0.128116327,-2.96447
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MW6,MW,Medium,4,2,1,Target,60,40,40,NA,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,24.95817119,34.01324,26.22237753,1.993,0.01726236,2.391327326,0.007218735,-7.11404
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MW7,MW,Medium,4,2,1,Target,60,40,40,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,25.70248086,34.01324,26.22237753,1.993,0.01726236,1.431241263,0.012061111,-6.37349
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,MW8,MW,Medium,4,2,1,Target,60,35.31082153,34.61149979,34.38155365,34.76795832,34.38155365,35.31082153,0.929267883,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,34.76796,24.59640071,34.01324,26.22237753,1.993,0.59945932,3.068964925,0.19532948,-2.35602
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LQ3,LQ,Late,5,1,1,Target,112,32.04373169,32.46821594,32.39174271,32.30123011,32.04373169,32.46821594,0.424484253,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,34.01324,26.22237753,34.01324,26.22237753,1.993,1,1,1,0
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LQ4,LQ,Late,5,1,1,Target,112,30.63370895,30.24212646,30.30926132,30.39503225,30.24212646,30.63370895,0.391582489,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,32.00601,25.57874041,34.01324,26.22237753,1.993,3.89996827,1.558741203,2.501998575,1.32308
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LQ5,LQ,Late,5,1,1,Target,112,30.3167305,30.52750778,30.07881546,30.30768458,30.07881546,30.52750778,0.448692322,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,31.91403,26.77102083,34.01324,26.22237753,1.993,4.1509285,0.684979114,6.059934416,2.5993
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LQ6,LQ,Late,5,2,1,Target,60,34.05728531,33.41775131,33.08649826,33.52051163,33.08649826,34.05728531,0.970787048,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,33.52051,26.04640784,34.01324,26.22237753,1.993,1.39665982,1.129026689,1.237047661,0.3069
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LQ7,LQ,Late,5,2,1,Target,60,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,28.86445,24.94259361,34.01324,26.22237753,1.993,32.82044961,2.417155701,13.57812804,3.76321
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LQ8,LQ,Late,5,2,1,Target,60,34.15986252,34.25359726,34.08116531,34.16487503,34.08116531,34.25359726,0.172431946,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,34.16488,25.43112605,34.01324,26.22237753,1.993,0.90229257,1.725780648,0.522831549,-0.93558
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LW3,LW,Late,6,1,1,Target,112,34.02347183,33.64465714,33.92927933,33.86580276,33.64465714,34.02347183,0.378814697,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,35.66073,26.28165682,34.01324,26.22237753,1.993,0.32724106,0.959942941,0.340896362,-1.55259
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LW4,LW,Late,6,1,1,Target,112,35.0999527,34.71051788,35.37016678,35.06021245,34.71051788,35.37016678,0.659648895,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,36.91845,25.96422745,34.01324,26.22237753,1.993,0.13948041,1.194862228,0.116733469,-3.09871
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LW5,LW,Late,6,1,1,Target,112,34.29133224,34.13283539,34.31762314,34.24726359,34.13283539,34.31762314,0.18478775,lq7,26.98518944,27.57252693,27.67707634,27.41159757,26.98518944,27.67707634,0.691886902,28.86444664,1.97,36.06241,25.98699537,34.01324,26.22237753,1.993,0.24922234,1.176247401,0.211879183,-2.23869
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LW6,LW,Late,6,2,1,Target,60,37.95858002,38.43982697,NA,38.19920349,37.95858002,38.43982697,0.481246948,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,38.1992,25.6146532,34.01324,26.22237753,1.993,0.05852938,1.520610107,0.03849072,-4.69935
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LW7,LW,Late,6,2,1,Target,60,36.98826218,NA,38.31073761,37.64949989,36.98826218,38.31073761,1.322475433,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,37.6495,25.49879946,34.01324,26.22237753,1.993,0.08496558,1.647088425,0.05158532,-4.2769
Cytochrome P450 6K1-like (LOC100648995),Cyt6k1.648995,LW8,LW,Late,6,2,1,Target,60,40,40,40,40,40,40,0,lq7,29.01778603,28.79082108,28.78473282,28.86444664,28.78473282,29.01778603,0.233053207,28.86444664,1.97,40,25.38152329,34.01324,26.22237753,1.993,0.01726236,1.78583768,0.009666253,-6.69283
Hexamerin,Hex,EQ3,EQ,Early,1,1,1,Target,113,25.0855732,25.96491814,25.9702816,25.67359098,25.0855732,25.9702816,0.884708405,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,25.67359,25.20347779,27.25065,26.22237753,1.993,2.68354564,2.019146866,1.329049256,0.41039
Hexamerin,Hex,EQ4,EQ,Early,1,1,1,Target,113,27.45504951,27.25394058,27.25693512,27.32197507,27.25394058,27.45504951,0.201108932,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,27.32198,25.71323328,27.25065,26.22237753,1.993,0.95633402,1.420667419,0.673158269,-0.57098
Hexamerin,Hex,EQ5,EQ,Early,1,1,1,Target,113,24.46604347,24.36461449,24.28368187,24.37144661,24.28368187,24.46604347,0.182361603,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,24.37145,25.41256952,27.25065,26.22237753,1.993,6.06297815,1.748007975,3.468507143,1.79431
Hexamerin,Hex,EQ6,EQ,Early,1,2,1,Target,27,26.71335602,26.29445648,26.23688698,26.41489983,26.23688698,26.71335602,0.47646904,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,26.52815,24.43307804,27.25065,26.22237753,1.993,1.57182044,3.434853817,0.457609121,-1.12781
Hexamerin,Hex,EQ7,EQ,Early,1,2,1,Target,27,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.66038,25.14049381,27.25065,26.22237753,1.993,5.05988923,2.108783874,2.399434713,1.26269
Hexamerin,Hex,EQ8,EQ,Early,1,2,1,Target,27,25.48739815,25.53453445,25.50214195,25.50802485,25.48739815,25.53453445,0.047136307,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,25.61739,24.92400707,27.25065,26.22237753,1.993,2.77962542,2.448338325,1.135310995,0.18309
Hexamerin,Hex,EW3,EW,Early,2,1,1,Target,113,24.96179581,24.96540833,24.92566872,24.95095762,24.92566872,24.96540833,0.039739609,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,24.95096,24.92159537,27.25065,26.22237753,1.993,4.21842954,2.452413802,1.720113278,0.7825
Hexamerin,Hex,EW4,EW,Early,2,1,1,Target,113,24.3417263,24.10678291,24.03609085,24.16153336,24.03609085,24.3417263,0.305635452,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,24.16153,24.71233996,27.25065,26.22237753,1.993,6.91431497,2.833135032,2.440517268,1.28719
Hexamerin,Hex,EW5,EW,Early,2,1,1,Target,113,25.58447075,25.77022552,25.71416473,25.68962034,25.58447075,25.77022552,0.185754776,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,25.68962,25.76160569,27.25065,26.22237753,1.993,2.6567552,1.374056303,1.933512618,0.95122
Hexamerin,Hex,EW6,EW,Early,2,2,1,Target,27,24.51926422,24.30394554,24.38354111,24.40225029,24.30394554,24.51926422,0.21531868,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.50688,24.47244454,27.25065,26.22237753,1.993,5.57019791,3.342856272,1.666298954,0.73665
Hexamerin,Hex,EW7,EW,Early,2,2,1,Target,27,25.34914017,25.23212814,25.16313934,25.24813588,25.16313934,25.34914017,0.186000824,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,25.35639,24.71929254,27.25065,26.22237753,1.993,3.2729383,2.819583273,1.160787955,0.2151
Hexamerin,Hex,EW8,EW,Early,2,2,1,Target,27,24.4489193,24.33283424,24.44674492,24.40949949,24.33283424,24.4489193,0.116085052,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.51416,24.24981237,27.25065,26.22237753,1.993,5.54487223,3.897603549,1.422636285,0.50857
Hexamerin,Hex,MQ3,MQ,Medium,3,1,1,Target,113,27.13161469,26.96114349,26.98467445,27.02581088,26.96114349,27.13161469,0.170471191,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,27.02581,25.89835631,27.25065,26.22237753,1.993,1.1511167,1.25039354,0.920603527,-0.11935
Hexamerin,Hex,MQ4,MQ,Medium,3,1,1,Target,113,27.11944389,27.1601181,26.79416275,27.02457492,26.79416275,27.1601181,0.365955353,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,27.02457,25.52299076,27.25065,26.22237753,1.993,1.15200759,1.619837474,0.711187148,-0.4917
Hexamerin,Hex,MQ5,MQ,Medium,3,1,1,Target,113,26.82928085,26.69985008,26.75338364,26.76083819,26.69985008,26.82928085,0.129430771,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,26.76084,25.12014958,27.25065,26.22237753,1.993,1.35878254,2.138579087,0.635366981,-0.65434
Hexamerin,Hex,MQ6,MQ,Medium,3,2,1,Target,27,24.83633423,24.80494881,24.59872246,24.7466685,24.59872246,24.83633423,0.237611771,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.85277,24.7511079,27.25065,26.22237753,1.993,4.48582158,2.75839205,1.626245109,0.70154
Hexamerin,Hex,MQ7,MQ,Medium,3,2,1,Target,27,24.81648445,24.71940041,24.66418648,24.73335711,24.66418648,24.81648445,0.152297974,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.8394,24.76237619,27.25065,26.22237753,1.993,4.52351567,2.737039456,1.652703859,0.72483
Hexamerin,Hex,MQ8,MQ,Medium,3,2,1,Target,27,24.93113708,24.96587563,24.97495651,24.95732307,24.93113708,24.97495651,0.043819427,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,25.06433,24.74641872,27.25065,26.22237753,1.993,3.92945023,2.767326721,1.419944453,0.50583
Hexamerin,Hex,MW3,MW,Medium,4,1,1,Target,113,24.98783302,24.88887215,24.9982338,24.95831299,24.88887215,24.9982338,0.109361649,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,24.95831,25.29041589,27.25065,26.22237753,1.993,4.19905249,1.901644416,2.208116541,1.14282
Hexamerin,Hex,MW4,MW,Medium,4,1,1,Target,113,25.63407898,25.534832,25.47891426,25.54927508,25.47891426,25.63407898,0.155164719,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,25.54928,24.63093495,27.25065,26.22237753,1.993,2.9007027,2.996737273,0.967953623,-0.04699
Hexamerin,Hex,MW5,MW,Medium,4,1,1,Target,113,26.20692444,26.24088287,26.23793411,26.22858047,26.20692444,26.24088287,0.033958435,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,26.22858,25.35247425,27.25065,26.22237753,1.993,1.89600618,1.821974874,1.040632454,0.05746
Hexamerin,Hex,MW6,MW,Medium,4,2,1,Target,27,27.62569618,27.37350273,27.34725571,27.44881821,27.34725571,27.62569618,0.278440475,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,27.56651,24.95817119,27.25065,26.22237753,1.993,0.82060816,2.391327326,0.343160114,-1.54305
Hexamerin,Hex,MW7,MW,Medium,4,2,1,Target,27,27.12716866,26.77307892,26.80742836,26.90255864,26.77307892,27.12716866,0.354089737,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,27.0179,25.70248086,27.25065,26.22237753,1.993,1.15682747,1.431241263,0.808268666,-0.30709
Hexamerin,Hex,MW8,MW,Medium,4,2,1,Target,27,25.17949486,25.13293266,25.04796982,25.12013245,25.04796982,25.17949486,0.13152504,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,25.22784,24.59640071,27.25065,26.22237753,1.993,3.54718363,3.068964925,1.155824101,0.20892
Hexamerin,Hex,LQ3,LQ,Late,5,1,1,Target,113,27.26163483,27.45355225,27.03674889,27.25064532,27.03674889,27.45355225,0.41680336,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,27.25065,26.22237753,27.25065,26.22237753,1.993,1,1,1,0
Hexamerin,Hex,LQ4,LQ,Late,5,1,1,Target,113,25.1722641,24.90525055,24.96551514,25.01434326,24.90525055,25.1722641,0.26701355,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,25.01434,25.57874041,27.25065,26.22237753,1.993,4.05433793,1.558741203,2.601033398,1.37908
Hexamerin,Hex,LQ5,LQ,Late,5,1,1,Target,113,27.88399506,27.95222092,27.97529793,27.9371713,27.88399506,27.97529793,0.091302872,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,27.93717,26.77102083,27.25065,26.22237753,1.993,0.65068931,0.684979114,0.949940372,-0.07409
Hexamerin,Hex,LQ6,LQ,Late,5,2,1,Target,27,26.61481285,26.80630302,26.79803276,26.73971621,26.61481285,26.80630302,0.191490173,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,26.85436,26.04640784,27.25065,26.22237753,1.993,1.28152103,1.129026689,1.135067085,0.18278
Hexamerin,Hex,LQ7,LQ,Late,5,2,1,Target,27,26.81216621,26.53416061,26.68416405,26.67683029,26.53416061,26.81216621,0.2780056,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,26.79121,24.94259361,27.25065,26.22237753,1.993,1.33319613,2.417155701,0.551555751,-0.85842
Hexamerin,Hex,LQ8,LQ,Late,5,2,1,Target,27,24.65836716,24.81998634,24.78478241,24.75437864,24.65836716,24.81998634,0.161619186,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.86051,25.43112605,27.25065,26.22237753,1.993,4.46413246,1.725780648,2.586732253,1.37113
Hexamerin,Hex,LW3,LW,Late,6,1,1,Target,113,25.03248215,24.80679131,24.77005386,24.86977577,24.77005386,25.03248215,0.262428284,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,24.86978,26.28165682,27.25065,26.22237753,1.993,4.43832809,0.959942941,4.62353323,2.209
Hexamerin,Hex,LW4,LW,Late,6,1,1,Target,113,25.90912056,25.77098846,25.75769424,25.81260109,25.75769424,25.90912056,0.151426315,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,25.8126,25.96422745,27.25065,26.22237753,1.993,2.45991587,1.194862228,2.058744355,1.04176
Hexamerin,Hex,LW5,LW,Late,6,1,1,Target,113,24.87003326,24.87427139,24.87405586,24.87278684,24.87003326,24.87427139,0.004238129,eq7,24.70294189,24.60329819,24.67490959,24.66038322,24.60329819,24.70294189,0.099643707,24.66038322,1.87,24.87279,25.98699537,27.25065,26.22237753,1.993,4.42997086,1.176247401,3.766189709,1.91311
Hexamerin,Hex,LW6,LW,Late,6,2,1,Target,27,24.96445656,24.91424751,24.89943695,24.92604701,24.89943695,24.96445656,0.065019608,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,25.03292,25.6146532,27.25065,26.22237753,1.993,4.00747094,1.520610107,2.635436212,1.39804
Hexamerin,Hex,LW7,LW,Late,6,2,1,Target,27,24.5117588,24.44429588,24.40208435,24.45271301,24.40208435,24.5117588,0.109674454,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,24.55756,25.49879946,27.25065,26.22237753,1.993,5.39627337,1.647088425,3.276249947,1.71205
Hexamerin,Hex,LW8,LW,Late,6,2,1,Target,27,25.60618591,25.32860947,25.52929878,25.48803139,25.32860947,25.60618591,0.277576447,eq7,24.60397339,24.71645164,24.34488106,24.55510203,24.34488106,24.71645164,0.371570587,24.66038322,1.87,25.59731,25.38152329,27.25065,26.22237753,1.993,2.81478114,1.78583768,1.576168527,0.65642
kruppel,Kruppel,EQ3,EQ,Early,1,1,1,Target,45,32.91544724,32.49277496,32.93657303,32.78159841,32.49277496,32.93657303,0.443798065,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,33.26634,25.20347779,34.66094,26.22237753,1.993,2.48546978,2.019146866,1.230950467,0.29977
kruppel,Kruppel,EQ4,EQ,Early,1,1,1,Target,45,31.89891052,31.67499352,31.29744148,31.62378184,31.29744148,31.89891052,0.60146904,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,32.0914,25.71323328,34.66094,26.22237753,1.993,5.35223687,1.420667419,3.767410163,1.91357
kruppel,Kruppel,EQ5,EQ,Early,1,1,1,Target,45,31.35998535,31.41131401,31.1414566,31.30425199,31.1414566,31.41131401,0.269857407,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,31.76715,25.41256952,34.66094,26.22237753,1.993,6.61409479,1.748007975,3.783789824,1.91983
kruppel,Kruppel,EQ6,EQ,Early,1,2,1,Target,110,31.3805027,31.17586327,31.02655029,31.19430542,31.02655029,31.3805027,0.353952408,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,31.19431,24.43307804,34.66094,26.22237753,1.993,9.61361402,3.434853817,2.798842261,1.48483
kruppel,Kruppel,EQ7,EQ,Early,1,2,1,Target,110,31.41065979,31.72662926,32.04882431,31.72870445,31.41065979,32.04882431,0.63816452,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,31.7287,25.14049381,34.66094,26.22237753,1.993,6.78218969,2.108783874,3.216161586,1.68534
kruppel,Kruppel,EQ8,EQ,Early,1,2,1,Target,110,31.45041656,31.18271446,31.01560593,31.21624565,31.01560593,31.45041656,0.434810638,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,31.21625,24.92400707,34.66094,26.22237753,1.993,9.47689406,2.448338325,3.870745297,1.95261
kruppel,Kruppel,EW3,EW,Early,2,1,1,Target,45,30.23582649,29.98538971,29.86480522,30.02867381,29.86480522,30.23582649,0.371021271,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,30.47271,24.92159537,34.66094,26.22237753,1.993,15.39854228,2.452413802,6.278933135,2.65052
kruppel,Kruppel,EW4,EW,Early,2,1,1,Target,45,29.98113251,30.22814178,29.85948944,30.02292124,29.85948944,30.22814178,0.368652344,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,30.46687,24.71233996,34.66094,26.22237753,1.993,15.45733918,2.833135032,5.455913327,2.44782
kruppel,Kruppel,EW5,EW,Early,2,1,1,Target,45,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,30.94604,25.76160569,34.66094,26.22237753,1.993,11.30515206,1.374056303,8.227575557,3.04047
kruppel,Kruppel,EW6,EW,Early,2,2,1,Target,110,30.99755096,30.54016876,30.39030075,30.64267349,30.39030075,30.99755096,0.607250214,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,30.64267,24.47244454,34.66094,26.22237753,1.993,13.78127593,3.342856272,4.122604986,2.04356
kruppel,Kruppel,EW7,EW,Early,2,2,1,Target,110,30.94189453,30.79017448,30.98571968,30.90592957,30.79017448,30.98571968,0.195545197,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,30.90593,24.71929254,34.66094,26.22237753,1.993,11.60510643,2.819583273,4.115894196,2.04121
kruppel,Kruppel,EW8,EW,Early,2,2,1,Target,110,30.09348869,30.0083046,30.23609352,30.11262894,30.0083046,30.23609352,0.227788925,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,30.11263,24.24981237,34.66094,26.22237753,1.993,19.47922009,3.897603549,4.997742803,2.32128
kruppel,Kruppel,MQ3,MQ,Medium,3,1,1,Target,45,30.83205605,31.08441162,31.06331253,30.99326007,30.83205605,31.08441162,0.252355576,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,31.45156,25.89835631,34.66094,26.22237753,1.993,8.12735017,1.25039354,6.499833779,2.7004
kruppel,Kruppel,MQ4,MQ,Medium,3,1,1,Target,45,30.92617798,30.7182827,30.8628788,30.83577983,30.7182827,30.92617798,0.207895279,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,31.29175,25.52299076,34.66094,26.22237753,1.993,9.02109314,1.619837474,5.569134735,2.47745
kruppel,Kruppel,MQ5,MQ,Medium,3,1,1,Target,45,30.4221859,30.80226517,30.6393261,30.62125905,30.4221859,30.80226517,0.380079269,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,31.07405,25.12014958,34.66094,26.22237753,1.993,10.39874807,2.138579087,4.862456636,2.28169
kruppel,Kruppel,MQ6,MQ,Medium,3,2,1,Target,110,31.48898506,31.21654129,31.44927597,31.38493411,31.21654129,31.48898506,0.272443771,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,31.38493,24.7511079,34.66094,26.22237753,1.993,8.48864228,2.75839205,3.077387885,1.62171
kruppel,Kruppel,MQ7,MQ,Medium,3,2,1,Target,110,30.89920425,30.29115868,30.21905327,30.4698054,30.21905327,30.89920425,0.680150986,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,30.46981,24.76237619,34.66094,26.22237753,1.993,15.4277383,2.737039456,5.636651771,2.49484
kruppel,Kruppel,MQ8,MQ,Medium,3,2,1,Target,110,29.67650986,29.74713516,29.52760124,29.65041542,29.52760124,29.74713516,0.21953392,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,29.65042,24.74641872,34.66094,26.22237753,1.993,26.34036234,2.767326721,9.518342064,3.25071
kruppel,Kruppel,MW3,MW,Medium,4,1,1,Target,45,32.49047852,32.20552063,31.66767693,32.12122536,31.66767693,32.49047852,0.82280159,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,32.5962,25.29041589,34.66094,26.22237753,1.993,3.84955863,1.901644416,2.024331468,1.01745
kruppel,Kruppel,MW4,MW,Medium,4,1,1,Target,45,30.52569962,30.23775673,30.51194954,30.42513529,30.23775673,30.52569962,0.287942886,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,30.87503,24.63093495,34.66094,26.22237753,1.993,11.84158291,2.996737273,3.951491851,1.9824
kruppel,Kruppel,MW5,MW,Medium,4,1,1,Target,45,31.53616905,31.22352219,NA,31.37984562,31.22352219,31.53616905,0.312646866,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,31.84386,25.35247425,34.66094,26.22237753,1.993,6.29101361,1.821974874,3.452854206,1.78779
kruppel,Kruppel,MW6,MW,Medium,4,2,1,Target,110,33.09719467,33.32415771,33.21341705,33.21158981,33.09719467,33.32415771,0.226963043,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,33.21159,24.95817119,34.66094,26.22237753,1.993,2.57591396,2.391327326,1.077190034,0.10727
kruppel,Kruppel,MW7,MW,Medium,4,2,1,Target,110,34.57814026,34.18028641,34.81772232,34.525383,34.18028641,34.81772232,0.637435913,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,34.52538,25.70248086,34.66094,26.22237753,1.993,1.09253418,1.431241263,0.763347321,-0.38959
kruppel,Kruppel,MW8,MW,Medium,4,2,1,Target,110,30.7471981,30.74035645,30.69045448,30.72600301,30.69045448,30.7471981,0.056743622,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,30.726,24.59640071,34.66094,26.22237753,1.993,13.0515842,3.068964925,4.252764213,2.0884
kruppel,Kruppel,LQ3,LQ,Late,5,1,1,Target,45,34.00530243,34.45130539,34.01103592,34.15588125,34.00530243,34.45130539,0.44600296,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,34.66094,26.22237753,34.66094,26.22237753,1.993,1,1,1,0
kruppel,Kruppel,LQ4,LQ,Late,5,1,1,Target,45,34.52342987,34.36223221,34.27435684,34.38667297,34.27435684,34.52342987,0.249073029,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,34.89515,25.57874041,34.66094,26.22237753,1.993,0.85821605,1.558741203,0.550582771,-0.86097
kruppel,Kruppel,LQ5,LQ,Late,5,1,1,Target,45,36.28503799,NA,36.51421738,36.39962769,36.28503799,36.51421738,0.229179382,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,36.93787,26.77102083,34.66094,26.22237753,1.993,0.22616716,0.684979114,0.330181111,-1.59867
kruppel,Kruppel,LQ6,LQ,Late,5,2,1,Target,110,35.30870438,35.89676285,NA,35.60273361,35.30870438,35.89676285,0.588058472,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,35.60273,26.04640784,34.66094,26.22237753,1.993,0.54072525,1.129026689,0.478930444,-1.06211
kruppel,Kruppel,LQ7,LQ,Late,5,2,1,Target,110,32.25599289,32.5526123,32.76990128,32.52616882,32.25599289,32.76990128,0.513908386,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,32.52617,24.94259361,34.66094,26.22237753,1.993,4.0296466,2.417155701,1.667102621,0.73734
kruppel,Kruppel,LQ8,LQ,Late,5,2,1,Target,110,33.46296692,33.62502289,33.5421257,33.54337184,33.46296692,33.62502289,0.162055969,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,33.54337,25.43112605,34.66094,26.22237753,1.993,2.07425464,1.725780648,1.20192253,0.26534
kruppel,Kruppel,LW3,LW,Late,6,1,1,Target,45,33.17163849,33.20108032,32.93103027,33.10124969,32.93103027,33.20108032,0.270050049,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,33.59072,26.28165682,34.66094,26.22237753,1.993,2.01112219,0.959942941,2.095043472,1.06698
kruppel,Kruppel,LW4,LW,Late,6,1,1,Target,45,33.1922493,33.56058121,33.03300858,33.26194636,33.03300858,33.56058121,0.527572632,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,33.75379,25.96422745,34.66094,26.22237753,1.993,1.80801835,1.194862228,1.51316052,0.59757
kruppel,Kruppel,LW5,LW,Late,6,1,1,Target,45,32.20344162,32.52826309,32.51473999,32.41548157,32.20344162,32.52826309,0.324821472,ew5,30.46265984,30.56473923,30.45793343,30.49511083,30.45793343,30.56473923,0.106805801,30.94604111,1.921,32.89481,25.98699537,34.66094,26.22237753,1.993,3.16772578,1.176247401,2.693077811,1.42926
kruppel,Kruppel,LW6,LW,Late,6,2,1,Target,110,35.69690323,35.31335831,34.87540817,35.29522324,34.87540817,35.69690323,0.821495056,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,35.29522,25.6146532,34.66094,26.22237753,1.993,0.66094338,1.520610107,0.434656708,-1.20205
kruppel,Kruppel,LW7,LW,Late,6,2,1,Target,110,NA,33.33543396,33.34042358,33.33792877,33.33543396,33.34042358,0.004989624,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,33.33793,25.49879946,34.66094,26.22237753,1.993,2.37197879,1.647088425,1.440104098,0.52617
kruppel,Kruppel,LW8,LW,Late,6,2,1,Target,110,31.25287819,31.45519829,31.29309464,31.3337237,31.25287819,31.45519829,0.202320099,ew5,30.79952621,30.90634537,31.13225174,30.94604111,30.79952621,31.13225174,0.332725525,30.94604111,1.921,31.33372,25.38152329,34.66094,26.22237753,1.993,8.77723615,1.78583768,4.914912618,2.29717
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EQ3,EQ,Early,1,1,1,Target,56,33.51232529,33.15142822,33.24727249,33.30367533,33.15142822,33.51232529,0.360897064,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,33.52877,25.20347779,30.76667,26.22237753,1.993,0.16172614,2.019146866,0.080096272,-3.64212
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EQ4,EQ,Early,1,1,1,Target,56,33.80469894,33.15742111,NA,33.48106003,33.15742111,33.80469894,0.647277832,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,33.70735,25.71323328,30.76667,26.22237753,1.993,0.14375527,1.420667419,0.101188547,-3.30488
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EQ5,EQ,Early,1,1,1,Target,56,34.0315361,34.717556,34.11260605,34.28723272,34.0315361,34.717556,0.686019897,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,34.51897,25.41256952,30.76667,26.22237753,1.993,0.08416462,1.748007975,0.048148877,-4.37635
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EQ6,EQ,Early,1,2,1,Target,126,33.75477982,34.18318558,33.32989883,33.75595474,33.32989883,34.18318558,0.853286743,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,33.75595,24.43307804,30.76667,26.22237753,1.993,0.13921981,3.434853817,0.040531508,-4.62481
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EQ7,EQ,Early,1,2,1,Target,126,33.70908356,34.02108002,33.81645966,33.84887441,33.70908356,34.02108002,0.31199646,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,33.84887,25.14049381,30.76667,26.22237753,1.993,0.13094339,2.108783874,0.062094268,-4.0094
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EQ8,EQ,Early,1,2,1,Target,126,32.5145607,32.87401199,32.69535446,32.69464238,32.5145607,32.87401199,0.359451294,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,32.69464,24.92400707,30.76667,26.22237753,1.993,0.28036307,2.448338325,0.114511573,-3.12643
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EW3,EW,Early,2,1,1,Target,56,36.19804382,36.53982162,NA,36.36893272,36.19804382,36.53982162,0.341777802,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,36.61474,24.92159537,30.76667,26.22237753,1.993,0.02112432,2.452413802,0.008613685,-6.85915
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EW4,EW,Early,2,1,1,Target,56,34.78593826,NA,34.41075134,34.5983448,34.41075134,34.78593826,0.37518692,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,34.83219,24.71233996,30.76667,26.22237753,1.993,0.06845531,2.833135032,0.02416239,-5.37109
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EW5,EW,Early,2,1,1,Target,56,34.11764145,34.19484711,34.0586586,34.12371572,34.0586586,34.19484711,0.136188507,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,34.35435,25.76160569,30.76667,26.22237753,1.993,0.09381811,1.374056303,0.06827821,-3.87243
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EW6,EW,Early,2,2,1,Target,126,33.57571411,33.76225281,33.84280777,33.7269249,33.57571411,33.84280777,0.267093658,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,33.72692,24.47244454,30.76667,26.22237753,1.993,0.14191124,3.342856272,0.042452093,-4.55802
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EW7,EW,Early,2,2,1,Target,126,33.31281662,33.41157532,34.28292847,33.6691068,33.31281662,34.28292847,0.970111847,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,33.66911,24.71929254,30.76667,26.22237753,1.993,0.14742773,2.819583273,0.052287063,-4.2574
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,EW8,EW,Early,2,2,1,Target,126,32.438694,32.03370667,32.3399086,32.27076976,32.03370667,32.438694,0.404987335,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,32.27077,24.24981237,30.76667,26.22237753,1.993,0.37080183,3.897603549,0.095135851,-3.39387
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MQ3,MQ,Medium,3,1,1,Target,56,33.09233093,33.88218307,33.95846176,33.64432526,33.09233093,33.95846176,0.866130829,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,33.87172,25.89835631,30.76667,26.22237753,1.993,0.12898503,1.25039354,0.103155546,-3.27711
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MQ4,MQ,Medium,3,1,1,Target,56,33.32409286,33.17493057,33.36858749,33.28920364,33.17493057,33.36858749,0.193656921,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,33.5142,25.52299076,30.76667,26.22237753,1.993,0.1632878,1.619837474,0.100805051,-3.31036
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MQ5,MQ,Medium,3,1,1,Target,56,33.55887222,33.79961014,33.36048126,33.57298787,33.36048126,33.79961014,0.439128876,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,33.7999,25.12014958,30.76667,26.22237753,1.993,0.13524228,2.138579087,0.063239315,-3.98303
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MQ6,MQ,Medium,3,2,1,Target,126,33.52081299,32.63476944,33.02601242,33.06053162,32.63476944,33.52081299,0.886043549,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,33.06053,24.7511079,30.76667,26.22237753,1.993,0.22024673,2.75839205,0.079846058,-3.64663
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MQ7,MQ,Medium,3,2,1,Target,126,32.73051834,33.61050797,33.53246689,33.2911644,32.73051834,33.61050797,0.879989624,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,33.29116,24.76237619,30.76667,26.22237753,1.993,0.18916606,2.737039456,0.069113385,-3.85489
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MQ8,MQ,Medium,3,2,1,Target,126,35.81987762,35.49755478,NA,35.6587162,35.49755478,35.81987762,0.322322845,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,35.65872,24.74641872,30.76667,26.22237753,1.993,0.03968649,2.767326721,0.014341093,-6.1237
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MW3,MW,Medium,4,1,1,Target,56,33.85360336,33.34550476,33.43572998,33.54494603,33.34550476,33.85360336,0.508098602,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,33.77167,25.29041589,30.76667,26.22237753,1.993,0.13778424,1.901644416,0.07245531,-3.78676
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MW4,MW,Medium,4,1,1,Target,56,36.69185638,36.11672211,NA,36.40428925,36.11672211,36.69185638,0.575134277,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,36.65034,24.63093495,30.76667,26.22237753,1.993,0.02063413,2.996737273,0.006885533,-7.18222
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MW5,MW,Medium,4,1,1,Target,56,34.2386055,34.40821075,NA,34.32340813,34.2386055,34.40821075,0.169605255,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,34.55539,25.35247425,30.76667,26.22237753,1.993,0.08216689,1.821974874,0.045097708,-4.4708
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MW6,MW,Medium,4,2,1,Target,126,NA,38.0697937,37.29785919,37.68382645,37.29785919,38.0697937,0.771934509,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,37.68383,24.95817119,30.76667,26.22237753,1.993,0.01043606,2.391327326,0.004364131,-7.84009
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MW7,MW,Medium,4,2,1,Target,126,38.9122963,38.15637589,NA,38.53433609,38.15637589,38.9122963,0.75592041,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,38.53434,25.70248086,30.76667,26.22237753,1.993,0.00595529,1.431241263,0.004160925,-7.90888
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,MW8,MW,Medium,4,2,1,Target,126,32.74189377,32.95613098,32.42763138,32.70855204,32.42763138,32.95613098,0.528499603,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,32.70855,24.59640071,30.76667,26.22237753,1.993,0.2778026,3.068964925,0.090519965,-3.46562
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LQ3,LQ,Late,5,1,1,Target,56,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,30.76667,26.22237753,30.76667,26.22237753,1.993,1,1,1,0
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LQ4,LQ,Late,5,1,1,Target,56,27.78791618,27.5602684,27.54999924,27.63272794,27.54999924,27.78791618,0.237916946,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,27.81949,25.57874041,30.76667,26.22237753,1.993,6.98617224,1.558741203,4.481932105,2.16412
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LQ5,LQ,Late,5,1,1,Target,56,30.59613991,30.62167931,30.31499672,30.51093864,30.31499672,30.62167931,0.306682587,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,30.71716,26.77102083,30.76667,26.22237753,1.993,1.03320075,0.684979114,1.508368259,0.59299
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LQ6,LQ,Late,5,2,1,Target,126,28.6329422,28.70641327,28.6814003,28.67358526,28.6329422,28.70641327,0.073471069,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,28.67359,26.04640784,30.76667,26.22237753,1.993,3.9772111,1.129026689,3.522690062,1.81668
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LQ7,LQ,Late,5,2,1,Target,126,27.43013954,27.37080765,27.2889595,27.36330223,27.2889595,27.43013954,0.141180038,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,27.3633,24.94259361,30.76667,26.22237753,1.993,9.43881974,2.417155701,3.904928316,1.9653
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LQ8,LQ,Late,5,2,1,Target,126,27.67730331,27.6722908,27.72517204,27.69158872,27.6722908,27.72517204,0.052881241,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,27.69159,25.43112605,30.76667,26.22237753,1.993,7.60112551,1.725780648,4.404456334,2.13896
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LW3,LW,Late,6,1,1,Target,56,35.97676468,NA,36.25333023,36.11504745,35.97676468,36.25333023,0.276565552,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,36.35914,26.28165682,30.76667,26.22237753,1.993,0.02500353,0.959942941,0.026046895,-5.26274
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LW4,LW,Late,6,1,1,Target,56,35.46699905,35.80855942,34.81481552,35.363458,34.81481552,35.80855942,0.993743896,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,35.60247,25.96422745,30.76667,26.22237753,1.993,0.04118643,1.194862228,0.03446961,-4.85853
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LW5,LW,Late,6,1,1,Target,56,34.9569397,NA,35.18067932,35.06880951,34.9569397,35.18067932,0.223739624,lq3,30.29889679,30.30713654,31.07433891,30.56012408,30.29889679,31.07433891,0.775442123,30.76667341,1.934,35.30583,25.98699537,30.76667,26.22237753,1.993,0.05008742,1.176247401,0.042582382,-4.5536
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LW6,LW,Late,6,2,1,Target,126,35.91853333,NA,36.51411819,36.21632576,35.91853333,36.51411819,0.595584869,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,36.21633,25.6146532,30.76667,26.22237753,1.993,0.02747336,1.520610107,0.018067329,-5.79047
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LW7,LW,Late,6,2,1,Target,126,36.31609344,36.15526199,NA,36.23567772,36.15526199,36.31609344,0.160831451,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,36.23568,25.49879946,30.76667,26.22237753,1.993,0.02712491,1.647088425,0.0164684,-5.92416
nose resistant to fluoxetine protein 6-like (LOC100640031),Nos.res.640031,LW8,LW,Late,6,2,1,Target,126,33.92953491,34.44836807,NA,34.18895149,33.92953491,34.44836807,0.51883316,lq3,30.86022186,30.65842247,30.78137589,30.76667341,30.65842247,30.86022186,0.201799393,30.76667341,1.934,34.18895,25.38152329,30.76667,26.22237753,1.993,0.10463243,1.78583768,0.058590111,-4.0932
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EQ3,EQ,Early,1,1,1,Target,55,31.72030449,31.58748817,31.91995811,31.74258359,31.58748817,31.91995811,0.33246994,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,31.74258,25.20347779,29.96419,26.22237753,1.993,0.3282053,2.019146866,0.162546522,-2.62108
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EQ4,EQ,Early,1,1,1,Target,55,30.19960785,30.07454109,29.9734211,30.08252335,29.9734211,30.19960785,0.226186752,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,30.08252,25.71323328,29.96419,26.22237753,1.993,0.92854824,1.420667419,0.653600009,-0.61352
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EQ5,EQ,Early,1,1,1,Target,55,34.01721954,34.54887772,33.7289505,34.09834925,33.7289505,34.54887772,0.819927216,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,34.09835,25.41256952,29.96419,26.22237753,1.993,0.07502457,1.748007975,0.042920038,-4.5422
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EQ6,EQ,Early,1,2,1,Target,115,30.67302513,29.9975872,30.16711807,30.27924347,29.9975872,30.67302513,0.675437927,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,30.36216,24.43307804,29.96419,26.22237753,1.993,0.77933245,3.434853817,0.226889554,-2.13994
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EQ7,EQ,Early,1,2,1,Target,115,30.9647522,31.00620651,30.94957733,30.97351201,30.94957733,31.00620651,0.056629181,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,31.05833,25.14049381,29.96419,26.22237753,1.993,0.50386388,2.108783874,0.238935763,-2.06531
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EQ8,EQ,Early,1,2,1,Target,115,30.65930939,30.55534935,30.46989441,30.56151772,30.46989441,30.65930939,0.189414978,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,30.64521,24.92400707,29.96419,26.22237753,1.993,0.65269913,2.448338325,0.266588616,-1.90731
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EW3,EW,Early,2,1,1,Target,55,36.56204224,NA,36.5328331,36.54743767,36.5328331,36.56204224,0.029209137,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,36.54744,24.92159537,29.96419,26.22237753,1.993,0.01617599,2.452413802,0.006595945,-7.24421
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EW4,EW,Early,2,1,1,Target,55,30.92895889,30.97549248,30.35422707,30.75289281,30.35422707,30.97549248,0.621265411,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,30.75289,24.71233996,29.96419,26.22237753,1.993,0.61011895,2.833135032,0.21535117,-2.21524
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EW5,EW,Early,2,1,1,Target,55,30.50347519,30.48266983,30.40060616,30.46225039,30.40060616,30.50347519,0.102869034,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,30.46225,25.76160569,29.96419,26.22237753,1.993,0.73196564,1.374056303,0.532704255,-0.90859
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EW6,EW,Early,2,2,1,Target,115,30.41295624,29.95759773,30.16943932,30.17999776,29.95759773,30.41295624,0.455358505,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,30.26264,24.47244454,29.96419,26.22237753,1.993,0.82946652,3.342856272,0.248131074,-2.01083
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EW7,EW,Early,2,2,1,Target,115,30.83395386,30.20711708,30.43040848,30.49049314,30.20711708,30.83395386,0.626836777,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,30.57399,24.71929254,29.96419,26.22237753,1.993,0.68247991,2.819583273,0.242049923,-2.04662
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,EW8,EW,Early,2,2,1,Target,115,30.21242523,29.95158577,30.00059319,30.05486806,29.95158577,30.21242523,0.260839462,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,30.13717,24.24981237,29.96419,26.22237753,1.993,0.89729779,3.897603549,0.230217819,-2.11893
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MQ3,MQ,Medium,3,1,1,Target,55,28.46281433,28.40192986,28.23742294,28.36738904,28.23742294,28.46281433,0.225391388,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,28.36739,25.89835631,29.96419,26.22237753,1.993,2.7192401,1.25039354,2.174707413,1.12082
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MQ4,MQ,Medium,3,1,1,Target,55,28.11639977,28.03235435,28.10876846,28.08584086,28.03235435,28.11639977,0.08404541,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,28.08584,25.52299076,29.96419,26.22237753,1.993,3.2437654,1.619837474,2.002525221,1.00182
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MQ5,MQ,Medium,3,1,1,Target,55,30.34142494,30.11122322,29.79947281,30.08404032,29.79947281,30.34142494,0.541952133,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,30.08404,25.12014958,29.96419,26.22237753,1.993,0.92766622,2.138579087,0.433776904,-1.20497
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MQ6,MQ,Medium,3,2,1,Target,115,28.11935997,28.03007126,28.87298393,28.34080505,28.03007126,28.87298393,0.842912674,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,28.41841,24.7511079,29.96419,26.22237753,1.993,2.63369368,2.75839205,0.954793094,-0.06674
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MQ7,MQ,Medium,3,2,1,Target,115,29.85857582,29.76485252,29.93108177,29.85150337,29.76485252,29.93108177,0.166229248,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,29.93325,24.76237619,29.96419,26.22237753,1.993,1.01957297,2.737039456,0.372509417,-1.42465
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MQ8,MQ,Medium,3,2,1,Target,115,35.25043869,34.46572113,34.3731308,34.69643021,34.3731308,35.25043869,0.877307892,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,34.79144,24.74641872,29.96419,26.22237753,1.993,0.04859942,2.767326721,0.017561864,-5.83141
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MW3,MW,Medium,4,1,1,Target,55,36.16529846,36.21193695,36.85498047,36.41073863,36.16529846,36.85498047,0.689682007,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,36.41074,25.29041589,29.96419,26.22237753,1.993,0.01762232,1.901644416,0.009266883,-6.7537
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MW4,MW,Medium,4,1,1,Target,55,36.08677292,35.51004028,36.23939896,35.94540405,35.51004028,36.23939896,0.729358673,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,35.9454,24.63093495,29.96419,26.22237753,1.993,0.02358676,2.996737273,0.007870814,-6.98927
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MW5,MW,Medium,4,1,1,Target,55,35.15299988,34.74593735,34.51652908,34.80515544,34.51652908,35.15299988,0.636470795,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,34.80516,25.35247425,29.96419,26.22237753,1.993,0.04818369,1.821974874,0.02644586,-5.24081
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MW6,MW,Medium,4,2,1,Target,115,39.02393341,40,NA,39.51196671,39.02393341,40,0.976066589,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,39.62017,24.95817119,29.96419,26.22237753,1.993,0.00235973,2.391327326,0.000986786,-9.98498
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MW7,MW,Medium,4,2,1,Target,115,40,40,NA,40,40,40,0,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,40.10954,25.70248086,29.96419,26.22237753,1.993,0.00173667,1.431241263,0.001213401,-9.68673
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,MW8,MW,Medium,4,2,1,Target,115,32.83042526,33.13126373,32.5880661,32.84991837,32.5880661,33.13126373,0.543197632,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,32.93987,24.59640071,29.96419,26.22237753,1.993,0.15502232,3.068964925,0.050512902,-4.3072
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LQ3,LQ,Late,5,1,1,Target,55,30.0441494,29.95932579,29.8890934,29.96418953,29.8890934,30.0441494,0.155056,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,29.96419,26.22237753,29.96419,26.22237753,1.993,1,1,1,0
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LQ4,LQ,Late,5,1,1,Target,55,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,27.75657,25.57874041,29.96419,26.22237753,1.993,3.98688902,1.558741203,2.557762003,1.35488
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LQ5,LQ,Late,5,1,1,Target,55,31.6458683,31.64800835,31.45392799,31.58260155,31.45392799,31.64800835,0.194080353,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,31.5826,26.77102083,29.96419,26.22237753,1.993,0.36280437,0.684979114,0.52965757,-0.91687
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LQ6,LQ,Late,5,2,1,Target,115,28.68220901,28.24271393,28.69131851,28.53874715,28.24271393,28.69131851,0.448604584,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,28.6169,26.04640784,29.96419,26.22237753,1.993,2.32574947,1.129026689,2.059959689,1.04262
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LQ7,LQ,Late,5,2,1,Target,115,27.99882317,27.89862251,27.96938324,27.95560964,27.89862251,27.99882317,0.100200653,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,28.03216,24.94259361,29.96419,26.22237753,1.993,3.35470104,2.417155701,1.387871308,0.47287
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LQ8,LQ,Late,5,2,1,Target,115,28.87068748,28.93596649,28.73880577,28.84848658,28.73880577,28.93596649,0.197160721,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,28.92748,25.43112605,29.96419,26.22237753,1.993,1.9145212,1.725780648,1.109365322,0.14973
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LW3,LW,Late,6,1,1,Target,55,40,40,40,40,40,40,0,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,40,26.28165682,29.96419,26.22237753,1.993,0.00186003,0.959942941,0.001937642,-9.01148
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LW4,LW,Late,6,1,1,Target,55,38.09527969,NA,37.36203003,37.72865486,37.36203003,38.09527969,0.733249664,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,37.72865,25.96422745,29.96419,26.22237753,1.993,0.00771778,1.194862228,0.00645914,-7.27444
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LW5,LW,Late,6,1,1,Target,55,35.12368774,NA,35.41534424,35.26951599,35.12368774,35.41534424,0.291656494,lq4,27.92819405,27.62391663,27.71761322,27.75657463,27.62391663,27.92819405,0.30427742,27.75657463,1.871,35.26952,25.98699537,29.96419,26.22237753,1.993,0.03602133,1.176247401,0.03062394,-5.0292
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LW6,LW,Late,6,2,1,Target,115,40,40,40,40,40,40,0,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,40.10954,25.6146532,29.96419,26.22237753,1.993,0.00173667,1.520610107,0.001142088,-9.77411
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LW7,LW,Late,6,2,1,Target,115,38.62946701,40,NA,39.31473351,38.62946701,40,1.37053299,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,39.42239,25.49879946,29.96419,26.22237753,1.993,0.00267098,1.647088425,0.001621638,-9.26833
nose resistant to fluoxetine protein 6-like (LOC100645614),Nos.res.645614,LW8,LW,Late,6,2,1,Target,115,NA,36.22451019,35.21463013,35.71957016,35.21463013,36.22451019,1.009880066,lq4,27.56847763,27.52442741,27.94941711,27.68077405,27.52442741,27.94941711,0.4249897,27.75657463,1.871,35.81738,25.38152329,29.96419,26.22237753,1.993,0.02555637,1.78583768,0.014310579,-6.12677
Novel1,Novel1,EQ3,EQ,Early,1,1,1,Target,128,29.33765602,29.56222153,29.72774315,29.54254023,29.33765602,29.72774315,0.390087128,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.54254,25.20347779,28.70086,26.22237753,1.993,0.56345106,2.019146866,0.279054023,-1.84138
Novel1,Novel1,EQ4,EQ,Early,1,1,1,Target,128,29.24455643,29.74517059,29.86601639,29.61858114,29.24455643,29.86601639,0.621459961,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.61858,25.71323328,28.70086,26.22237753,1.993,0.53499236,1.420667419,0.376578189,-1.40898
Novel1,Novel1,EQ5,EQ,Early,1,1,1,Target,128,29.95094681,29.94473648,29.5643177,29.82000033,29.5643177,29.95094681,0.386629105,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.82,25.41256952,28.70086,26.22237753,1.993,0.46636522,1.748007975,0.26679811,-1.90618
Novel1,Novel1,EQ6,EQ,Early,1,2,1,Target,101,29.53241539,29.66643715,30.19103432,29.79662895,29.53241539,30.19103432,0.658618927,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,30.03821,24.43307804,28.70086,26.22237753,1.993,0.40191504,3.434853817,0.117010813,-3.09529
Novel1,Novel1,EQ7,EQ,Early,1,2,1,Target,101,29.713377,29.82469559,29.60666084,29.71491114,29.60666084,29.82469559,0.218034744,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.95583,25.14049381,28.70086,26.22237753,1.993,0.42512769,2.108783874,0.201598509,-2.31044
Novel1,Novel1,EQ8,EQ,Early,1,2,1,Target,101,28.46555328,28.69827652,28.50713921,28.55698967,28.46555328,28.69827652,0.232723236,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.78852,24.92400707,28.70086,26.22237753,1.993,0.94199978,2.448338325,0.384750659,-1.378
Novel1,Novel1,EW3,EW,Early,2,1,1,Target,128,30.04463959,29.69267082,30.53115654,30.08948898,29.69267082,30.53115654,0.838485718,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,30.08949,24.92159537,28.70086,26.22237753,1.993,0.38811079,2.452413802,0.15825665,-2.65966
Novel1,Novel1,EW4,EW,Early,2,1,1,Target,128,28.86208534,28.85228729,28.83214569,28.84883944,28.83214569,28.86208534,0.029939651,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,28.84884,24.71233996,28.70086,26.22237753,1.993,0.90405793,2.833135032,0.319101602,-1.64791
Novel1,Novel1,EW5,EW,Early,2,1,1,Target,128,29.75973129,30.01478577,29.93586159,29.90345955,29.75973129,30.01478577,0.255054474,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.90346,25.76160569,28.70086,26.22237753,1.993,0.44057685,1.374056303,0.32063959,-1.64098
Novel1,Novel1,EW6,EW,Early,2,2,1,Target,101,29.26968384,29.05784225,29.0649929,29.13083967,29.05784225,29.26968384,0.211841583,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.36702,24.47244454,28.70086,26.22237753,1.993,0.6350533,3.342856272,0.189973258,-2.39613
Novel1,Novel1,EW7,EW,Early,2,2,1,Target,101,29.03287697,29.08844566,29.26088905,29.1274039,29.03287697,29.26088905,0.228012085,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.36356,24.71929254,28.70086,26.22237753,1.993,0.63655426,2.819583273,0.225761824,-2.14713
Novel1,Novel1,EW8,EW,Early,2,2,1,Target,101,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.58331,24.24981237,28.70086,26.22237753,1.993,1.08341677,3.897603549,0.277969977,-1.847
Novel1,Novel1,MQ3,MQ,Medium,3,1,1,Target,128,29.62715721,29.99591064,30.09660912,29.90655899,29.62715721,30.09660912,0.469451904,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.90656,25.89835631,28.70086,26.22237753,1.993,0.4396471,1.25039354,0.351606987,-1.50796
Novel1,Novel1,MQ4,MQ,Medium,3,1,1,Target,128,29.55331039,29.3257103,29.87277603,29.58393224,29.3257103,29.87277603,0.547065735,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.58393,25.52299076,28.70086,26.22237753,1.993,0.54777712,1.619837474,0.33816795,-1.56419
Novel1,Novel1,MQ5,MQ,Medium,3,1,1,Target,128,29.33036995,29.58771324,29.39211655,29.43673325,29.33036995,29.58771324,0.257343292,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.43673,25.12014958,28.70086,26.22237753,1.993,0.60558592,2.138579087,0.283172096,-1.82025
Novel1,Novel1,MQ6,MQ,Medium,3,2,1,Target,101,28.4405117,28.58866119,28.85115623,28.62677638,28.4405117,28.85115623,0.410644531,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.85887,24.7511079,28.70086,26.22237753,1.993,0.8978959,2.75839205,0.325514243,-1.61921
Novel1,Novel1,MQ7,MQ,Medium,3,2,1,Target,101,28.02380562,28.23121262,27.84005737,28.03169187,27.84005737,28.23121262,0.391155243,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.25896,24.76237619,28.70086,26.22237753,1.993,1.35145982,2.737039456,0.49376702,-1.0181
Novel1,Novel1,MQ8,MQ,Medium,3,2,1,Target,101,28.64712715,28.96368217,28.54988289,28.72023074,28.54988289,28.96368217,0.413799286,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.95309,24.74641872,28.70086,26.22237753,1.993,0.84205138,2.767326721,0.304283325,-1.71651
Novel1,Novel1,MW3,MW,Medium,4,1,1,Target,128,29.49322319,28.9376297,29.10123253,29.17736181,28.9376297,29.49322319,0.555593491,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.17736,25.29041589,28.70086,26.22237753,1.993,0.72268901,1.901644416,0.380033722,-1.3958
Novel1,Novel1,MW4,MW,Medium,4,1,1,Target,128,29.17041969,29.23989487,28.68761063,29.03264173,28.68761063,29.23989487,0.552284241,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.03264,24.63093495,28.70086,26.22237753,1.993,0.7976081,2.996737273,0.266158834,-1.90964
Novel1,Novel1,MW5,MW,Medium,4,1,1,Target,128,30.06271553,30.00787926,29.73320389,29.93459956,29.73320389,30.06271553,0.329511642,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.9346,25.35247425,28.70086,26.22237753,1.993,0.4313244,1.821974874,0.236734547,-2.07866
Novel1,Novel1,MW6,MW,Medium,4,2,1,Target,101,29.27854538,29.38697624,29.1701355,29.27855237,29.1701355,29.38697624,0.216840744,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.51593,24.95817119,28.70086,26.22237753,1.993,0.57376193,2.391327326,0.239934501,-2.05929
Novel1,Novel1,MW7,MW,Medium,4,2,1,Target,101,29.96234131,30.19815254,30.46119118,30.20722834,29.96234131,30.46119118,0.498849869,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,30.45214,25.70248086,28.70086,26.22237753,1.993,0.30311588,1.431241263,0.211785311,-2.23933
Novel1,Novel1,MW8,MW,Medium,4,2,1,Target,101,29.43494034,NA,29.22647285,29.3307066,29.22647285,29.43494034,0.208467484,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.56851,24.59640071,28.70086,26.22237753,1.993,0.55356492,3.068964925,0.180375122,-2.47093
Novel1,Novel1,LQ3,LQ,Late,5,1,1,Target,128,28.5940609,28.67529106,28.83321953,28.70085716,28.5940609,28.83321953,0.23915863,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,28.70086,26.22237753,28.70086,26.22237753,1.993,1,1,1,0
Novel1,Novel1,LQ4,LQ,Late,5,1,1,Target,128,28.18196297,28.10457039,28.32246971,28.20300102,28.10457039,28.32246971,0.217899323,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,28.203,25.57874041,28.70086,26.22237753,1.993,1.40400527,1.558741203,0.900730196,-0.15083
Novel1,Novel1,LQ5,LQ,Late,5,1,1,Target,128,29.45370483,29.09542465,29.06522942,29.2047863,29.06522942,29.45370483,0.388475418,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.20479,26.77102083,28.70086,26.22237753,1.993,0.70930597,0.684979114,1.03551474,0.05035
Novel1,Novel1,LQ6,LQ,Late,5,2,1,Target,101,29.34432983,29.58346176,29.62094116,29.51624425,29.34432983,29.62094116,0.276611328,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.75555,26.04640784,28.70086,26.22237753,1.993,0.48730728,1.129026689,0.43161715,-1.21218
Novel1,Novel1,LQ7,LQ,Late,5,2,1,Target,101,27.94956589,28.00923538,27.82619286,27.92833138,27.82619286,28.00923538,0.183042526,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.15477,24.94259361,28.70086,26.22237753,1.993,1.45093044,2.417155701,0.600263539,-0.73633
Novel1,Novel1,LQ8,LQ,Late,5,2,1,Target,101,28.43688965,28.15149879,28.2317791,28.27338918,28.15149879,28.43688965,0.285390854,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.50262,25.43112605,28.70086,26.22237753,1.993,1.14466674,1.725780648,0.663274755,-0.59232
Novel1,Novel1,LW3,LW,Late,6,1,1,Target,128,29.0945816,28.95560074,28.77554321,28.94190852,28.77554321,29.0945816,0.319038391,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,28.94191,26.28165682,28.70086,26.22237753,1.993,0.8484909,0.959942941,0.883897227,-0.17805
Novel1,Novel1,LW4,LW,Late,6,1,1,Target,128,28.64810944,28.67015648,28.60848808,28.64225133,28.60848808,28.67015648,0.061668396,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,28.64225,25.96422745,28.70086,26.22237753,1.993,1.04075311,1.194862228,0.871023524,-0.19922
Novel1,Novel1,LW5,LW,Late,6,1,1,Target,128,28.99827766,29.345541,29.23271179,29.19217682,28.99827766,29.345541,0.347263336,ew8,28.69133568,28.60344315,28.45514297,28.58330727,28.45514297,28.69133568,0.236192703,28.58330727,1.977,29.19218,25.98699537,28.70086,26.22237753,1.993,0.71542828,1.176247401,0.608229427,-0.71731
Novel1,Novel1,LW6,LW,Late,6,2,1,Target,101,30.44787598,30.14489555,30.28164673,30.29147275,30.14489555,30.44787598,0.302980423,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,30.53707,25.6146532,28.70086,26.22237753,1.993,0.28606821,1.520610107,0.188127257,-2.41022
Novel1,Novel1,LW7,LW,Late,6,2,1,Target,101,28.27869415,28.25797844,28.10863495,28.21510251,28.10863495,28.27869415,0.170059204,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,28.44386,25.49879946,28.70086,26.22237753,1.993,1.19144003,1.647088425,0.723361306,-0.46721
Novel1,Novel1,LW8,LW,Late,6,2,1,Target,101,29.50792694,29.33071518,29.25491333,29.36451848,29.25491333,29.50792694,0.253013611,ew8,28.4170208,28.33353043,28.30972672,28.35342598,28.30972672,28.4170208,0.107294083,28.58330727,1.977,29.6026,25.38152329,28.70086,26.22237753,1.993,0.54085254,1.78583768,0.302856494,-1.72329
P17/29C-like protein DDB_G0287399,p17.29C,EQ3,EQ,Early,1,1,1,Target,109,40,40,40,40,40,40,0,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,41.12867,25.20347779,25.50903,26.22237753,1.993,0.00002287,2.019146866,1.13E-05,-16.42966
P17/29C-like protein DDB_G0287399,p17.29C,EQ4,EQ,Early,1,1,1,Target,109,NA,40,40,40,40,40,0,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,41.12867,25.71323328,25.50903,26.22237753,1.993,0.00002287,1.420667419,1.61E-05,-15.92248
P17/29C-like protein DDB_G0287399,p17.29C,EQ5,EQ,Early,1,1,1,Target,109,40,39.7858429,39.37300873,39.71961721,39.37300873,40,0.626991272,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,40.84037,25.41256952,25.50903,26.22237753,1.993,0.00002786,1.748007975,1.59E-05,-15.93709
P17/29C-like protein DDB_G0287399,p17.29C,EQ6,EQ,Early,1,2,1,Target,67,40,NA,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.43307804,25.50903,26.22237753,1.993,0.00004951,3.434853817,1.44E-05,-16.08222
P17/29C-like protein DDB_G0287399,p17.29C,EQ7,EQ,Early,1,2,1,Target,67,40,40,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,25.14049381,25.50903,26.22237753,1.993,0.00004951,2.108783874,2.35E-05,-15.37838
P17/29C-like protein DDB_G0287399,p17.29C,EQ8,EQ,Early,1,2,1,Target,67,40,40,39.45818329,39.81939443,39.45818329,40,0.541816711,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,39.81939,24.92400707,25.50903,26.22237753,1.993,0.00005602,2.448338325,2.29E-05,-15.41552
P17/29C-like protein DDB_G0287399,p17.29C,EW3,EW,Early,2,1,1,Target,109,40,40,40,40,40,40,0,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,41.12867,24.92159537,25.50903,26.22237753,1.993,0.00002287,2.452413802,9.33E-06,-16.71012
P17/29C-like protein DDB_G0287399,p17.29C,EW4,EW,Early,2,1,1,Target,109,40,NA,40,40,40,40,0,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,41.12867,24.71233996,25.50903,26.22237753,1.993,0.00002287,2.833135032,8.07E-06,-16.91831
P17/29C-like protein DDB_G0287399,p17.29C,EW5,EW,Early,2,1,1,Target,109,40,40,39.59983444,39.86661148,39.59983444,40,0.400165558,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,40.99151,25.76160569,25.50903,26.22237753,1.993,0.00002512,1.374056303,1.83E-05,-15.73899
P17/29C-like protein DDB_G0287399,p17.29C,EW6,EW,Early,2,2,1,Target,67,40,40,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.47244454,25.50903,26.22237753,1.993,0.00004951,3.342856272,1.48E-05,-16.04305
P17/29C-like protein DDB_G0287399,p17.29C,EW7,EW,Early,2,2,1,Target,67,40,40,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.71929254,25.50903,26.22237753,1.993,0.00004951,2.819583273,1.76E-05,-15.79745
P17/29C-like protein DDB_G0287399,p17.29C,EW8,EW,Early,2,2,1,Target,67,40,40,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.24981237,25.50903,26.22237753,1.993,0.00004951,3.897603549,1.27E-05,-16.26456
P17/29C-like protein DDB_G0287399,p17.29C,MQ3,MQ,Medium,3,1,1,Target,109,37.76021576,37.36593628,37.37937927,37.50184377,37.36593628,37.76021576,0.39427948,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,38.56002,25.89835631,25.50903,26.22237753,1.993,0.00013259,1.25039354,0.000106036,-13.20315
P17/29C-like protein DDB_G0287399,p17.29C,MQ4,MQ,Medium,3,1,1,Target,109,36.07408524,36.33766556,36.41482544,36.27552541,36.07408524,36.41482544,0.340740204,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,37.2991,25.52299076,25.50903,26.22237753,1.993,0.00031414,1.619837474,0.000193934,-12.33214
P17/29C-like protein DDB_G0287399,p17.29C,MQ5,MQ,Medium,3,1,1,Target,109,40,39.17267227,NA,39.58633614,39.17267227,40,0.827327728,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,40.70333,25.12014958,25.50903,26.22237753,1.993,0.0000306,2.138579087,1.43E-05,-16.09278
P17/29C-like protein DDB_G0287399,p17.29C,MQ6,MQ,Medium,3,2,1,Target,67,40,NA,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.7511079,25.50903,26.22237753,1.993,0.00004951,2.75839205,1.79E-05,-15.7658
P17/29C-like protein DDB_G0287399,p17.29C,MQ7,MQ,Medium,3,2,1,Target,67,NA,40,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.76237619,25.50903,26.22237753,1.993,0.00004951,2.737039456,1.81E-05,-15.75458
P17/29C-like protein DDB_G0287399,p17.29C,MQ8,MQ,Medium,3,2,1,Target,67,40,40,40,40,40,40,0,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,40,24.74641872,25.50903,26.22237753,1.993,0.00004951,2.767326721,1.79E-05,-15.77046
P17/29C-like protein DDB_G0287399,p17.29C,MW3,MW,Medium,4,1,1,Target,109,31.06063843,31.38168526,31.44715118,31.29649162,31.06063843,31.44715118,0.386512756,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,32.17957,25.29041589,25.50903,26.22237753,1.993,0.01042691,1.901644416,0.005483104,-7.51079
P17/29C-like protein DDB_G0287399,p17.29C,MW4,MW,Medium,4,1,1,Target,109,NA,34.28630066,33.91117096,34.09873581,33.91117096,34.28630066,0.3751297,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,35.06089,24.63093495,25.50903,26.22237753,1.993,0.00145247,2.996737273,0.000484683,-11.01067
P17/29C-like protein DDB_G0287399,p17.29C,MW5,MW,Medium,4,1,1,Target,109,32.21186829,32.72827148,32.69904327,32.54639435,32.21186829,32.72827148,0.516403198,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,33.46474,25.35247425,25.50903,26.22237753,1.993,0.00432839,1.821974874,0.002375661,-8.71746
P17/29C-like protein DDB_G0287399,p17.29C,MW6,MW,Medium,4,2,1,Target,67,33.72460938,34.19991684,33.7023468,33.87562434,33.7023468,34.19991684,0.497570038,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,33.87562,24.95817119,25.50903,26.22237753,1.993,0.00326778,2.391327326,0.001366514,-9.51528
P17/29C-like protein DDB_G0287399,p17.29C,MW7,MW,Medium,4,2,1,Target,67,31.96368217,31.5583992,31.47920609,31.66709582,31.47920609,31.96368217,0.484476089,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,31.6671,25.70248086,25.50903,26.22237753,1.993,0.01480522,1.431241263,0.010344322,-6.59502
P17/29C-like protein DDB_G0287399,p17.29C,MW8,MW,Medium,4,2,1,Target,67,33.25507355,NA,32.88374329,33.06940842,32.88374329,33.25507355,0.371330261,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,33.06941,24.59640071,25.50903,26.22237753,1.993,0.00567261,3.068964925,0.001848378,-9.07952
P17/29C-like protein DDB_G0287399,p17.29C,LQ3,LQ,Late,5,1,1,Target,109,24.80003548,24.84260941,24.78435135,24.80899874,24.78435135,24.84260941,0.058258057,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,25.50903,26.22237753,25.50903,26.22237753,1.993,1,1,1,0
P17/29C-like protein DDB_G0287399,p17.29C,LQ4,LQ,Late,5,1,1,Target,109,24.5723877,24.52814102,24.54966545,24.55006472,24.52814102,24.5723877,0.044246674,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,25.24279,25.57874041,25.50903,26.22237753,1.993,1.19977817,1.558741203,0.769709667,-0.37761
P17/29C-like protein DDB_G0287399,p17.29C,LQ5,LQ,Late,5,1,1,Target,109,24.54262924,24.33414841,24.38482666,24.42053477,24.33414841,24.54262924,0.208480835,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,25.1096,26.77102083,25.50903,26.22237753,1.993,1.31422784,0.684979114,1.918639291,0.94008
P17/29C-like protein DDB_G0287399,p17.29C,LQ6,LQ,Late,5,2,1,Target,67,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,26.33378,26.04640784,25.50903,26.22237753,1.993,0.5688034,1.129026689,0.503799784,-0.98908
P17/29C-like protein DDB_G0287399,p17.29C,LQ7,LQ,Late,5,2,1,Target,67,28.00056458,28.13371658,27.95699692,28.03042603,27.95699692,28.13371658,0.176719666,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,28.03043,24.94259361,25.50903,26.22237753,1.993,0.17818976,2.417155701,0.073718774,-3.76182
P17/29C-like protein DDB_G0287399,p17.29C,LQ8,LQ,Late,5,2,1,Target,67,24.2549057,24.16784668,24.11349678,24.17874972,24.11349678,24.2549057,0.14140892,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,24.17875,25.43112605,25.50903,26.22237753,1.993,2.48444749,1.725780648,1.43960792,0.52568
P17/29C-like protein DDB_G0287399,p17.29C,LW3,LW,Late,6,1,1,Target,109,23.69442368,23.64039421,23.66204262,23.66562017,23.64039421,23.69442368,0.054029465,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,24.33338,26.28165682,25.50903,26.22237753,1.993,2.23504978,0.959942941,2.328315237,1.21929
P17/29C-like protein DDB_G0287399,p17.29C,LW4,LW,Late,6,1,1,Target,109,25.65944099,25.6203537,25.661726,25.64717356,25.6203537,25.661726,0.041372299,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,26.37085,25.96422745,25.50903,26.22237753,1.993,0.55456002,1.194862228,0.46412047,-1.10743
P17/29C-like protein DDB_G0287399,p17.29C,LW5,LW,Late,6,1,1,Target,109,22.76906013,22.78360176,22.77729797,22.77665329,22.76906013,22.78360176,0.014541626,lq6,25.62688065,25.6199913,25.58649063,25.61112086,25.58649063,25.62688065,0.040390015,26.33378092,1.982,23.41933,25.98699537,25.50903,26.22237753,1.993,4.17690842,1.176247401,3.551045823,1.82824
P17/29C-like protein DDB_G0287399,p17.29C,LW6,LW,Late,6,2,1,Target,67,29.50849533,29.58605194,29.52743721,29.54066149,29.50849533,29.58605194,0.07755661,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,29.54066,25.6146532,25.50903,26.22237753,1.993,0.06341416,1.520610107,0.041703105,-4.5837
P17/29C-like protein DDB_G0287399,p17.29C,LW7,LW,Late,6,2,1,Target,67,24.17797852,24.19199371,24.26605606,24.21200943,24.17797852,24.26605606,0.088077545,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,24.21201,25.49879946,25.50903,26.22237753,1.993,2.42855666,1.647088425,1.474454332,0.56018
P17/29C-like protein DDB_G0287399,p17.29C,LW8,LW,Late,6,2,1,Target,67,24.0725708,24.16090584,24.09529686,24.10959117,24.0725708,24.16090584,0.088335037,lq6,26.33047295,26.4234848,26.24738503,26.33378092,26.24738503,26.4234848,0.176099777,26.33378092,1.982,24.10959,25.38152329,25.50903,26.22237753,1.993,2.60481618,1.78583768,1.458596271,0.54458
pancreatic lipase-related protein 2-like,plrp2,EQ3,EQ,Early,1,1,1,Target,107,29.94006348,29.79022408,29.83121109,29.85383288,29.79022408,29.94006348,0.149839401,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,31.77914,25.20347779,27.4859,26.22237753,1.993,0.08192171,2.019146866,0.040572438,-4.62336
pancreatic lipase-related protein 2-like,plrp2,EQ4,EQ,Early,1,1,1,Target,107,30.93658447,30.89717484,31.80890274,31.21422068,30.89717484,31.80890274,0.911727905,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,33.22726,25.71323328,27.4859,26.22237753,1.993,0.03522784,1.420667419,0.024796684,-5.33371
pancreatic lipase-related protein 2-like,plrp2,EQ5,EQ,Early,1,1,1,Target,107,27.184021,27.1894474,27.3369503,27.23680623,27.184021,27.3369503,0.152929306,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,28.99334,25.41256952,27.4859,26.22237753,1.993,0.41540637,1.748007975,0.237645581,-2.07312
pancreatic lipase-related protein 2-like,plrp2,EQ6,EQ,Early,1,2,1,Target,62,30.6581707,30.76387978,30.50537682,30.64247576,30.50537682,30.76387978,0.25850296,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,30.64248,24.43307804,27.4859,26.22237753,1.993,0.15888548,3.434853817,0.046256838,-4.43419
pancreatic lipase-related protein 2-like,plrp2,EQ7,EQ,Early,1,2,1,Target,62,29.97072411,29.80316353,29.52957726,29.76782163,29.52957726,29.97072411,0.441146851,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,29.76782,25.14049381,27.4859,26.22237753,1.993,0.26451804,2.108783874,0.125436299,-2.99497
pancreatic lipase-related protein 2-like,plrp2,EQ8,EQ,Early,1,2,1,Target,62,31.75392342,31.90766335,31.85250664,31.83803113,31.75392342,31.90766335,0.153739929,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,31.83803,24.92400707,27.4859,26.22237753,1.993,0.07915793,2.448338325,0.032331286,-4.95093
pancreatic lipase-related protein 2-like,plrp2,EW3,EW,Early,2,1,1,Target,107,29.89856911,29.7565155,29.51758957,29.72422473,29.51758957,29.89856911,0.380979538,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,31.64118,24.92159537,27.4859,26.22237753,1.993,0.08878054,2.452413802,0.036201289,-4.78782
pancreatic lipase-related protein 2-like,plrp2,EW4,EW,Early,2,1,1,Target,107,31.49752998,30.97646904,31.195961,31.22332001,30.97646904,31.49752998,0.521060944,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,33.23695,24.71233996,27.4859,26.22237753,1.993,0.03502955,2.833135032,0.012364234,-6.33768
pancreatic lipase-related protein 2-like,plrp2,EW5,EW,Early,2,1,1,Target,107,29.25846481,29.28753853,29.20938492,29.25179609,29.20938492,29.28753853,0.07815361,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,31.13828,25.76160569,27.4859,26.22237753,1.993,0.11901414,1.374056303,0.086615187,-3.52924
pancreatic lipase-related protein 2-like,plrp2,EW6,EW,Early,2,2,1,Target,62,32.22034836,32.52253342,32.46898651,32.4039561,32.22034836,32.52253342,0.302185059,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,32.40396,24.47244454,27.4859,26.22237753,1.993,0.05691953,3.342856272,0.017027214,-5.87601
pancreatic lipase-related protein 2-like,plrp2,EW7,EW,Early,2,2,1,Target,62,33.8240509,33.76291275,33.94200134,33.84298833,33.76291275,33.94200134,0.179088593,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,33.84299,24.71929254,27.4859,26.22237753,1.993,0.02460643,2.819583273,0.008726974,-6.8403
pancreatic lipase-related protein 2-like,plrp2,EW8,EW,Early,2,2,1,Target,62,33.65235519,32.94918442,32.9662323,33.1892573,32.94918442,33.65235519,0.703170776,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,33.18926,24.24981237,27.4859,26.22237753,1.993,0.03601679,3.897603549,0.009240753,-6.75777
pancreatic lipase-related protein 2-like,plrp2,MQ3,MQ,Medium,3,1,1,Target,107,23.97333908,24.27947044,24.27672768,24.1765124,23.97333908,24.27947044,0.306131363,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,25.73568,25.89835631,27.4859,26.22237753,1.993,2.77313784,1.25039354,2.217812037,1.14914
pancreatic lipase-related protein 2-like,plrp2,MQ4,MQ,Medium,3,1,1,Target,107,24.65542221,24.63486099,24.62573051,24.63867124,24.62573051,24.65542221,0.029691696,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,26.22765,25.52299076,27.4859,26.22237753,1.993,2.08188885,1.619837474,1.285245513,0.36204
pancreatic lipase-related protein 2-like,plrp2,MQ5,MQ,Medium,3,1,1,Target,107,24.9774189,25.36630249,25.05871201,25.13414447,24.9774189,25.36630249,0.388883591,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,26.75508,25.12014958,27.4859,26.22237753,1.993,1.53097508,2.138579087,0.71588425,-0.4822
pancreatic lipase-related protein 2-like,plrp2,MQ6,MQ,Medium,3,2,1,Target,62,24.82210922,24.75085449,24.81159592,24.79485321,24.75085449,24.82210922,0.07125473,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,24.79485,24.7511079,27.4859,26.22237753,1.993,4.79834707,2.75839205,1.73954499,0.79871
pancreatic lipase-related protein 2-like,plrp2,MQ7,MQ,Medium,3,2,1,Target,62,26.33414078,26.33599281,26.1270752,26.26573626,26.1270752,26.33599281,0.208917618,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,26.26574,24.76237619,27.4859,26.22237753,1.993,2.03618695,2.737039456,0.743937742,-0.42675
pancreatic lipase-related protein 2-like,plrp2,MQ8,MQ,Medium,3,2,1,Target,62,25.94652367,26.27581215,26.10267639,26.1083374,25.94652367,26.27581215,0.329288483,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,26.10834,24.74641872,27.4859,26.22237753,1.993,2.23179658,2.767326721,0.806481056,-0.31029
pancreatic lipase-related protein 2-like,plrp2,MW3,MW,Medium,4,1,1,Target,107,29.40669632,29.30700874,29.24738121,29.32036209,29.24738121,29.40669632,0.159315109,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,31.21127,25.29041589,27.4859,26.22237753,1.993,0.11405798,1.901644416,0.059978603,-4.05941
pancreatic lipase-related protein 2-like,plrp2,MW4,MW,Medium,4,1,1,Target,107,28.91125107,29.046175,29.06869125,29.00870577,28.91125107,29.06869125,0.157440186,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,30.87951,24.63093495,27.4859,26.22237753,1.993,0.13838579,2.996737273,0.046178819,-4.43662
pancreatic lipase-related protein 2-like,plrp2,MW5,MW,Medium,4,1,1,Target,107,30.04778671,29.85526085,NA,29.95152378,29.85526085,30.04778671,0.192525864,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,31.88313,25.35247425,27.4859,26.22237753,1.993,0.07710443,1.821974874,0.042319151,-4.56255
pancreatic lipase-related protein 2-like,plrp2,MW6,MW,Medium,4,2,1,Target,62,34.58460617,33.67692566,34.12309647,34.12820943,33.67692566,34.58460617,0.907680511,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,34.12821,24.95817119,27.4859,26.22237753,1.993,0.02083821,2.391327326,0.008714077,-6.84244
pancreatic lipase-related protein 2-like,plrp2,MW7,MW,Medium,4,2,1,Target,62,35.15827179,34.80487823,NA,34.98157501,34.80487823,35.15827179,0.353393555,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,34.98158,25.70248086,27.4859,26.22237753,1.993,0.01267294,1.431241263,0.008854509,-6.81937
pancreatic lipase-related protein 2-like,plrp2,MW8,MW,Medium,4,2,1,Target,62,31.97226143,32.26981354,32.65617752,32.2994175,31.97226143,32.65617752,0.683916092,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,32.29942,24.59640071,27.4859,26.22237753,1.993,0.06049501,3.068964925,0.019711861,-5.66479
pancreatic lipase-related protein 2-like,plrp2,LQ3,LQ,Late,5,1,1,Target,107,25.92534828,25.78577614,25.75094986,25.82069143,25.75094986,25.92534828,0.174398422,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,27.4859,26.22237753,27.4859,26.22237753,1.993,1,1,1,0
pancreatic lipase-related protein 2-like,plrp2,LQ4,LQ,Late,5,1,1,Target,107,27.79585266,27.95255852,27.94540024,27.89793714,27.79585266,27.95255852,0.156705856,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,29.69711,25.57874041,27.4859,26.22237753,1.993,0.27564656,1.558741203,0.176839207,-2.49949
pancreatic lipase-related protein 2-like,plrp2,LQ5,LQ,Late,5,1,1,Target,107,26.5132122,26.97921753,26.70077515,26.73106829,26.5132122,26.97921753,0.466005325,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,28.45499,26.77102083,27.4859,26.22237753,1.993,0.56849686,0.684979114,0.829947724,-0.26891
pancreatic lipase-related protein 2-like,plrp2,LQ6,LQ,Late,5,2,1,Target,62,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,30.12724,26.04640784,27.4859,26.22237753,1.993,0.21452973,1.129026689,0.190012986,-2.39583
pancreatic lipase-related protein 2-like,plrp2,LQ7,LQ,Late,5,2,1,Target,62,27.57475471,27.3762207,27.16704178,27.3726724,27.16704178,27.57475471,0.407712936,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,27.37267,24.94259361,27.4859,26.22237753,1.993,1.06821089,2.417155701,0.44192887,-1.17811
pancreatic lipase-related protein 2-like,plrp2,LQ8,LQ,Late,5,2,1,Target,62,30.94593811,30.62662697,30.62400055,30.73218854,30.62400055,30.94593811,0.321937561,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,30.73219,25.43112605,27.4859,26.22237753,1.993,0.150792,1.725780648,0.087376109,-3.51662
pancreatic lipase-related protein 2-like,plrp2,LW3,LW,Late,6,1,1,Target,107,33.42406845,32.72757339,33.01130295,33.05431493,32.72757339,33.42406845,0.696495056,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,35.18603,26.28165682,27.4859,26.22237753,1.993,0.01124945,0.959942941,0.011718879,-6.41502
pancreatic lipase-related protein 2-like,plrp2,LW4,LW,Late,6,1,1,Target,107,33.96800232,33.69799805,NA,33.83300018,33.69799805,33.96800232,0.270004272,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,36.01493,25.96422745,27.4859,26.22237753,1.993,0.00693968,1.194862228,0.005807936,-7.42776
pancreatic lipase-related protein 2-like,plrp2,LW5,LW,Late,6,1,1,Target,107,33.23567963,32.50676727,32.92673111,32.889726,32.50676727,33.23567963,0.728912354,lq6,28.42149353,28.22753334,28.25700951,28.30201213,28.22753334,28.42149353,0.19396019,30.12724241,1.791,35.01082,25.98699537,27.4859,26.22237753,1.993,0.01245876,1.176247401,0.010591953,-6.56089
pancreatic lipase-related protein 2-like,plrp2,LW6,LW,Late,6,2,1,Target,62,36.60736084,NA,35.6905632,36.14896202,35.6905632,36.60736084,0.916797638,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,36.14896,25.6146532,27.4859,26.22237753,1.993,0.00641825,1.520610107,0.004220841,-7.88825
pancreatic lipase-related protein 2-like,plrp2,LW7,LW,Late,6,2,1,Target,62,33.86731339,NA,34.26167297,34.06449318,33.86731339,34.26167297,0.394359589,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,34.06449,25.49879946,27.4859,26.22237753,1.993,0.02162652,1.647088425,0.013130153,-6.25097
pancreatic lipase-related protein 2-like,plrp2,LW8,LW,Late,6,2,1,Target,62,35.96867752,36.00766754,NA,35.98817253,35.96867752,36.00766754,0.038990021,lq6,30.43196678,29.94535065,30.00440979,30.12724241,29.94535065,30.43196678,0.486616135,30.12724241,1.791,35.98817,25.38152329,27.4859,26.22237753,1.993,0.00704875,1.78583768,0.003947026,-7.98502
putative uncharacterized protein DDB_G0271606,unch.G0271606,EQ3,EQ,Early,1,1,1,Target,108,NA,40,40,40,40,40,0,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,40,25.20347779,35.40148,26.22237753,1.993,0.06117516,2.019146866,0.03029753,-5.04466
putative uncharacterized protein DDB_G0271606,unch.G0271606,EQ4,EQ,Early,1,1,1,Target,108,NA,34.95376205,34.79792404,34.87584305,34.79792404,34.95376205,0.155838013,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,34.87584,25.71323328,35.40148,26.22237753,1.993,1.37625893,1.420667419,0.968741107,-0.04582
putative uncharacterized protein DDB_G0271606,unch.G0271606,EQ5,EQ,Early,1,1,1,Target,108,36.17298889,NA,36.73625183,36.45462036,36.17298889,36.73625183,0.563262939,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,36.45462,25.41256952,35.40148,26.22237753,1.993,0.52735608,1.748007975,0.301689742,-1.72886
putative uncharacterized protein DDB_G0271606,unch.G0271606,EQ6,EQ,Early,1,2,1,Target,70,34.05243683,33.88676834,33.35561371,33.76493963,33.35561371,34.05243683,0.69682312,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,34.21448,24.43307804,35.40148,26.22237753,1.993,2.0569165,3.434853817,0.598836693,-0.73977
putative uncharacterized protein DDB_G0271606,unch.G0271606,EQ7,EQ,Early,1,2,1,Target,70,33.46159744,33.97229004,33.57038116,33.66808955,33.46159744,33.97229004,0.510692596,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,34.11634,25.14049381,35.40148,26.22237753,1.993,2.18329791,2.108783874,1.035335073,0.0501
putative uncharacterized protein DDB_G0271606,unch.G0271606,EQ8,EQ,Early,1,2,1,Target,70,33.97843552,33.99346161,34.95211029,34.30800247,33.97843552,34.95211029,0.973674774,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,34.76477,24.92400707,35.40148,26.22237753,1.993,1.47234426,2.448338325,0.601364705,-0.73369
putative uncharacterized protein DDB_G0271606,unch.G0271606,EW3,EW,Early,2,1,1,Target,108,35.72597122,35.07698059,NA,35.40147591,35.07698059,35.72597122,0.648990631,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,35.40148,24.92159537,35.40148,26.22237753,1.993,1,2.452413802,0.407761528,-1.2942
putative uncharacterized protein DDB_G0271606,unch.G0271606,EW4,EW,Early,2,1,1,Target,108,40,NA,40,40,40,40,0,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,40,24.71233996,35.40148,26.22237753,1.993,0.06117516,2.833135032,0.021592745,-5.53331
putative uncharacterized protein DDB_G0271606,unch.G0271606,EW5,EW,Early,2,1,1,Target,108,40,40,40,40,40,40,0,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,40,25.76160569,35.40148,26.22237753,1.993,0.06117516,1.374056303,0.044521584,-4.48935
putative uncharacterized protein DDB_G0271606,unch.G0271606,EW6,EW,Early,2,2,1,Target,70,26.91562271,26.87051964,26.76812172,26.85142136,26.76812172,26.91562271,0.147500992,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,27.20891,24.47244454,35.40148,26.22237753,1.993,145.1420684,3.342856272,43.41857877,5.44024
putative uncharacterized protein DDB_G0271606,unch.G0271606,EW7,EW,Early,2,2,1,Target,70,34.42337799,34.67292786,35.22082901,34.77237829,34.42337799,35.22082901,0.797451019,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,35.23533,24.71929254,35.40148,26.22237753,1.993,1.10622092,2.819583273,0.392334898,-1.34984
putative uncharacterized protein DDB_G0271606,unch.G0271606,EW8,EW,Early,2,2,1,Target,70,26.32988548,26.21231842,26.34430504,26.29550298,26.21231842,26.34430504,0.131986618,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,26.64559,24.24981237,35.40148,26.22237753,1.993,204.3798095,3.897603549,52.43730074,5.71252
putative uncharacterized protein DDB_G0271606,unch.G0271606,MQ3,MQ,Medium,3,1,1,Target,108,35.09823608,35.39440918,NA,35.24632263,35.09823608,35.39440918,0.296173096,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,35.24632,25.89835631,35.40148,26.22237753,1.993,1.09885581,1.25039354,0.878807973,-0.18638
putative uncharacterized protein DDB_G0271606,unch.G0271606,MQ4,MQ,Medium,3,1,1,Target,108,36.10544586,NA,36.56730652,36.33637619,36.10544586,36.56730652,0.461860657,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,36.33638,25.52299076,35.40148,26.22237753,1.993,0.56663756,1.619837474,0.349811367,-1.51535
putative uncharacterized protein DDB_G0271606,unch.G0271606,MQ5,MQ,Medium,3,1,1,Target,108,25.49248695,25.43110847,25.46072769,25.46144104,25.43110847,25.49248695,0.061378479,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,25.46144,25.12014958,35.40148,26.22237753,1.993,419.6658288,2.138579087,196.2358238,7.61644
putative uncharacterized protein DDB_G0271606,unch.G0271606,MQ6,MQ,Medium,3,2,1,Target,70,34.37551117,NA,34.95735931,34.66643524,34.37551117,34.95735931,0.581848145,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,35.12797,24.7511079,35.40148,26.22237753,1.993,1.18078148,2.75839205,0.428068765,-1.22409
putative uncharacterized protein DDB_G0271606,unch.G0271606,MQ7,MQ,Medium,3,2,1,Target,70,25.39226723,25.33432198,25.35524368,25.36061096,25.33432198,25.39226723,0.057945251,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,25.69826,24.76237619,35.40148,26.22237753,1.993,363.425038,2.737039456,132.7803431,7.0529
putative uncharacterized protein DDB_G0271606,unch.G0271606,MQ8,MQ,Medium,3,2,1,Target,70,29.73832893,29.75028801,29.69672585,29.7284476,29.69672585,29.75028801,0.053562164,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,30.12424,24.74641872,35.40148,26.22237753,1.993,24.68978634,2.767326721,8.921890627,3.15735
putative uncharacterized protein DDB_G0271606,unch.G0271606,MW3,MW,Medium,4,1,1,Target,108,34.40450668,NA,34.70655823,34.55553246,34.40450668,34.70655823,0.302051544,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,34.55553,25.29041589,35.40148,26.22237753,1.993,1.6719426,1.901644416,0.879208849,-0.18572
putative uncharacterized protein DDB_G0271606,unch.G0271606,MW4,MW,Medium,4,1,1,Target,108,NA,34.44362259,35.42546082,34.9345417,34.44362259,35.42546082,0.981838226,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,34.93454,24.63093495,35.40148,26.22237753,1.993,1.32804007,2.996737273,0.443161995,-1.17409
putative uncharacterized protein DDB_G0271606,unch.G0271606,MW5,MW,Medium,4,1,1,Target,108,33.93773651,33.69471359,34.2388382,33.9570961,33.69471359,34.2388382,0.544124603,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,33.9571,25.35247425,35.40148,26.22237753,1.993,2.40509572,1.821974874,1.320048787,0.40059
putative uncharacterized protein DDB_G0271606,unch.G0271606,MW6,MW,Medium,4,2,1,Target,70,34.45319748,34.31275558,35.1967392,34.65423075,34.31275558,35.1967392,0.883983612,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,35.11561,24.95817119,35.40148,26.22237753,1.993,1.18968734,2.391327326,0.497500833,-1.00723
putative uncharacterized protein DDB_G0271606,unch.G0271606,MW7,MW,Medium,4,2,1,Target,70,34.73561859,33.85109711,34.778759,34.45515823,33.85109711,34.778759,0.927661896,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,34.91388,25.70248086,35.40148,26.22237753,1.993,1.34481333,1.431241263,0.939613301,-0.08986
putative uncharacterized protein DDB_G0271606,unch.G0271606,MW8,MW,Medium,4,2,1,Target,70,32.67571259,33.28761673,NA,32.98166466,32.67571259,33.28761673,0.611904144,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,33.42077,24.59640071,35.40148,26.22237753,1.993,3.3316032,3.068964925,1.085578779,0.11846
putative uncharacterized protein DDB_G0271606,unch.G0271606,LQ3,LQ,Late,5,1,1,Target,108,35.54227448,35.54550552,35.74685287,35.61154429,35.54227448,35.74685287,0.2045784,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,35.61154,26.22237753,35.40148,26.22237753,1.993,0.88017432,1,0.880174319,-0.18414
putative uncharacterized protein DDB_G0271606,unch.G0271606,LQ4,LQ,Late,5,1,1,Target,108,37.13014603,37.05053711,NA,37.09034157,37.05053711,37.13014603,0.079608917,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,37.09034,25.57874041,35.40148,26.22237753,1.993,0.35838886,1.558741203,0.229921975,-2.12078
putative uncharacterized protein DDB_G0271606,unch.G0271606,LQ5,LQ,Late,5,1,1,Target,108,36.87276077,35.66982269,NA,36.27129173,35.66982269,36.87276077,1.20293808,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,36.27129,26.77102083,35.40148,26.22237753,1.993,0.58949397,0.684979114,0.860601384,-0.21658
putative uncharacterized protein DDB_G0271606,unch.G0271606,LQ6,LQ,Late,5,2,1,Target,70,32.58878326,32.31523514,32.9633522,32.62245687,32.31523514,32.9633522,0.648117065,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,33.05678,26.04640784,35.40148,26.22237753,1.993,4.15623728,1.129026689,3.681256895,1.8802
putative uncharacterized protein DDB_G0271606,unch.G0271606,LQ7,LQ,Late,5,2,1,Target,70,30.6366024,31.01807785,30.43519974,30.69662666,30.43519974,31.01807785,0.582878113,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,31.10531,24.94259361,35.40148,26.22237753,1.993,13.60316584,2.417155701,5.627757384,2.49256
putative uncharacterized protein DDB_G0271606,unch.G0271606,LQ8,LQ,Late,5,2,1,Target,70,34.02661896,NA,33.70149231,33.86405563,33.70149231,34.02661896,0.325126648,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,34.31491,25.43112605,35.40148,26.22237753,1.993,1.93514917,1.725780648,1.121318156,0.1652
putative uncharacterized protein DDB_G0271606,unch.G0271606,LW3,LW,Late,6,1,1,Target,108,NA,34.12520218,34.53728867,34.33124542,34.12520218,34.53728867,0.412086487,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,34.33125,26.28165682,35.40148,26.22237753,1.993,1.91604005,0.959942941,1.995993686,0.99711
putative uncharacterized protein DDB_G0271606,unch.G0271606,LW4,LW,Late,6,1,1,Target,108,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,25.83436,25.96422745,35.40148,26.22237753,1.993,334.5799781,1.194862228,280.0155284,8.12936
putative uncharacterized protein DDB_G0271606,unch.G0271606,LW5,LW,Late,6,1,1,Target,108,25.95375252,25.88402939,25.91108131,25.91628774,25.88402939,25.95375252,0.069723129,mq7,26.2063961,26.12724876,26.20129585,26.17831357,26.12724876,26.2063961,0.079147339,26.17831357,1.836,25.91629,25.98699537,35.40148,26.22237753,1.993,318.3332418,1.176247401,270.6345973,8.0802
putative uncharacterized protein DDB_G0271606,unch.G0271606,LW6,LW,Late,6,2,1,Target,70,27.74509621,27.66337204,27.82838058,27.74561628,27.66337204,27.82838058,0.165008545,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,28.11501,25.6146532,35.40148,26.22237753,1.993,83.69476212,1.520610107,55.04025109,5.78242
putative uncharacterized protein DDB_G0271606,unch.G0271606,LW7,LW,Late,6,2,1,Target,70,28.931427,28.87220764,28.84319687,28.88227717,28.84319687,28.931427,0.088230133,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,29.26681,25.49879946,35.40148,26.22237753,1.993,41.56916411,1.647088425,25.23796748,4.65752
putative uncharacterized protein DDB_G0271606,unch.G0271606,LW8,LW,Late,6,2,1,Target,70,31.94409561,32.07461548,31.99014091,32.00295067,31.94409561,32.07461548,0.130519867,mq7,25.72631645,25.91050529,25.86626434,25.83436203,25.72631645,25.91050529,0.184188843,26.17831357,1.836,32.42903,25.38152329,35.40148,26.22237753,1.993,6.08621791,1.78583768,3.40804653,1.76895
Takeout,Takeout,EQ3,EQ,Early,1,1,1,Target,111,24.07610512,24.05616188,24.02879524,24.05368741,24.02879524,24.07610512,0.047309875,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,24.05369,25.20347779,30.38801,26.22237753,1.993,67.84587566,2.019146866,33.60125844,5.07044
Takeout,Takeout,EQ4,EQ,Early,1,1,1,Target,111,23.65994453,23.9507885,23.57883644,23.72985649,23.57883644,23.9507885,0.371952057,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,23.72986,25.71323328,30.38801,26.22237753,1.993,84.16990624,1.420667419,59.24673511,5.88866
Takeout,Takeout,EQ5,EQ,Early,1,1,1,Target,111,29.92913055,29.91212463,29.65251732,29.8312575,29.65251732,29.92913055,0.276613235,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,29.83126,25.41256952,30.38801,26.22237753,1.993,1.44870664,1.748007975,0.828775763,-0.27095
Takeout,Takeout,EQ6,EQ,Early,1,2,1,Target,30,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,26.25342,24.43307804,30.38801,26.22237753,1.993,15.68510711,3.434853817,4.566455501,2.19107
Takeout,Takeout,EQ7,EQ,Early,1,2,1,Target,30,24.08180237,24.00608444,23.92531776,24.00440152,23.92531776,24.08180237,0.156484604,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,24.76287,25.14049381,30.38801,26.22237753,1.993,42.31259015,2.108783874,20.06492494,4.3266
Takeout,Takeout,EQ8,EQ,Early,1,2,1,Target,30,25.49451256,25.33481789,25.29070663,25.37334569,25.29070663,25.49451256,0.203805923,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,26.17506,24.92400707,30.38801,26.22237753,1.993,16.52504009,2.448338325,6.749492063,2.75478
Takeout,Takeout,EW3,EW,Early,2,1,1,Target,111,31.04732132,30.55207062,30.90287209,30.83408801,30.55207062,31.04732132,0.495250702,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,30.83409,24.92159537,30.38801,26.22237753,1.993,0.74305198,2.452413802,0.302988011,-1.72267
Takeout,Takeout,EW4,EW,Early,2,1,1,Target,111,24.35485077,24.38161278,24.2753315,24.33726501,24.2753315,24.38161278,0.106281281,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,24.33727,24.71233996,30.38801,26.22237753,1.993,56.1731813,2.833135032,19.82721638,4.30941
Takeout,Takeout,EW5,EW,Early,2,1,1,Target,111,23.94415092,24.08091545,23.94436455,23.98981031,23.94415092,24.08091545,0.136764526,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,23.98981,25.76160569,30.38801,26.22237753,1.993,70.79344709,1.374056303,51.52150383,5.6871
Takeout,Takeout,EW6,EW,Early,2,2,1,Target,30,24.08661652,23.90739059,23.97562981,23.98987897,23.90739059,24.08661652,0.179225922,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,24.74788,24.47244454,30.38801,26.22237753,1.993,42.7367391,3.342856272,12.78449793,3.67632
Takeout,Takeout,EW7,EW,Early,2,2,1,Target,30,25.20582008,25.23754311,25.04992485,25.16442935,25.04992485,25.23754311,0.187618256,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,25.95955,24.71929254,30.38801,26.22237753,1.993,19.07470583,2.819583273,6.765079794,2.75811
Takeout,Takeout,EW8,EW,Early,2,2,1,Target,30,22.94506454,22.8886528,22.98787308,22.94053014,22.8886528,22.98787308,0.099220276,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,23.66538,24.24981237,30.38801,26.22237753,1.993,87.8617688,3.897603549,22.54251047,4.49458
Takeout,Takeout,MQ3,MQ,Medium,3,1,1,Target,111,25.71823311,25.64544106,25.90685463,25.75684293,25.64544106,25.90685463,0.261413574,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,25.75684,25.89835631,30.38801,26.22237753,1.993,21.83072886,1.25039354,17.45908641,4.12591
Takeout,Takeout,MQ4,MQ,Medium,3,1,1,Target,111,25.52759171,25.64836502,25.44029427,25.53875033,25.44029427,25.64836502,0.208070755,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,25.53875,25.52299076,30.38801,26.22237753,1.993,25.24225494,1.619837474,15.58320223,3.96192
Takeout,Takeout,MQ5,MQ,Medium,3,1,1,Target,111,26.44077873,26.33103752,26.30272865,26.35818164,26.30272865,26.44077873,0.138050079,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,26.35818,25.12014958,30.38801,26.22237753,1.993,14.62835386,2.138579087,6.840221132,2.77404
Takeout,Takeout,MQ6,MQ,Medium,3,2,1,Target,30,25.15771484,25.5080204,25.3429184,25.33621788,25.15771484,25.5080204,0.350305557,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,26.13676,24.7511079,30.38801,26.22237753,1.993,16.95184462,2.75839205,6.145553031,2.61954
Takeout,Takeout,MQ7,MQ,Medium,3,2,1,Target,30,26.24522591,26.28288841,26.25796127,26.2620252,26.24522591,26.28288841,0.037662506,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,27.09182,24.76237619,30.38801,26.22237753,1.993,8.97569782,2.737039456,3.279345426,1.71341
Takeout,Takeout,MQ8,MQ,Medium,3,2,1,Target,30,28.32736397,27.94505692,27.76078224,28.01106771,27.76078224,28.32736397,0.566581726,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,28.89613,24.74641872,30.38801,26.22237753,1.993,2.70001289,2.767326721,0.975675503,-0.03553
Takeout,Takeout,MW3,MW,Medium,4,1,1,Target,111,34.67337036,35.22759628,NA,34.95048332,34.67337036,35.22759628,0.554225922,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,34.95048,25.29041589,30.38801,26.22237753,1.993,0.04795048,1.901644416,0.025215273,-5.30956
Takeout,Takeout,MW4,MW,Medium,4,1,1,Target,111,27.59459686,27.5298214,27.47597504,27.53346443,27.47597504,27.59459686,0.118621826,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,27.53346,24.63093495,30.38801,26.22237753,1.993,6.6891488,2.996737273,2.23214389,1.15843
Takeout,Takeout,MW5,MW,Medium,4,1,1,Target,111,27.61826324,27.63729668,27.70756149,27.6543738,27.61826324,27.70756149,0.089298248,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,27.65437,25.35247425,30.38801,26.22237753,1.993,6.17178484,1.821974874,3.387414905,1.76018
Takeout,Takeout,MW6,MW,Medium,4,2,1,Target,30,28.12052727,27.92637062,27.9241581,27.99035199,27.9241581,28.12052727,0.196369171,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,28.87476,24.95817119,30.38801,26.22237753,1.993,2.73870274,2.391327326,1.145264688,0.19568
Takeout,Takeout,MW7,MW,Medium,4,2,1,Target,30,29.42882156,29.31240654,29.3497715,29.36366653,29.31240654,29.42882156,0.116415024,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,30.29147,25.70248086,30.38801,26.22237753,1.993,1.06638547,1.431241263,0.745077364,-0.42454
Takeout,Takeout,MW8,MW,Medium,4,2,1,Target,30,25.97312927,25.90720177,25.93549728,25.93860944,25.90720177,25.97312927,0.065927505,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,26.75819,24.59640071,30.38801,26.22237753,1.993,11.20821111,3.068964925,3.652114436,1.86873
Takeout,Takeout,LQ3,LQ,Late,5,1,1,Target,111,30.25567436,30.6843853,30.22396469,30.38800812,30.22396469,30.6843853,0.460420609,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,30.38801,26.22237753,30.38801,26.22237753,1.993,1,1,1,0
Takeout,Takeout,LQ4,LQ,Late,5,1,1,Target,111,32.54963684,32.99985886,32.23540115,32.59496562,32.23540115,32.99985886,0.764457703,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,32.59497,25.57874041,30.38801,26.22237753,1.993,0.23007745,1.558741203,0.147604648,-2.76019
Takeout,Takeout,LQ5,LQ,Late,5,1,1,Target,111,31.21447754,31.17439079,31.07109642,31.15332158,31.07109642,31.21447754,0.143381119,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,31.15332,26.77102083,30.38801,26.22237753,1.993,0.60077948,0.684979114,0.877077085,-0.18922
Takeout,Takeout,LQ6,LQ,Late,5,2,1,Target,30,30.55867767,30.98149681,30.93977165,30.82664871,30.55867767,30.98149681,0.422819138,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,31.80067,26.04640784,30.38801,26.22237753,1.993,0.39042496,1.129026689,0.345806674,-1.53196
Takeout,Takeout,LQ7,LQ,Late,5,2,1,Target,30,28.0724659,27.921978,27.94652176,27.98032188,27.921978,28.0724659,0.1504879,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,28.86441,24.94259361,30.38801,26.22237753,1.993,2.75763426,2.417155701,1.140859178,0.19012
Takeout,Takeout,LQ8,LQ,Late,5,2,1,Target,30,NA,32.77684784,32.28572083,32.53128433,32.28572083,32.77684784,0.491127014,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,33.55917,25.43112605,30.38801,26.22237753,1.993,0.12108235,1.725780648,0.070160914,-3.83319
Takeout,Takeout,LW3,LW,Late,6,1,1,Target,111,34.564888,34.34884262,33.9713974,34.29504267,33.9713974,34.564888,0.593490601,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,34.29504,26.28165682,30.38801,26.22237753,1.993,0.07418378,0.959942941,0.077279362,-3.69377
Takeout,Takeout,LW4,LW,Late,6,1,1,Target,111,NA,35.39223099,35.39956665,35.39589882,35.39223099,35.39956665,0.007335663,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,35.3959,25.96422745,30.38801,26.22237753,1.993,0.03564546,1.194862228,0.02983228,-5.06698
Takeout,Takeout,LW5,LW,Late,6,1,1,Target,111,29.20304489,29.18982124,29.0049057,29.13259061,29.0049057,29.20304489,0.198139191,eq6,26.22455788,26.15439224,26.38129997,26.2534167,26.15439224,26.38129997,0.22690773,26.2534167,1.946,29.13259,25.98699537,30.38801,26.22237753,1.993,2.30672051,1.176247401,1.961084472,0.97165
Takeout,Takeout,LW6,LW,Late,6,2,1,Target,30,34.64723969,NA,34.14969635,34.39846802,34.14969635,34.64723969,0.497543335,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,35.48535,25.6146532,30.38801,26.22237753,1.993,0.03358454,1.520610107,0.022086228,-5.50071
Takeout,Takeout,LW7,LW,Late,6,2,1,Target,30,34.08621216,NA,34.47887421,34.28254318,34.08621216,34.47887421,0.392662048,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,35.36576,25.49879946,30.38801,26.22237753,1.993,0.03636783,1.647088425,0.02208007,-5.50111
Takeout,Takeout,LW8,LW,Late,6,2,1,Target,30,34.44932556,34.09794998,34.92835999,34.49187851,34.09794998,34.92835999,0.830410004,eq6,25.49835968,25.50338364,25.34615135,25.44929822,25.34615135,25.50338364,0.157232285,26.2534167,1.946,35.58171,25.38152329,30.38801,26.22237753,1.993,0.03149757,1.78583768,0.017637421,-5.82522
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EQ3,EQ,Early,1,1,1,Target,123,28.6014328,28.70879364,28.55858612,28.62293752,28.55858612,28.70879364,0.15020752,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.62294,25.20347779,30.19988,26.22237753,1.993,2.75840253,2.019146866,1.366122779,0.45009
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EQ4,EQ,Early,1,1,1,Target,123,28.9482193,29.00561142,28.95704269,28.97029114,28.9482193,29.00561142,0.05739212,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.97029,25.71323328,30.19988,26.22237753,1.993,2.20593935,1.420667419,1.552748604,0.63482
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EQ5,EQ,Early,1,1,1,Target,123,28.13858795,28.37014198,28.18920708,28.23264567,28.13858795,28.37014198,0.231554031,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.23265,25.41256952,30.19988,26.22237753,1.993,3.54584971,1.748007975,2.028508884,1.02042
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EQ6,EQ,Early,1,2,1,Target,22,27.88492584,27.74898148,27.874403,27.83610344,27.74898148,27.88492584,0.135944366,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.62811,24.43307804,30.19988,26.22237753,1.993,2.74923348,3.434853817,0.800393153,-0.32122
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EQ7,EQ,Early,1,2,1,Target,22,28.45560646,28.5626564,28.51124573,28.5098362,28.45560646,28.5626564,0.107049942,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.32101,25.14049381,30.19988,26.22237753,1.993,1.76030481,2.108783874,0.834748801,-0.26059
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EQ8,EQ,Early,1,2,1,Target,22,27.82684135,27.95485687,27.80119133,27.86096319,27.80119133,27.95485687,0.153665543,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.65368,24.92400707,30.19988,26.22237753,1.993,2.70437674,2.448338325,1.104576404,0.14349
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EW3,EW,Early,2,1,1,Target,123,28.35087013,28.60048676,28.42769814,28.45968501,28.35087013,28.60048676,0.249616623,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.45969,24.92159537,30.19988,26.22237753,1.993,3.0639151,2.452413802,1.249346702,0.32117
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EW4,EW,Early,2,1,1,Target,123,28.6805687,28.9895401,28.97146797,28.88052559,28.6805687,28.9895401,0.308971405,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.88053,24.71233996,30.19988,26.22237753,1.993,2.33710134,2.833135032,0.824917031,-0.27768
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EW5,EW,Early,2,1,1,Target,123,28.79049683,28.91910744,28.7161274,28.80857722,28.7161274,28.91910744,0.202980042,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.80858,25.76160569,30.19988,26.22237753,1.993,2.44783821,1.374056303,1.781468636,0.83307
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EW6,EW,Early,2,2,1,Target,22,27.62192726,27.70663071,27.73731422,27.68862406,27.62192726,27.73731422,0.115386963,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.47644,24.47244454,30.19988,26.22237753,1.993,3.03106782,3.342856272,0.906729926,-0.14126
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EW7,EW,Early,2,2,1,Target,22,27.95337677,28.05502129,28.05893517,28.02244441,27.95337677,28.05893517,0.105558395,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.81976,24.71929254,30.19988,26.22237753,1.993,2.43029604,2.819583273,0.861934477,-0.21435
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,EW8,EW,Early,2,2,1,Target,22,27.68656731,27.82829094,27.6459465,27.72026825,27.6459465,27.82829094,0.182344437,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.50898,24.24981237,30.19988,26.22237753,1.993,2.96825662,3.897603549,0.761559401,-0.39297
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MQ3,MQ,Medium,3,1,1,Target,123,29.26251793,29.27577972,29.42864227,29.32231331,29.26251793,29.42864227,0.166124344,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.32231,25.89835631,30.19988,26.22237753,1.993,1.75883427,1.25039354,1.406624566,0.49224
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MQ4,MQ,Medium,3,1,1,Target,123,29.75576782,29.6662178,29.73276138,29.718249,29.6662178,29.75576782,0.089550018,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.71825,25.52299076,30.19988,26.22237753,1.993,1.3632802,1.619837474,0.841615425,-0.24877
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MQ5,MQ,Medium,3,1,1,Target,123,29.65295029,29.29558754,29.49570847,29.48141543,29.29558754,29.65295029,0.357362747,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.48142,25.12014958,30.19988,26.22237753,1.993,1.58768965,2.138579087,0.742403994,-0.42972
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MQ6,MQ,Medium,3,2,1,Target,22,27.67863655,27.62573624,27.5476532,27.617342,27.5476532,27.67863655,0.130983353,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.40313,24.7511079,30.19988,26.22237753,1.993,3.17746938,2.75839205,1.151928126,0.20405
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MQ7,MQ,Medium,3,2,1,Target,22,28.13877106,28.24881172,28.15006065,28.17921448,28.13877106,28.24881172,0.110040665,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.98099,24.76237619,30.19988,26.22237753,1.993,2.19081194,2.737039456,0.800431261,-0.32115
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MQ8,MQ,Medium,3,2,1,Target,22,28.24473,28.32196426,28.20976639,28.25882022,28.20976639,28.32196426,0.112197876,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.06286,24.74641872,30.19988,26.22237753,1.993,2.07839095,2.767326721,0.751046465,-0.41303
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MW3,MW,Medium,4,1,1,Target,123,29.06729507,29.32739067,29.36234093,29.25234222,29.06729507,29.36234093,0.295045853,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.25234,25.29041589,30.19988,26.22237753,1.993,1.8398294,1.901644416,0.967493914,-0.04768
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MW4,MW,Medium,4,1,1,Target,123,28.88662529,28.703022,29.1055336,28.89839363,28.703022,29.1055336,0.402511597,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.89839,24.63093495,30.19988,26.22237753,1.993,2.31038587,2.996737273,0.770967109,-0.37526
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MW5,MW,Medium,4,1,1,Target,123,29.35879135,29.32422638,29.20778465,29.29693413,29.20778465,29.35879135,0.151006699,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.29693,25.35247425,30.19988,26.22237753,1.993,1.78779141,1.821974874,0.981238236,-0.02732
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MW6,MW,Medium,4,2,1,Target,22,28.43301201,28.38591957,28.6189003,28.47927729,28.38591957,28.6189003,0.232980728,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.28959,24.95817119,30.19988,26.22237753,1.993,1.79626408,2.391327326,0.75115776,-0.41281
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MW7,MW,Medium,4,2,1,Target,22,29.00188255,28.88304901,29.02052879,28.96848679,28.88304901,29.02052879,0.137479782,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.79271,25.70248086,30.19988,26.22237753,1.993,1.29950063,1.431241263,0.907953579,-0.13931
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,MW8,MW,Medium,4,2,1,Target,22,28.00435257,27.99392128,28.19446564,28.0642465,27.99392128,28.19446564,0.200544357,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,28.86275,24.59640071,30.19988,26.22237753,1.993,2.36399033,3.068964925,0.770289132,-0.37653
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LQ3,LQ,Late,5,1,1,Target,123,30.35534859,30.00671005,30.2375679,30.19987551,30.00671005,30.35534859,0.348638535,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,30.19988,26.22237753,30.19988,26.22237753,1.993,1,1,1,0
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LQ4,LQ,Late,5,1,1,Target,123,28.87875557,28.96676064,29.08296394,28.97616005,28.87875557,29.08296394,0.204208374,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,28.97616,25.57874041,30.19988,26.22237753,1.993,2.1976249,1.558741203,1.409871563,0.49556
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LQ5,LQ,Late,5,1,1,Target,123,31.39684677,31.13653755,31.18632507,31.23990313,31.13653755,31.39684677,0.260309219,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,31.2399,26.77102083,30.19988,26.22237753,1.993,0.51212496,0.684979114,0.747650475,-0.41956
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LQ6,LQ,Late,5,2,1,Target,22,29.30273247,30.02470589,29.70923233,29.67889023,29.30273247,30.02470589,0.721973419,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,30.52333,26.04640784,30.19988,26.22237753,1.993,0.81210835,1.129026689,0.719299512,-0.47534
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LQ7,LQ,Late,5,2,1,Target,22,28.20962143,28.22476959,28.22494125,28.21977743,28.20962143,28.22494125,0.015319824,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.0227,24.94259361,30.19988,26.22237753,1.993,2.13278822,2.417155701,0.882354503,-0.18057
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LQ8,LQ,Late,5,2,1,Target,22,28.6051445,28.69741631,29.0939827,28.79884783,28.6051445,29.0939827,0.488838196,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.61825,25.43112605,30.19988,26.22237753,1.993,1.45388131,1.725780648,0.842448495,-0.24734
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LW3,LW,Late,6,1,1,Target,123,29.58460999,29.75499344,29.78768921,29.70909754,29.58460999,29.78768921,0.203079224,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.7091,26.28165682,30.19988,26.22237753,1.993,1.37133134,0.959942941,1.428555053,0.51456
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LW4,LW,Late,6,1,1,Target,123,28.91005325,29.1772747,29.25733757,29.11488851,28.91005325,29.25733757,0.347284317,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.11489,25.96422745,30.19988,26.22237753,1.993,2.00996018,1.194862228,1.682168985,0.75032
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LW5,LW,Late,6,1,1,Target,123,29.7652359,29.65978241,29.72977638,29.7182649,29.65978241,29.7652359,0.105453491,lw7,29.5663414,29.26194382,29.52311325,29.45046616,29.26194382,29.5663414,0.304397583,29.45046616,1.903,29.71826,25.98699537,30.19988,26.22237753,1.993,1.36326626,1.176247401,1.158996195,0.21288
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LW6,LW,Late,6,2,1,Target,22,28.67745209,28.8692131,28.98265266,28.84310595,28.67745209,28.98265266,0.305200577,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.66377,25.6146532,30.19988,26.22237753,1.993,1.41191852,1.520610107,0.928521065,-0.10699
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LW7,LW,Late,6,2,1,Target,22,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.45047,25.49879946,30.19988,26.22237753,1.993,1.6196234,1.647088425,0.983325106,-0.02426
uncharacterized LOC100648272%2C transcript variant X2 (dnmt3),dnmt3,LW8,LW,Late,6,2,1,Target,22,28.69783974,28.63710785,28.61913872,28.6513621,28.61913872,28.69783974,0.078701019,lw7,28.6186924,28.58201981,28.70640755,28.63570658,28.58201981,28.70640755,0.124387741,29.45046616,1.903,29.46657,25.38152329,30.19988,26.22237753,1.993,1.60293094,1.78583768,0.897579301,-0.15589
VHDL,VHDL,EQ3,EQ,Early,1,1,1,Target,124,40,40,NA,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.20347779,27.84772,26.22237753,1.993,0.00035944,2.019146866,0.000178016,-12.45571
VHDL,VHDL,EQ4,EQ,Early,1,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.71323328,27.84772,26.22237753,1.993,0.00035944,1.420667419,0.000253008,-11.94853
VHDL,VHDL,EQ5,EQ,Early,1,1,1,Target,124,39.4877243,40,40,39.82924143,39.4877243,40,0.512275696,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,39.8543,25.41256952,27.84772,26.22237753,1.993,0.00040175,1.748007975,0.000229833,-12.08713
VHDL,VHDL,EQ6,EQ,Early,1,2,1,Target,39,40,40,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.43307804,27.84772,26.22237753,1.993,0.00036538,3.434853817,0.000106374,-13.19857
VHDL,VHDL,EQ7,EQ,Early,1,2,1,Target,39,40,40,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,25.14049381,27.84772,26.22237753,1.993,0.00036538,2.108783874,0.000173265,-12.49473
VHDL,VHDL,EQ8,EQ,Early,1,2,1,Target,39,40,NA,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.92400707,27.84772,26.22237753,1.993,0.00036538,2.448338325,0.000149236,-12.71012
VHDL,VHDL,EW3,EW,Early,2,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,24.92159537,27.84772,26.22237753,1.993,0.00035944,2.452413802,0.000146566,-12.73617
VHDL,VHDL,EW4,EW,Early,2,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,24.71233996,27.84772,26.22237753,1.993,0.00035944,2.833135032,0.00012687,-12.94436
VHDL,VHDL,EW5,EW,Early,2,1,1,Target,124,40,NA,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.76160569,27.84772,26.22237753,1.993,0.00035944,1.374056303,0.00026159,-11.9004
VHDL,VHDL,EW6,EW,Early,2,2,1,Target,39,40,40,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.47244454,27.84772,26.22237753,1.993,0.00036538,3.342856272,0.000109302,-13.1594
VHDL,VHDL,EW7,EW,Early,2,2,1,Target,39,40,40,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.71929254,27.84772,26.22237753,1.993,0.00036538,2.819583273,0.000129586,-12.9138
VHDL,VHDL,EW8,EW,Early,2,2,1,Target,39,40,40,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.24981237,27.84772,26.22237753,1.993,0.00036538,3.897603549,9.37E-05,-13.38091
VHDL,VHDL,MQ3,MQ,Medium,3,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.89835631,27.84772,26.22237753,1.993,0.00035944,1.25039354,0.000287461,-11.76435
VHDL,VHDL,MQ4,MQ,Medium,3,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.52299076,27.84772,26.22237753,1.993,0.00035944,1.619837474,0.000221899,-12.13781
VHDL,VHDL,MQ5,MQ,Medium,3,1,1,Target,124,39.91457748,NA,39.19132996,39.55295372,39.19132996,39.91457748,0.723247528,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,39.57784,25.12014958,27.84772,26.22237753,1.993,0.00048101,2.138579087,0.000224919,-12.11831
VHDL,VHDL,MQ6,MQ,Medium,3,2,1,Target,39,40,40,NA,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.7511079,27.84772,26.22237753,1.993,0.00036538,2.75839205,0.000132461,-12.88215
VHDL,VHDL,MQ7,MQ,Medium,3,2,1,Target,39,40,40,NA,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.76237619,27.84772,26.22237753,1.993,0.00036538,2.737039456,0.000133494,-12.87093
VHDL,VHDL,MQ8,MQ,Medium,3,2,1,Target,39,40,39.12120056,40,39.70706685,39.12120056,40,0.878799438,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,39.70707,24.74641872,27.84772,26.22237753,1.993,0.00044218,2.767326721,0.000159786,-12.61157
VHDL,VHDL,MW3,MW,Medium,4,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.29041589,27.84772,26.22237753,1.993,0.00035944,1.901644416,0.000189015,-12.36921
VHDL,VHDL,MW4,MW,Medium,4,1,1,Target,124,40,40,40,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,24.63093495,27.84772,26.22237753,1.993,0.00035944,2.996737273,0.000119944,-13.02536
VHDL,VHDL,MW5,MW,Medium,4,1,1,Target,124,40,40,NA,40,40,40,0,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,40.02517,25.35247425,27.84772,26.22237753,1.993,0.00035944,1.821974874,0.00019728,-12.30747
VHDL,VHDL,MW6,MW,Medium,4,2,1,Target,39,40,NA,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,24.95817119,27.84772,26.22237753,1.993,0.00036538,2.391327326,0.000152793,-12.67613
VHDL,VHDL,MW7,MW,Medium,4,2,1,Target,39,40,40,40,40,40,40,0,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,40,25.70248086,27.84772,26.22237753,1.993,0.00036538,1.431241263,0.000255288,-11.93558
VHDL,VHDL,MW8,MW,Medium,4,2,1,Target,39,37.46625519,NA,37.58779144,37.52702332,37.46625519,37.58779144,0.121536255,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,37.52702,24.59640071,27.84772,26.22237753,1.993,0.00182903,3.068964925,0.000595978,-10.71245
VHDL,VHDL,LQ3,LQ,Late,5,1,1,Target,124,27.65834236,27.95804024,27.87424278,27.83020846,27.65834236,27.95804024,0.299697876,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,27.84772,26.22237753,27.84772,26.22237753,1.993,1,1,1,0
VHDL,VHDL,LQ4,LQ,Late,5,1,1,Target,124,23.01096535,22.98855782,23.07968521,23.02640279,22.98855782,23.07968521,0.091127396,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,23.04089,25.57874041,27.84772,26.22237753,1.993,22.88782705,1.558741203,14.68353246,3.87613
VHDL,VHDL,LQ5,LQ,Late,5,1,1,Target,124,28.22291756,28.4294548,28.31814003,28.32350413,28.22291756,28.4294548,0.206537247,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,28.34132,26.77102083,27.84772,26.22237753,1.993,0.72507716,0.684979114,1.058539083,0.08207
VHDL,VHDL,LQ6,LQ,Late,5,2,1,Target,39,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,27.59153,26.04640784,27.84772,26.22237753,1.993,1.18157486,1.129026689,1.046542897,0.06563
VHDL,VHDL,LQ7,LQ,Late,5,2,1,Target,39,30.51218987,30.44260597,30.23641014,30.39706866,30.23641014,30.51218987,0.275779724,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,30.39707,24.94259361,27.84772,26.22237753,1.993,0.19007263,2.417155701,0.078634832,-3.66869
VHDL,VHDL,LQ8,LQ,Late,5,2,1,Target,39,23.59091759,23.77414513,23.46027184,23.60844485,23.46027184,23.77414513,0.313873291,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,23.60844,25.43112605,27.84772,26.22237753,1.993,15.81511477,1.725780648,9.164035294,3.19598
VHDL,VHDL,LW3,LW,Late,6,1,1,Target,124,22.55075836,22.6348629,22.55518532,22.58026886,22.55075836,22.6348629,0.084104538,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,22.59447,26.28165682,27.84772,26.22237753,1.993,30.61063312,0.959942941,31.88797146,4.99494
VHDL,VHDL,LW4,LW,Late,6,1,1,Target,124,33.39291382,34.02282333,NA,33.70786858,33.39291382,34.02282333,0.629909515,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,33.72908,25.96422745,27.84772,26.22237753,1.993,0.02170032,1.194862228,0.018161356,-5.78298
VHDL,VHDL,LW5,LW,Late,6,1,1,Target,124,24.07267189,24.0261898,24.16275597,24.08720589,24.0261898,24.16275597,0.136566162,lq6,27.55011559,27.57976532,27.59267616,27.57418569,27.55011559,27.59267616,0.042560577,27.59153366,1.918,24.10236,25.98699537,27.84772,26.22237753,1.993,11.46486958,1.176247401,9.746988235,3.28496
VHDL,VHDL,LW6,LW,Late,6,2,1,Target,39,30.60309029,30.80060768,30.23450089,30.54606628,30.23450089,30.80060768,0.566106796,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,30.54607,25.6146532,27.84772,26.22237753,1.993,0.17249472,1.520610107,0.113437839,-3.14003
VHDL,VHDL,LW7,LW,Late,6,2,1,Target,39,25.11412621,25.11701012,24.99587631,25.07567088,24.99587631,25.11701012,0.121133804,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,25.07567,25.49879946,27.84772,26.22237753,1.993,6.08232351,1.647088425,3.692772903,1.8847
VHDL,VHDL,LW8,LW,Late,6,2,1,Target,39,NA,26.96730614,26.95545006,26.9613781,26.95545006,26.96730614,0.011856079,lq6,27.41981888,27.57100487,27.78377724,27.59153366,27.41981888,27.78377724,0.363958359,27.59153366,1.918,26.96138,25.38152329,27.84772,26.22237753,1.993,1.78114741,1.78583768,0.997373629,-0.00379
Yellow protein,Yellow,EQ3,EQ,Early,1,1,1,Target,47,31.03388405,31.24312592,31.00163078,31.09288025,31.00163078,31.24312592,0.241495132,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,31.17298,25.20347779,33.089,26.22237753,1.993,3.9308,2.019146866,1.946762799,0.96108
Yellow protein,Yellow,EQ4,EQ,Early,1,1,1,Target,47,30.75880432,30.78666687,30.40151405,30.64899508,30.40151405,30.78666687,0.385152817,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,30.72795,25.71323328,33.089,26.22237753,1.993,5.40205836,1.420667419,3.802479235,1.92694
Yellow protein,Yellow,EQ5,EQ,Early,1,1,1,Target,47,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,29.05278,25.41256952,33.089,26.22237753,1.993,17.87766015,1.748007975,10.2274477,3.35437
Yellow protein,Yellow,EQ6,EQ,Early,1,2,1,Target,125,29.94227219,29.72767448,29.5554924,29.74181302,29.5554924,29.94227219,0.386779785,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,29.74181,24.43307804,33.089,26.22237753,1.993,10.92765446,3.434853817,3.181403066,1.66966
Yellow protein,Yellow,EQ7,EQ,Early,1,2,1,Target,125,30.93462372,31.01896286,30.96034813,30.97131157,30.93462372,31.01896286,0.084339142,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,30.97131,25.14049381,33.089,26.22237753,1.993,4.53996441,2.108783874,2.152882744,1.10627
Yellow protein,Yellow,EQ8,EQ,Early,1,2,1,Target,125,29.94722366,29.86710739,30.00184441,29.93872515,29.86710739,30.00184441,0.134737015,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,29.93873,24.92400707,33.089,26.22237753,1.993,9.49360918,2.448338325,3.877572425,1.95515
Yellow protein,Yellow,EW3,EW,Early,2,1,1,Target,47,29.27580261,29.27066422,29.10559654,29.21735446,29.10559654,29.27580261,0.17020607,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,29.29262,24.92159537,33.089,26.22237753,1.993,15.06247651,2.452413802,6.141898441,2.61868
Yellow protein,Yellow,EW4,EW,Early,2,1,1,Target,47,29.83900261,29.8153019,29.64127922,29.76519457,29.64127922,29.83900261,0.197723389,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,29.84187,24.71233996,33.089,26.22237753,1.993,10.17374936,2.833135032,3.590986398,1.84438
Yellow protein,Yellow,EW5,EW,Early,2,1,1,Target,47,29.94815826,29.84382057,29.89508057,29.89568647,29.84382057,29.94815826,0.104337692,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,29.9727,25.76160569,33.089,26.22237753,1.993,9.2659345,1.374056303,6.743489676,2.7535
Yellow protein,Yellow,EW6,EW,Early,2,2,1,Target,125,29.91416359,29.75338364,29.63669968,29.7680823,29.63669968,29.91416359,0.277463913,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,29.76808,24.47244454,33.089,26.22237753,1.993,10.72448456,3.342856272,3.208179978,1.68176
Yellow protein,Yellow,EW7,EW,Early,2,2,1,Target,125,29.97549438,30.13504601,30.20489883,30.10514641,29.97549438,30.20489883,0.229404449,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,30.10515,24.71929254,33.089,26.22237753,1.993,8.42938904,2.819583273,2.989586838,1.57995
Yellow protein,Yellow,EW8,EW,Early,2,2,1,Target,125,30.04839325,30.13326454,29.9674511,30.04970296,29.9674511,30.13326454,0.165813446,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,30.0497,24.24981237,33.089,26.22237753,1.993,8.76997683,3.897603549,2.250094633,1.16999
Yellow protein,Yellow,MQ3,MQ,Medium,3,1,1,Target,47,31.81419754,32.15846634,31.80954361,31.9274025,31.80954361,32.15846634,0.348922729,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,32.00965,25.89835631,33.089,26.22237753,1.993,2.16216114,1.25039354,1.729184512,0.79009
Yellow protein,Yellow,MQ4,MQ,Medium,3,1,1,Target,47,31.44861603,31.36015892,31.46415138,31.42430878,31.36015892,31.46415138,0.103992462,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,31.50526,25.52299076,33.089,26.22237753,1.993,3.100161,1.619837474,1.913871638,0.93649
Yellow protein,Yellow,MQ5,MQ,Medium,3,1,1,Target,47,30.73417854,30.66982841,30.75903511,30.72101402,30.66982841,30.75903511,0.089206696,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,30.80016,25.12014958,33.089,26.22237753,1.993,5.13046238,2.138579087,2.399005215,1.26244
Yellow protein,Yellow,MQ6,MQ,Medium,3,2,1,Target,125,31.23072243,31.59446335,31.60165787,31.47561455,31.23072243,31.60165787,0.37093544,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,31.47561,24.7511079,33.089,26.22237753,1.993,3.16652694,2.75839205,1.147961161,0.19907
Yellow protein,Yellow,MQ7,MQ,Medium,3,2,1,Target,125,30.97146034,30.96819115,31.2552166,31.06495603,30.96819115,31.2552166,0.287025452,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,31.06496,24.76237619,33.089,26.22237753,1.993,4.24617154,2.737039456,1.551373886,0.63355
Yellow protein,Yellow,MQ8,MQ,Medium,3,2,1,Target,125,31.59382057,31.22221565,31.48219872,31.43274498,31.22221565,31.59382057,0.371604919,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,31.43274,24.74641872,33.089,26.22237753,1.993,3.26500806,2.767326721,1.179841916,0.23859
Yellow protein,Yellow,MW3,MW,Medium,4,1,1,Target,47,29.48028946,28.96175003,29.49831772,29.3134524,28.96175003,29.49831772,0.536567688,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,29.38897,25.29041589,33.089,26.22237753,1.993,14.06058592,1.901644416,7.39390908,2.88634
Yellow protein,Yellow,MW4,MW,Medium,4,1,1,Target,47,29.92942238,29.66550827,29.67867279,29.75786781,29.66550827,29.92942238,0.263914108,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,29.83453,24.63093495,33.089,26.22237753,1.993,10.22728016,2.996737273,3.412805072,1.77096
Yellow protein,Yellow,MW5,MW,Medium,4,1,1,Target,47,31.00597,31.02233315,31.37689018,31.13506444,31.00597,31.37689018,0.370920181,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,31.21527,25.35247425,33.089,26.22237753,1.993,3.81380786,1.821974874,2.09322747,1.06573
Yellow protein,Yellow,MW6,MW,Medium,4,2,1,Target,125,29.84990692,30.01970673,29.51557922,29.79506429,29.51557922,30.01970673,0.504127502,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,29.79506,24.95817119,33.089,26.22237753,1.993,10.51973429,2.391327326,4.399119342,2.13721
Yellow protein,Yellow,MW7,MW,Medium,4,2,1,Target,125,32.70870972,32.96828842,32.66572189,32.78090668,32.66572189,32.96828842,0.302566528,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,32.78091,25.70248086,33.089,26.22237753,1.993,1.24621315,1.431241263,0.870721928,-0.19972
Yellow protein,Yellow,MW8,MW,Medium,4,2,1,Target,125,30.65320778,30.28611374,30.36319923,30.43417358,30.28611374,30.65320778,0.36709404,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,30.43417,24.59640071,33.089,26.22237753,1.993,6.66360697,3.068964925,2.171288084,1.11855
Yellow protein,Yellow,LQ3,LQ,Late,5,1,1,Target,47,32.68360138,32.87354279,33.45479202,33.00397873,32.68360138,33.45479202,0.771190643,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,33.089,26.22237753,33.089,26.22237753,1.993,1,1,1,0
Yellow protein,Yellow,LQ4,LQ,Late,5,1,1,Target,47,35.13106537,34.36470413,34.77955627,34.75844193,34.36470413,35.13106537,0.766361237,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,34.84799,25.57874041,33.089,26.22237753,1.993,0.28460546,1.558741203,0.182586729,-2.45335
Yellow protein,Yellow,LQ5,LQ,Late,5,1,1,Target,47,35.12080383,35.3639946,NA,35.24239922,35.12080383,35.3639946,0.243190765,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,35.33319,26.77102083,33.089,26.22237753,1.993,0.20123326,0.684979114,0.293780138,-1.76719
Yellow protein,Yellow,LQ6,LQ,Late,5,2,1,Target,125,33.98675156,34.49580383,34.21722031,34.23325857,33.98675156,34.49580383,0.509052277,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,34.23326,26.04640784,33.089,26.22237753,1.993,0.44154349,1.129026689,0.391083302,-1.35445
Yellow protein,Yellow,LQ7,LQ,Late,5,2,1,Target,125,30.31091118,30.14108849,29.70967865,30.05389277,29.70967865,30.31091118,0.601232529,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,30.05389,24.94259361,33.089,26.22237753,1.993,8.74376507,2.417155701,3.617377677,1.85494
Yellow protein,Yellow,LQ8,LQ,Late,5,2,1,Target,125,33.19939041,33.42623901,33.1549263,33.26018524,33.1549263,33.42623901,0.271312714,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,33.26019,25.43112605,33.089,26.22237753,1.993,0.88488606,1.725780648,0.512745383,-0.96369
Yellow protein,Yellow,LW3,LW,Late,6,1,1,Target,47,34.37603378,NA,33.93990707,34.15797043,33.93990707,34.37603378,0.436126709,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,34.24597,26.28165682,33.089,26.22237753,1.993,0.43755283,0.959942941,0.455811285,-1.13349
Yellow protein,Yellow,LW4,LW,Late,6,1,1,Target,47,NA,34.93208313,34.22764969,34.57986641,34.22764969,34.93208313,0.704433441,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,34.66895,25.96422745,33.089,26.22237753,1.993,0.32343893,1.194862228,0.270691402,-1.88528
Yellow protein,Yellow,LW5,LW,Late,6,1,1,Target,47,33.00194931,33.42869568,32.94427109,33.12497203,32.94427109,33.42869568,0.484424591,eq5,29.06079483,28.97247124,28.90113068,28.97813225,28.90113068,29.06079483,0.159664154,29.0527846,2.043,33.21031,25.98699537,33.089,26.22237753,1.993,0.91698641,1.176247401,0.77958634,-0.35922
Yellow protein,Yellow,LW6,LW,Late,6,2,1,Target,125,33.42180252,33.52858734,33.57720566,33.50919851,33.42180252,33.57720566,0.155403137,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,33.5092,25.6146532,33.089,26.22237753,1.993,0.74067283,1.520610107,0.487089246,-1.03774
Yellow protein,Yellow,LW7,LW,Late,6,2,1,Target,125,34.96672058,35.27972031,35.67523193,35.30722427,34.96672058,35.67523193,0.708511353,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,35.30722,25.49879946,33.089,26.22237753,1.993,0.20500096,1.647088425,0.124462632,-3.00622
Yellow protein,Yellow,LW8,LW,Late,6,2,1,Target,125,34.46694183,34.92765427,34.45300293,34.61586634,34.45300293,34.92765427,0.474651337,eq5,29.01953125,29.01187515,29.1269474,29.0527846,29.01187515,29.1269474,0.11507225,29.0527846,2.043,34.61587,25.38152329,33.089,26.22237753,1.993,0.33594052,1.78583768,0.188113693,-2.41032")


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

##################################################################################################################
##################################################################################################################
###########################################Apis comparison analysis ##############################################

########### Cameron et al comparison##########

# This script is a modification of David Prince's script to compare the mRNAseq results from Bombus terrestris with a microarray
# on Apis mellifera published by Cameron et al. 2013
#-----------------------------------------------------------------------------------

# DESCRIPTION

# Purpose: to calculate whether the overlap between two gene lists is statistically significant

## Situation - the two lists to be compared are from different species
## AND
## one of the lists is a microarray study
### Necessary data conversion - convert the gene list from one of the species into a list of orthologs in the other species
### Necessary inputs: 
#### input1 = a list of genes from this study as .txt file - For example "DC_EQ_AllSig.txt" (10 files)
#### input2 = a list of genes from the study to be compared to as a .txt file - "Cameron_108_worker" and "Cameron_004_queen" (14 files)
#### input3 = the number of genes with orthologs in both species - "new_conversion_sheet_both_part2.csv"
#### input4 = a list of the total number of genes that appear in the study to be compared to as a .txt file -
####                                                        "Cameron_et_al_2013_microarray.csv"

## Output will be:
## the number of genes that overlap between the two lists
## a p-value for how significant the overlap between the two lists is
## an odds ratio (representing the strength and direction of the association)

# The dataset is from a microarray, therefore not all genes may be present on the array
# Load the list of genes present on the microarray. 
#The original authors and full citation is Cameron, RC, Duncan, EL, and Dearden, PK. 2013. Biased gene expression in early honeybee larval development. BMC Genomics 14 (903).
## File taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52288 #The authors of the present study were not involved in producing these data, and do not claim credit for them. They are reproduced here for comparison purposes only.
microarray = read.csv("Cameron_et_al_2013_microarray.csv", header = TRUE) 

microarray <- read_csv("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/RAnalysis/comparative_analyses/BombusApisCaste/Cameron_et_al_2013_microarray.csv")

# Isolate the Reporter.Comment column to find the names of the genes present on the microarray
gene.names = as.character(microarray$Beebase)

# Some genes on the microarray are represented by more than one spot
# Therefore identify the unique gene Ids (removing redundancy)
unique.gene.names = unique(gene.names)

# Calculate which of the unique.gene.names have an RBH BLAST ortholog in Bombus terrestris
# Input Bombus terrestris to Apis mellifera RBH BLAST conversion list
conversion.df = read.csv("new_conversion_sheet_both_part2.csv", header = TRUE) 

# The Cameron et al 2013 data set uses OGS_v1_1 Ids, therefore remove the other Apis Ids
conversion.df2 = conversion.df[,c(4,6)]
# and the instances where no RBH BLAST ortholog exists between OGS_v1_1 and Bter Gene_Ids
conversion.df3 = na.omit(conversion.df2)
conversion.df3 = unique(conversion.df3) #removes redundancy
# Calculate the overlap between the RBH BLAST ortholog list (the column OGS_v1_1 in the 'new conversion sheet 2') 
#and the microarray gene list (the values labelled 'unique.gene.names')
microarray.RBH = intersect(conversion.df3$OGS_v1_1, unique.gene.names)
# count the number of ortholog pairs between Bombus terrestris and Apis mellifera on the microarray
total.gene.list.size = length(microarray.RBH)
#Note the total number of orthologs is 7598

## Define a function for comparing the gene lists

compare.bt.amel.lists = function (x, y)
{
  # Import the necessary files
  bombus.list = readLines(x)
  apis.list = readLines(y)
  
  # Calculate the Bombus terrestris genes with Apis mellifera orthologs
  bombus.present.in.df = conversion.df3[match(bombus.list, conversion.df3$Bombus_Gene_ID),]
  bombus.present.in.df2 = na.omit(bombus.present.in.df)
  count.removed.NAs.bombus = length(bombus.present.in.df$Bombus_Gene_ID) - length(bombus.present.in.df2$Bombus_Gene_ID)
  print("Number of Bombus terrestris genes in the list with no Apis mellifera homolog:")
  print(count.removed.NAs.bombus)
  
  # Place the Bombus terrestris genes with orthologs in Apis mellifera into a list
  bombus.ortho.list = as.character(bombus.present.in.df2$Bombus_Gene_ID)
  
  # Calculate the Apis mellifera genes with Bombus terrestis orthologs
  apis.present.in.df = conversion.df3[match(apis.list, conversion.df3$OGS_v1_1),]
  apis.present.in.df2 = na.omit(apis.present.in.df)
  count.removed.NAs.apis = length(apis.present.in.df$OGS_v1_1) - length(apis.present.in.df2$OGS_v1_1)
  print("Number of Apis mellifera genes in the list with no Bombus terrestris homolog:")
  print(count.removed.NAs.apis)
  
  # Place the Apis mellifera genes with orthologs in Bombus terrestris into a list
  apis.ortho.list = as.character(apis.present.in.df2$Bombus_Gene_ID)
  
  # Determine the genes present in both list1 and list2, and store the number of genes separately
  both.lists = intersect(bombus.ortho.list, apis.ortho.list)
  print("Overlapping genes are:")
  print(both.lists)
  number.both.lists = length(both.lists)
  
  # Determine the genes present in list1 but not list2, and store the number of genes separately
  only.in.bombus = setdiff(bombus.ortho.list, apis.ortho.list)
  number.only.in.bombus = length(only.in.bombus)
  
  # Determine the genes present in list2 but not list1, and store the number of genes separately
  only.in.apis = setdiff(apis.ortho.list, bombus.ortho.list)
  number.only.in.apis = length(only.in.apis)
  
  # Determine the number of genes not present in either list
  neither.list = total.gene.list.size - (number.both.lists + number.only.in.bombus + number.only.in.apis)
  
  # Create a matrix representing the numbers of genes in both lists, only list1, only list2, and neither list
  contingency.table = matrix(c(number.both.lists, number.only.in.bombus, number.only.in.apis, neither.list))
  # Change the dimensions to 2 rows and 2 columns
  dim(contingency.table) <- c(2,2)
  
  # Print the matrix 
  print("Contingency table")
  print(contingency.table)
  
  # Conduct two-tailed Fisher's exact test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance
  print(fisher.test(contingency.table))
}

# Run comparisons between gene lists with the compare.bt.amel.lists function
# Form is comparelists("bombus list", "apis list")
# Gene list name needs to be a string


####Early instar significant (LFC>0) comparisons

#Queen
EQAllsig_006_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_006_queen.txt")
EQAllsig_012_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_012_queen.txt")
EQAllsig_036_queen = compare.bt.amel.lists("DC_EQ_AllSig.txt","Cameron_036_queen.txt")

#Worker
EWAllsig_006_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_006_worker.txt")
EWAllsig_012_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_012_worker.txt")
EWAllsig_036_worker = compare.bt.amel.lists("DC_EW_AllSig.txt","Cameron_036_worker.txt")


####Mid instar significant (LFC>0) comparisons

#Queen
MQAllsig_036_queen = compare.bt.amel.lists("DC_MQ_AllSig.txt","Cameron_036_queen.txt")
MQAllsig_060_queen = compare.bt.amel.lists("DC_MQ_AllSig.txt","Cameron_060_queen.txt")
MQAllsig_084_queen = compare.bt.amel.lists("DC_MQ_AllSig.txt","Cameron_084_queen.txt")

#Worker
MWAllsig_036_worker = compare.bt.amel.lists("DC_MW_AllSig.txt","Cameron_036_worker.txt")
MWAllsig_060_worker = compare.bt.amel.lists("DC_MW_AllSig.txt","Cameron_060_worker.txt")
MWAllsig_084_worker = compare.bt.amel.lists("DC_MW_AllSig.txt","Cameron_084_worker.txt")


####Late instar significant (LFC>0) comparisons

#Queen
LQAllsig_084_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_084_queen.txt")
LQAllsig_108_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_108_queen.txt")
LQAllsig_132_queen = compare.bt.amel.lists("DC_LQ_AllSig.txt","Cameron_132_queen.txt")

#Worker
LWAllsig_084_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_084_worker.txt")
LWAllsig_108_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_108_worker.txt")
LWAllsig_132_worker = compare.bt.amel.lists("DC_LW_AllSig.txt","Cameron_132_worker.txt")



####Mid instar highly significant (LFC>1) comparisons

#Queen
MQHighsig_036_queen = compare.bt.amel.lists("DC_MQ_HighSig.txt","Cameron_036_queen.txt")
MQHighsig_060_queen = compare.bt.amel.lists("DC_MQ_HighSig.txt","Cameron_060_queen.txt")
MQHighsig_084_queen = compare.bt.amel.lists("DC_MQ_HighSig.txt","Cameron_084_queen.txt")

#Worker
MWHighsig_036_worker = compare.bt.amel.lists("DC_MW_HighSig.txt","Cameron_036_worker.txt")
MWHighsig_060_worker = compare.bt.amel.lists("DC_MW_HighSig.txt","Cameron_060_worker.txt")
MWHighsig_084_worker = compare.bt.amel.lists("DC_MW_HighSig.txt","Cameron_084_worker.txt")


####Late instar highly significant (LFC>1) comparisons

#Queen
LQHighsig_084_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_084_queen.txt")
LQHighsig_108_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_108_queen.txt")
LQHighsig_132_queen = compare.bt.amel.lists("DC_LQ_HighSig.txt","Cameron_132_queen.txt")

#Worker
LWHighsig_084_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_084_worker.txt")
LWHighsig_108_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_108_worker.txt")
LWHighsig_132_worker = compare.bt.amel.lists("DC_LW_HighSig.txt","Cameron_132_worker.txt")



########### He et al comparison##########

# This script is a modification of David Prince's script to compare the mRNAseq results from He et al 2017 with the results
# from this study.
#-----------------------------------------------------------------------------------

# DESCRIPTION

# Purpose: to calculate whether the overlap between two gene lists is statistically significant

## Situation - the two lists to be compared are from different species
### Necessary data conversion - convert the gene list from one of the species into a list of orthologs in the other species
### Necessary inputs: 
#### input1 = a list of genes from this study as .txt file - Files labelled for example "DC_EQ_AllSig.txt" (10 files)
#### input2 = a list of genes from the study to be compared to as a .txt file - Files labelled for example "He_2d_queen","He_2d_worker" (4 files)
#### input3 = the number of genes with orthologs in both species - "new_conversion_sheet_both_part2.csv"
#### input4 = a list of the total number of genes that appear in the study to be compared to as a .txt file -
####                                                        "He_et_al_2017_All_genes_count.csv"

## Output will be:
## the number of genes that overlap between the two lists
## a p-value for how significant the overlap between the two lists is
## an odds ratio (representing the strength and direction of the association)

# The dataset is from mRNAseq, but includes several unique identifiers that are not comparable between studies (for their novel genes)
# We will only look at the genes that are present in both studies.
# I recieved the total number of genes sequenced by mRNAseq from the authors themselves. They come from He XJ, Jiang WJ, Zhou M, Barron AB, Zeng ZJ. 2017. A comparison of honeybee (Apis mellifera) queen, worker and drone larvae by RNA-Seq. Insect Sci.26 (3) 499
mRNAseq = read.csv("He_et_al_2017_All_genes_count.csv", header = TRUE) 

# Isolate the BeeBase column to find the names of the genes present in the mRNAseq
He.gene.names = as.character(mRNAseq[,1])

# Some genes might be represented twice so remove these
# Therefore identify the unique gene Ids (removing redundancy)
He.unique.gene.names = unique(He.gene.names)

# Calculate which of the unique.gene.names have an RBH BLAST ortholog in Bombus terrestris
# Input Bombus terrestris to Apis mellifera RBH BLAST conversion list
He.conversion.df = read.csv("new_conversion_sheet_both_part2.csv", header = TRUE) 

# The He et al data is annotated using the OGS_v2 on beebase (i.e. after 2014)
#therefore match to beebase ids, to do this remove the other Apis Ids.
He.conversion.df2 = He.conversion.df[,c(3,6)]
# and the instances where no RBH BLAST ortholog exists between OGS_v1_1 and Bter Gene_Ids
He.conversion.df3 = na.omit(He.conversion.df2)
He.conversion.df3 = unique(He.conversion.df3) #removes redundancy
# Calculate the overlap between the RBH BLAST ortholog list (the column BEEBASE in the 'new conversion sheet 2') 
#and the gene list (the values labelled 'unique.gene.names')
mRNAseq.RBH = intersect(He.conversion.df3$BEEBASE, He.unique.gene.names)
# count the number of ortholog pairs between Bombus terrestris and Apis mellifera
He.total.gene.list.size = length(mRNAseq.RBH)
#Note the total number of orthologs is 6677

## Define a function for comparing the gene lists

compare.bt.He.amel.lists = function (x, y)
{
  # Import the necessary files
  bombus.list = readLines(x)
  apis.list = readLines(y)
  
  # Calculate the Bombus terrestris genes with Apis mellifera orthologs
  bombus.present.in.df = He.conversion.df3[match(bombus.list, He.conversion.df3$Bombus_Gene_ID),]
  bombus.present.in.df2 = na.omit(bombus.present.in.df)
  count.removed.NAs.bombus = length(bombus.present.in.df$Bombus_Gene_ID) - length(bombus.present.in.df2$Bombus_Gene_ID)
  print("Number of Bombus terrestris genes in the list with no Apis mellifera homolog:")
  print(count.removed.NAs.bombus)
  
  # Place the Bombus terrestris genes with orthologs in Apis mellifera into a list
  bombus.ortho.list = as.character(bombus.present.in.df2$Bombus_Gene_ID)
  
  # Calculate the Apis mellifera genes with Bombus terrestis orthologs
  apis.present.in.df = He.conversion.df3[match(apis.list, He.conversion.df3$BEEBASE),]
  apis.present.in.df2 = na.omit(apis.present.in.df)
  count.removed.NAs.apis = length(apis.present.in.df$BEEBASE) - length(apis.present.in.df2$BEEBASE)
  print("Number of Apis mellifera genes in the list with no Bombus terrestris homolog:")
  print(count.removed.NAs.apis)
  
  # Place the Apis mellifera genes with orthologs in Bombus terrestris into a list
  apis.ortho.list = as.character(apis.present.in.df2$Bombus_Gene_ID)
  
  # Determine the genes present in both list1 and list2, and store the number of genes separately
  both.lists = intersect(bombus.ortho.list, apis.ortho.list)
  print("Overlapping genes are:")
  print(both.lists)
  number.both.lists = length(both.lists)
  
  # Determine the genes present in list1 but not list2, and store the number of genes separately
  only.in.bombus = setdiff(bombus.ortho.list, apis.ortho.list)
  number.only.in.bombus = length(only.in.bombus)
  
  # Determine the genes present in list2 but not list1, and store the number of genes separately
  only.in.apis = setdiff(apis.ortho.list, bombus.ortho.list)
  number.only.in.apis = length(only.in.apis)
  
  # Determine the number of genes not present in either list
  neither.list = He.total.gene.list.size - (number.both.lists + number.only.in.bombus + number.only.in.apis)
  
  # Create a matrix representing the numbers of genes in both lists, only list1, only list2, and neither list
  contingency.table = matrix(c(number.both.lists, number.only.in.bombus, number.only.in.apis, neither.list))
  # Change the dimensions to 2 rows and 2 columns
  dim(contingency.table) <- c(2,2)
  
  # Print the matrix 
  print("Contingency table")
  print(contingency.table)
  
  # Conduct two-tailed Fisher's exact test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance
  print(fisher.test(contingency.table))
}

# Run comparisons between gene lists with the compare.bt.amel.lists function
# Form is comparelists("bombus list", "apis list")
# Bombus lists are available in "Supplementary datafile 4", Apis lists are located in the supplementary material for 
# He et al. 2017. Insect Sci. 00, 1-11, each gene list has been split into separate queen and worker lists at each time point (2d and 4d), and then stored
#as separate text files
## Gene list name needs to be a string

#### Early instar significant comparisons

#Within caste
EWsig_2dw = compare.bt.He.amel.lists("DC_EW_AllSig.txt","He_2d_worker.txt")
EQsig_2dq = compare.bt.He.amel.lists("DC_EQ_AllSig.txt","He_2d_queen.txt")

#Between caste
EWsig_2dq = compare.bt.He.amel.lists("DC_EW_AllSig.txt","He_2d_queen.txt")
EQsig_2dw = compare.bt.He.amel.lists("DC_EQ_AllSig.txt","He_2d_worker.txt")



##### Mid instar significant comparisons

#Within caste
MWsig_2dw = compare.bt.He.amel.lists("DC_MW_AllSig.txt","He_2d_worker.txt")
MQsig_2dq = compare.bt.He.amel.lists("DC_MQ_AllSig.txt","He_2d_queen.txt")
MWsig_4dw = compare.bt.He.amel.lists("DC_MW_AllSig.txt","He_4d_worker.txt")
MQsig_4dq = compare.bt.He.amel.lists("DC_MQ_AllSig.txt","He_4d_queen.txt")

#Between caste
MWsig_2dq = compare.bt.He.amel.lists("DC_MW_AllSig.txt","He_2d_queen.txt")
MQsig_2dw = compare.bt.He.amel.lists("DC_MQ_AllSig.txt","He_2d_worker.txt")
MWsig_4dq = compare.bt.He.amel.lists("DC_MW_AllSig.txt","He_4d_queen.txt")
MQsig_4dw = compare.bt.He.amel.lists("DC_MQ_AllSig.txt","He_4d_worker.txt")




#### Late instar significant comparisons

#Within caste
LWsig_2dw = compare.bt.He.amel.lists("DC_LW_AllSig.txt","He_2d_worker.txt")
LQsig_2dq = compare.bt.He.amel.lists("DC_LQ_AllSig.txt","He_2d_queen.txt")
LWsig_4dw = compare.bt.He.amel.lists("DC_LW_AllSig.txt","He_4d_worker.txt")
LQsig_4dq = compare.bt.He.amel.lists("DC_LQ_AllSig.txt","He_4d_queen.txt")

#Between caste
LWsig_2dq = compare.bt.He.amel.lists("DC_LW_AllSig.txt","He_2d_queen.txt")
LQsig_2dw = compare.bt.He.amel.lists("DC_LQ_AllSig.txt","He_2d_worker.txt")
LWsig_4dq = compare.bt.He.amel.lists("DC_LW_AllSig.txt","He_4d_queen.txt")
LQsig_4dw = compare.bt.He.amel.lists("DC_LQ_AllSig.txt","He_4d_worker.txt")



#### Medium instar highly significant (LFC>1) comparisons 

#Within caste
MWHighsig_2dw = compare.bt.He.amel.lists("DC_MW_HighSig.txt","He_2d_worker.txt")
MQHighsig_2dq = compare.bt.He.amel.lists("DC_MQ_HighSig.txt","He_2d_queen.txt")
MWHighsig_4dw = compare.bt.He.amel.lists("DC_MW_HighSig.txt","He_4d_worker.txt")
MQHighsig_4dq = compare.bt.He.amel.lists("DC_MQ_HighSig.txt","He_4d_queen.txt")

#Between caste
MWHighsig_2dq = compare.bt.He.amel.lists("DC_MW_HighSig.txt","He_2d_queen.txt")
MQHighsig_2dw = compare.bt.He.amel.lists("DC_MQ_HighSig.txt","He_2d_worker.txt")
MWHighsig_4dq = compare.bt.He.amel.lists("DC_MW_HighSig.txt","He_4d_queen.txt")
MQHighsig_4dw = compare.bt.He.amel.lists("DC_MQ_HighSig.txt","He_4d_worker.txt")


####Late instar highly significant (LFC>1) comparisons

#Within caste
LWHighsig_2dw = compare.bt.He.amel.lists("DC_LW_HighSig.txt","He_2d_worker.txt")
LQHighsig_2dq = compare.bt.He.amel.lists("DC_LQ_HighSig.txt","He_2d_queen.txt")
LWHighsig_4dw = compare.bt.He.amel.lists("DC_LW_HighSig.txt","He_4d_worker.txt")
LQHighsig_4dq = compare.bt.He.amel.lists("DC_LQ_HighSig.txt","He_4d_queen.txt")

#Between caste
LWHighsig_2dq = compare.bt.He.amel.lists("DC_LW_HighSig.txt","He_2d_queen.txt")
LQHighsig_2dw = compare.bt.He.amel.lists("DC_LQ_HighSig.txt","He_2d_worker.txt")
LWHighsig_4dq = compare.bt.He.amel.lists("DC_LW_HighSig.txt","He_4d_queen.txt")
LQHighsig_4dw = compare.bt.He.amel.lists("DC_LQ_HighSig.txt","He_4d_worker.txt")


##################################################################################################################
##################################################################################################################
###########################################Euler diagrams for figure 5############################################

#Below are the Euler diagrams for figure 5 in the paper

#Euler diagram for comparison with Cameron 60 hr old queen-destined larvae with MQ  
Cameron.60hr.vs.MQ <- euler(c(Bombus = 2346, 
                              Apis = 204, 
                              "Apis&Bombus" = 51))
ePlot.Cameron.60hr.vs.MQ <- as.ggplot(plot(Cameron.60hr.vs.MQ, quantities = TRUE, main = "A"))

#Euler diagram for comparison with Cameron 108 hr old queen-destined larvae with LQ  
Cameron.108hr.vs.LQ <- euler(c(Bombus = 1166, 
                               Apis = 210, 
                               "Apis&Bombus" = 8))
ePlot.Cameron.108hr.vs.LQ <- as.ggplot(plot(Cameron.108hr.vs.LQ, quantities = TRUE, main = "B"))

#Euler diagram for comparison with He 4 day old queen-destined larvae with MQ  
He.4d.Queen.vs.MQ <- euler(c(Bombus = 2506, 
                             Apis = 100, 
                             "Apis&Bombus" = 33))
ePlot.He.4d.Queen.vs.MQ <- as.ggplot(plot(He.4d.Queen.vs.MQ, quantities = TRUE, main = "C"))

#Euler diagram for comparison with He 4 day old worker-destined larvae with MW  
He.4d.Worker.vs.MW <- euler(c(Bombus = 830, 
                              Apis = 129, 
                              "Apis&Bombus" = 37))
eplot.He.4d.Worker.vs.MW <- as.ggplot(plot(He.4d.Worker.vs.MW, quantities = TRUE, main = "D"))

#Euler diagram for comparison with He 2 day old queen-destined larvae with MQ using just the HDEGs
He.2d.Queen.vs.MQ.HDEG <- euler(c(Bombus = 18, 
                                  Apis = 40, 
                                  "Apis&Bombus" = 3))
eplot.He.2d.Queen.vs.MQ.HDEG <- as.ggplot(plot(He.2d.Queen.vs.MQ.HDEG, quantities = TRUE, main = "E"))

#Euler diagram for comparison with He 2 day old worker-destined larvae with MW using just the HDEGs
He.2d.Worker.vs.MW.HDEG <- euler(c(Bombus = 14, 
                                   Apis = 162, 
                                   "Apis&Bombus" = 4))
ePlot.He.2d.Worker.vs.MW.HDEG <- as.ggplot(plot(He.2d.Worker.vs.MW.HDEG, quantities = TRUE, main = "F"))


ePlot.Cameron.60hr.vs.MQ
ePlot.Cameron.108hr.vs.LQ
ePlot.He.4d.Queen.vs.MQ
eplot.He.4d.Worker.vs.MW
eplot.He.2d.Queen.vs.MQ.HDEG
ePlot.He.2d.Worker.vs.MW.HDEG

Euler.panel <- gridExtra::arrangeGrob(ePlot.Cameron.60hr.vs.MQ,
                                      ePlot.Cameron.108hr.vs.LQ,
                                      ePlot.He.4d.Queen.vs.MQ,
                                      eplot.He.4d.Worker.vs.MW,
                                      eplot.He.2d.Queen.vs.MQ.HDEG,
                                      ePlot.He.2d.Worker.vs.MW.HDEG,
                                      ncol = 3) #
grid::grid.newpage() #clear current graphic and start a new page (not sure why this is necessary)
grid::grid.draw(Euler.panel) #draw your plot on the newpage



