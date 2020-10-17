#DavidHCollins
#6 08 2020

#This script is relevant to the following publication:

#David H. Collins1*, Anders Wirén1, 2*, Marjorie Labédan1, 3, Michael Smith1, David C. Prince1, Irina Mohorianu1, 4, Tamas Dalmay1, and Andrew F. G. Bourke1. Gene expression in larval caste determination of bumblebees and the transition to advanced eusociality.

#This script is to use the gene expression data to produce a Heatmap that summarises the  results across all genes.

#File required: 'all_genes.csv' - this file is a single tab csv file of Table S9 containing all of the information on that table, an extra column that shows the phenotype that each gene is expressed in 


###############################################################################################

#First you must alwys reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(readr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)
library(ExPosition)



#Set working directory
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Academic writing/My papers/Collins et al_BB-M001482-1_obj1 paper/R scripts/")

#view 'all_genes' file
all_genes <- read_csv("all_genes.csv")
glimpse(all_genes)

#Select the columns that will be most relevant for this analysis
df <- all_genes %>% 
  select(phenotype,feature_ID,"EQ1":"LW3",DEG,LOFC,HDEG,type,gene_symbol,accession_number,product_ID)

#Name the unannotated genes Unannotated (genome location)
df <- df %>%  
  mutate(name=ifelse(is.na (product_ID), "unannotated",product_ID)) %>% 
  mutate(brackets = paste0(" (", feature_ID, ")")) %>% 
  mutate(newname = ifelse(name== "unannotated", 
                          paste0(name, brackets),
                          paste0(name))) %>% 
  select(-(name),-(brackets))


df2 <- df %>% 
  arrange(desc(LOFC)) %>% 
  distinct(feature_ID, .keep_all = TRUE) %>% 
    slice_max(order_by = LOFC, n = 50) %>% 
  select(newname, feature_ID,accession_number, "EQ1":"LW3")   

#Find all the rows which have the same product name and rename them so they don't confuse the figure
 df2 %>%
  group_by(newname) %>% 
  filter(n()>1) %>% 
   select(newname,feature_ID,accession_number)
 
 
 df2 <- df2 %>% 
   mutate(newname=replace(newname, feature_ID=="NC_015771.1_3137226_3141484","plrp2-like (XM_003398387.2)"),
          newname=replace(newname, feature_ID=="NC_015772.1_3772353_3775563","plrp2-like (XM_012313415.1)"),
          newname=replace(newname, feature_ID=="NC_015764.1_10241255_10243195","Nrf-6 (XM_012321236.1)"),
          newname=replace(newname, feature_ID=="NC_015764.1_10150985_10155037","Nrf-6 (XM_003394457.2)"))
 
#Find all uncharacterised genes and make them consistent in their names with other genes
 df2 %>%
   filter(grepl("uncharacterized",newname)) %>% 
   select(newname,feature_ID,accession_number)

 df2<- df2 %>% 
   mutate(newname=replace(newname, newname=="uncharacterized LOC100645614","uncharacterized (XM_012319905.1)"),
          newname=replace(newname, newname=="putative uncharacterized protein DDB_G0271606","uncharacterized (XM_003399878.3)"),
          newname=replace(newname, newname=="uncharacterized LOC100646546","uncharacterized (XM_003397986.2)"),
          newname=replace(newname, newname=="uncharacterized LOC100645518","uncharacterized (XM_012321073.1)"),
          newname=replace(newname, newname=="uncharacterized LOC100642732","uncharacterized (XM_012312992.1)"))
 
#Find all names, and make the overly long ones shorter
df2$newname

df2 <- df2 %>% 
  mutate(newname= replace(newname, newname== "retinal guanylyl cyclase 2-like%2C transcript variant X2", "retinal guanylyl cyclase 2-like"),
         newname=replace(newname, newname=="esterase E4-like%2C transcript variant X1", "esterase E4-like"),
         newname=replace(newname, newname=="short-chain dehydrogenase/reductase family 16C member 6%2C transcript variant X1", "short-chain dehydrogenase/reductase family 16C"),
         newname=replace(newname, newname=="probable phospholipase A1 magnifin%2C transcript variant X1", "probable phospholipase A1 magnifin"),
         newname=replace(newname, newname=="phospholipase B1%2C membrane-associated-like", "phospholipase B12C membrane-associated-like"),
         newname=replace(newname, newname=="26S proteasome non-ATPase regulatory subunit 1-like", "26S proteasome non-ATPase reg.s1-like"),
         newname=replace(newname, newname=="glucose dehydrogenase [FAD%2C quinone]-like", "glucose dehydrogenase [FAD quinone]-like"),
         newname=replace(newname, newname=="short-chain dehydrogenase/reductase family 16C", "short-chain dehydrogenase 16C"),
         newname=replace(newname, newname=="endocuticle structural glycoprotein SgAbd-8-like", "endocuticle glycoprotein SgAbd-8-like")) %>% 
  select(newname,"EQ1":"LW3") 

#Make the df a matrix
mx2 <- as.matrix(df2[,-1])
rownames(mx2) <- df2$newname

#Create names for annotation rows (if developmental stages wanted)
#my_sample_col <- data.frame(Pathway = rep(c("Q", "W","Q","W","Q","W"), c(3,3,3,2,3,3)),
#                            Developmental_stage = rep(c("Early", "Mid", "Late"), c(6,5,6))
#)

#Specify desired colours for annotation columns (if development stage wanted)
#annotation_colors = list(
#  Pathway = c(Q="black", W="white"),
#  Developmental_stage = c(Early="azure4", Mid="azure4", Late="azure4"))

#Otherwise just use castes for annotation column names
my_sample_col2 <- data.frame(Pathway = rep(c("Q", "W","Q","W","Q","W"), c(3,3,3,2,3,3)))
row.names(my_sample_col2) <- colnames(mx2)

#And annotation colours
annotation_colors2 = list(
  Pathway = c(Q="black", W="white"))

#This code scales the gene expression values in a custom way - but makes little difference to the overall results
#x <- rowNorms(mx2, type="z", scale = T)

#top50 <- pheatmap(x,
#                  cluster_rows=TRUE,
#                 cluster_cols=FALSE,
#                 cutree_rows=4,
#                  gaps_col = c(6,11),
#                  annotation_col = my_sample_col2,
#                  annotation_colors = annotation_colors,
#                  annotation_legend=F,
#                  legend = T)

#Make heatmap and save it
(top50 <- pheatmap(mx2,
                  scale="row",
                  cluster_rows=T,
                  cluster_cols=F,
                  annotation_names_col=F,
                  cutree_rows=4,
                  gaps_col = c(6,11),
                  annotation_col = my_sample_col2,
                  annotation_colors = annotation_colors2,
                  annotation_legend=T,
                  legend = T,
                  cellwidth=15,
                  cellheight=10,
                  fontsize=8,
                  filename = "HDEG_heatmap1.pdf"))

ggsave("HDEG_heatmap1.pdf", top50, width = 10, height = 16)

### Other options for displaying Heatmaps are hashed in below

### Option 2 : top 10 in each phenotype across all phenotypes

#df3 <- df %>% 
#  group_by(phenotype) %>% 
#  slice_max(order_by = LOFC, n = 10) %>% 
#  ungroup() %>% 
#  select(newname,"EQ1":"LW3") 

#Find all names, and make the overly long ones shorter
#df3$newname
#df3 <- df3 %>% 
#  mutate(newname= replace(newname, newname=="uncharacterized LOC100644337%2C transcript variant #X1", "uncharacterized LOC100644337"),
#         newname=replace(newname, newname=="esterase E4-like%2C transcript variant X1", #"esterase E4-like"),
#         newname=replace(newname, newname=="short-chain dehydrogenase/reductase family 16C #member 6%2C transcript variant X1", "short-chain dehydrogenase/reductase family 16C"),
#         newname=replace(newname, newname=="probable phospholipase A1 magnifin%2C transcript #variant X1", "probable phospholipase A1 magnifin"),
#         newname=replace(newname, newname=="phospholipase B1%2C membrane-associated-like", #"phospholipase B12C membrane-associated-like"),
#         newname=replace(newname, newname=="26S proteasome non-ATPase regulatory subunit 1-like", "26S proteasome non-ATPase reg.s1-like"))

#Turn df3 into a numerical matrix then add names
#mx3 <- as.matrix(df3[,-1])
#rownames(mx3) <- df3$newname




#topeach10 <- pheatmap(mx3,scale="column", cluster_rows=T, cluster_cols=F,
#                      gaps_col = c(6,11),
#                      annotation_col = my_sample_col2,
#                      annotation_colors = annotation_colors,
#                      annotation_legend=F,
#                      legend = T)


### Option 3: top 50 across both within a developmental stage

#early
#df4E <- df %>% 
#  filter(phenotype=="EQ"|phenotype=="EW") %>% 
#  slice_max(order_by = LOFC, n = 50) %>% 
#  select(feature_ID,"EQ1":"EW3")   

#mx4E <- as.matrix(df4E[,-1])
#rownames(mx4E) <- df4E$feature_ID

#top50E <- pheatmap(mx4E,scale="row", cluster_rows=T, cluster_cols=T, legend=F)

#mid
#df4M <- df %>% 
#  filter(phenotype=="MQ"|phenotype=="MW") %>% 
#  slice_max(order_by = LOFC, n = 50) %>% 
#  select(feature_ID,"MQ1":"MW3")   

#mx4M <- as.matrix(df4M[,-1])
#rownames(mx4M) <- df4M$feature_ID

#top50M <- pheatmap(mx4M,scale="row", cluster_rows=T, cluster_cols=T, legend=F)

#late
#df4L <- df %>% 
#  filter(phenotype=="LQ"|phenotype=="LW") %>% 
#  slice_max(order_by = LOFC, n = 50) %>% 
#  select(feature_ID,"LQ1":"LW3")   

#mx4L <- as.matrix(df4L[,-1])
#rownames(mx4L) <- df4L$feature_ID

#top50L <- pheatmap(mx4L,scale="row", cluster_rows=T, cluster_cols=T, legend=F)


#make heatmaps into graphical objects
#top50E.g <- as.grob(top50E)
#top50M.g <- as.grob(top50M)
#top50L.g <- as.grob(top50L)

#grid.arrange(top50E.g,top50M.g,top50L.g, nrow=1, top="
               ")
# setHook("grid.newpage", NULL, "replace")
# grid.text(expression(""~bold(A)~"Early-instar larvae"), y=0.985, x=0.1, gp=gpar(fontsize=12))
# grid.text(expression(""~bold(B)~"Mid-instar larvae"), y=0.985, x=0.42, gp=gpar(fontsize=12))
# grid.text(expression(""~bold(C)~"Late-instar larvae"), y=0.985, x=0.74, gp=gpar(fontsize=12))








