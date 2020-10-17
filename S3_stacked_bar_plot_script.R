#DavidHCollins
#04 August 2017
#This script is to make a stacked barplot of the read types in my mRNA-seq data

########################################################################################################################

#First you must alwys reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(dplyr)
library(ggplot2)



#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor

#Set wording directory
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Academic writing/My papers/Collins et al_BB-M001482-1_obj1 paper/R scripts/")


#Enter table

library(readr)
transcriptstatus <- read_csv("transcriptstatus.csv")
View(transcriptstatus)

#Reorder the phenotype factors so they are in the correct order for the x axis.
transcriptstatus$Phenotype <- factor(transcriptstatus$Phenotype,levels = c("MQ", "MW", "LQ", "LW"))
transcriptstatus$Age <- factor(transcriptstatus$Age, levels = c("mid", "late"), 
                  labels = c("Mid", "Late"))

#Proportional bar plot

ggplot(transcriptstatus, #the data you want to use
       aes(x = factor(Phenotype), #x axis
       fill = reorder(Status,Order)))+ #y axis ordered by the "order" column
      geom_bar(position = "fill") + #makes figure a proportional barplot
      scale_fill_grey(name="Status")+ #colour scale and name of the legend
      xlab("Phenotype")+ #x label
      ylab("Proportion") #y label

# Count bar plot

#barplot without gaps
ggplot(transcriptstatus, #the data you want to use
       aes(x = factor (Phenotype), #x axis
           fill = reorder(Status,Order))) + #y axis ordered by the "order" column
        geom_bar() + #type of chart
        scale_fill_grey(name="Status") + #colour scale and name of the legend
  xlab("Phenotype")+ #x label      
  ylab("Count") #y label

#barplot with gaps but also additional problems
(stacked_plot <- ggplot(transcriptstatus, #the data you want to use
       aes(x = factor (Phenotype), #x axis
           fill = reorder(Status,Order))) + #y axis ordered by the "order" column
  geom_bar() + #type of chart
  scale_fill_grey(name="Status") + #colour scale and name of the legend
  facet_grid(~ Age,scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(),
        axis.line= element_line(),
        legend.position = "right",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(1,"line")) +
  xlab("Phenotype")+ #x label      
  ylab("Count") )#y label

ggsave("stacked_plot.png", stacked_plot, width = 8, height = 5)
ggsave("stacked_plot.pdf", stacked_plot, width = 8, height = 5)
