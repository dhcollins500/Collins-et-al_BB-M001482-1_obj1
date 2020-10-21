#First you must always reset R
rm(list=ls())


# LOADING PACKAGES ----

library(eulerr)  # eulerr()
library(ggplotify)  # as.ggplot()
library(grid)
library(ggpubr) #Requires cars to be installed to work


# Make a Euler diagram (proportional venn diagram).This can be done by making lists and finding the commonalities between them (as shown below). 


onlyInGeneList1<- 1:10
onlyInGeneList2<- 1:50
overlappingGenesList<- 1:5


# euler is the function with the euler package, it can read in numbers or the length of a given listobject. In this example onlyInGeneList1 is a list object, Treatment2 is the name of the first list being compared. Treatment2&Treatment3 I think is the name of the overlap list but need to make sure. overlappingGenesList is the name of value of the overlapping gene list we made before.
eulerData <- euler(c(Treatment2 = length(onlyInGeneList1), 
                     Treatment3 = length(onlyInGeneList2), 
                     "Treatment2&Treatment3" = length(overlappingGenesList)))
plot(eulerData, fill_opacity = .1)


#Make into an easy to use ggplot object  
ePlot <- as.ggplot(plot(eulerData, quantities = TRUE, main = "Treatment2&Treatment3"))

#It can also be done by inputting the numbers directly.

eulerData2 <- euler(c(Treatment2 = 5, 
                      Treatment3 = 3, 
                      "Treatment2&Treatment3" = 1))
plot(eulerData, fill_opacity = .1)

----------------------------------------------------------------------------------------------------

#Can also create custom plots  

  fit <- euler(c("A" = 10, "B" = 5, "A&B" = 3))

# Customize colors, remove borders, bump alpha, color labels white
plot(fit,
     fills = list(fill = c("red", "steelblue4"), alpha = 0.5),
     labels = list(col = "white", font = 4))

# Add quantities to the plot
plot(fit, quantities = TRUE)

# Add a custom legend and retain quantities
plot(fit, quantities = TRUE, legend = list(labels = c("foo", "bar")))

# Plot without fills and distinguish sets with border types instead
plot(fit, fills = "transparent", lty = 1:2)
  
  
---------------------------------------------------------------------------------------------------
#Below are the Euler diagrams for figure 5 in the paper
  


#Euler diagram for comparison with Cameron 60 hr old queen-destined larvae with MQ  

Cameron.60hr.vs.MQ <- euler(c(Bombus = 2346, 
                               Apis = 204, 
                               "Apis&Bombus" = 51))

ePlot.Cameron.60hr.vs.MQ<- as.ggplot(plot(Cameron.60hr.vs.MQ,labels=list(cex=0.8), quantities = TRUE))


#Euler diagram for comparison with Cameron 108 hr old queen-destined larvae with LQ  

Cameron.108hr.vs.LQ <- euler(c(Bombus = 1166, 
                             Apis = 210, 
                             "Apis&Bombus" = 8))
ePlot.Cameron.108hr.vs.LQ <- as.ggplot(plot(Cameron.108hr.vs.LQ,labels=list(cex=0.8), quantities = TRUE))

#Euler diagram for comparison with He 4 day old queen-destined larvae with MQ  
He.4d.Queen.vs.MQ <- euler(c(Bombus = 2506, 
                               Apis = 100, 
                               "Apis&Bombus" = 33))
ePlot.He.4d.Queen.vs.MQ <- as.ggplot(plot(He.4d.Queen.vs.MQ,labels=list(cex=0.8), quantities = TRUE))

#Euler diagram for comparison with He 4 day old worker-destined larvae with MW  
He.4d.Worker.vs.MW <- euler(c(Bombus = 830, 
                               Apis = 129, 
                               "Apis&Bombus" = 37))
eplot.He.4d.Worker.vs.MW <- as.ggplot(plot(He.4d.Worker.vs.MW,labels=list(cex=0.8), quantities = TRUE))

#Euler diagram for comparison with He 2 day old queen-destined larvae with MQ using just the HDEGs
He.2d.Queen.vs.MQ.HDEG <- euler(c(Bombus = 18, 
                               Apis = 40, 
                               "Apis&Bombus" = 3))
eplot.He.2d.Queen.vs.MQ.HDEG <- as.ggplot(plot(He.2d.Queen.vs.MQ.HDEG,labels=list(cex=0.8), quantities = TRUE))

#Euler diagram for comparison with He 2 day old worker-destined larvae with MW using just the HDEGs
He.2d.Worker.vs.MW.HDEG <- euler(c(Bombus = 14, 
                               Apis = 162, 
                               "Apis&Bombus" = 4))
ePlot.He.2d.Worker.vs.MW.HDEG <- as.ggplot(plot(He.2d.Worker.vs.MW.HDEG, labels=list(cex=0.8), quantities = TRUE))







Euler.panel <- ggpubr::ggarrange(ePlot.Cameron.60hr.vs.MQ,
                  ePlot.Cameron.108hr.vs.LQ,
                  ePlot.He.4d.Queen.vs.MQ,
                  eplot.He.4d.Worker.vs.MW,
                  eplot.He.2d.Queen.vs.MQ.HDEG,
                  ePlot.He.2d.Worker.vs.MW.HDEG,
          labels = c("A", "B", "C", "D", "E", "F"),
          nrow = 2, ncol = 3)


#Set WD
setwd("C:/Users/dhcollins500/OneDrive - University of East Anglia/Documents/Academic writing/My papers/Collins et al_BB-M001482-1_obj1 paper/R scripts/")

ggsave("Euler.panel.svg", Euler.panel, width = 8, height = 6)

