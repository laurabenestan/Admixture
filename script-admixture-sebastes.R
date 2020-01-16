###################### CROSS VALIDATON RESULTS #####################

### Remove last features
rm(list=ls())
ls()

### Download librairies
library(stringr)
library(ggplot2)
library(dplyr)

### Download the cross-validation
cv <- read.table("cross_validation.txt", header=TRUE)

### Add one K-cluster column
cv$K <- seq(1,29,1)

### Select only two columns
CV <- select(cv, V4,K)

### Rename your two columns CV and K-cluster
colnames(CV) <- c("CV","K")

### Do a graph showing the cross validation results (the optimal number of clusters has the lowest cross validation error)
graph_title=""
x_title="Number of genetic groups (K)"
y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV))
graph_1+geom_line()+scale_x_continuous(breaks=seq(1,16,1))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))+
  theme_classic()

### Save the graph
ggsave("Admixture_cross-validation_fasciatus.pdf",width=7,height=5,dpi=600)
dev.off()

###################### Q ESTIMATES RESULTS #####################

### Remove last features
rm(list=ls())
ls()

### Download libraries
library(reshape2)
library(plyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

### Read file.Q with K values the most likely (the smallest CV value)
admixture <- read.table("24603snps_444ind.5.Q")
admixture <- read.table("24603snps_444ind.8.Q")

### Add one column with the individuals names and the population they belong to
id <- read.table("444ind_admixture.txt",header=FALSE)
admixture <- cbind(id,admixture)
admixture$POP <- substr(admixture$V1, 1,5)

### Rename columns
colnames(admixture) <- c("IND", "K1","K2","K3","K4","K5","POP")
colnames(admixture) <- c("IND", "K1","K2","K3","K4","K5","K6","K7","K8","POP")

### Gather the colum "K1" and "K2" into the column "ANCESTRY"
admixture_long <- melt(admixture,id.vars=c("IND","POP"),variable.name="ANCESTRY",value.name="PERC")
names(admixture_long)
class(admixture_long$ANCESTRY)
levels(admixture_long$ANCESTRY)

### subset only 50%
admixture_long_50 <- subset(admixture_long, subset=admixture_long$PERC>=0.50)

### Make a ggpot graph with ADMIXTURE results
graph_title="Stacked barplot of Admixture analysis in species"
x_title="Individuals"
y_title="Ancestry"
graph_1<-ggplot(admixture_long_50,aes(x=POP,y=PERC,fill=ANCESTRY))
graph_1+geom_bar(stat="identity")+
scale_fill_manual(values=col5, name= "K", labels=c("I","II","III","IV","V" ))+ 
  labs(y=y_title)+
  labs(x=x_title)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey", linetype="dashed"),
        axis.title.x=element_text(size=14,family="Helvetica",face="bold"),
        axis.text.x=element_text(size=6,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=14,family="Helvetica",face="bold"),
        axis.text.y=element_text(size=14,family="Helvetica",face="bold"))

### Save the map
ggsave("Sebastes_cercles.pdf",width=15,height=10,dpi=600,units="cm",useDingbats=F)

### Check percent of ancestry per sampling locations
aggregate(admixture[, 2:6], list(admixture$POP), mean)
aggregate(admixture[, 2:6], list(admixture$IND), max)

### Check the variation of % in coancestry for each genetic group
admixture[, "max"] <- apply(admixture[, 2:6], 1, max)
summary(admixture)

### Check the indivdiuals that could not be clearly attributed to one genetic cluster
admixture_subset <- subset(admixture, subset=admixture$max >= 0.5)

### Report how many individuals per cluster
admixture_cluster <- select(admixture, K1,K2,K3,K4,K5)
admixture %>% 
  gather(POP,cnt, K1:K5) %>% 
  group_by(POP) %>% 
  slice(which.max(cnt)) 

### Cretae a pop map regarding thh cluster found
admixture_subset$CLUSTER <- colnames(admixture_subset)[apply(admixture_subset,1,which.max)]
admixture_results <- select(admixture_subset, IND, CLUSTER)

### Save the files
write.table(admixture_results, 'Admixture_results_fas_K4.txt',quote=FALSE, row.names=FALSE, sep="\t", dec=".")
table(admixture_subset$CLUSTER, admixture_subset$POP)

### Subset individuals belonging to GSL
admixture$GSL <- ifelse(admixture$K4==admixture$max, 'GSL','no GSL')
gsl_subset <- subset(admixture, subset=admixture$GSL=='GSL')
table(gsl_subset$POP)
aggregate(gsl_subset[, 5], list(gsl_subset$POP), mean)

### Subset only Bonne Bay
admixture_bonne_bay <- subset(admixture, subset=admixture$POP=='fas12')
admixture_bonne_bay$RESULTS <- ifelse(admixture_bonne_bay$max==admixture_bonne_bay$K2, 'good','no')
bad_assignment <- subset(admixture_bonne_bay, subset=admixture_bonne_bay$RESULTS=='no')

### Number of individuals bellow 80%
quantile(admixture$max)
group_unknown <- subset(admixture, subset=admixture$max<0.5)
group_unknown$POP <- substr(group_unknown$IND,1,5)
table(group_unknown$POP)


############### ANALYSE RESULTS IN POP HELPER

### Install the dependency packages
install.packages(c("Cairo","devtools","ggplot2","gridExtra","gtable","tidyr"),dependencies=T)

### Install the current version of pophelper
devtools::install_github('royfrancis/pophelper', force=TRUE)

### Load library
library(pophelper)

### Check version
packageDescription("pophelper", fields="Version")

### Load population map
popmap =  read.delim("444ind_admixture.txt",header=FALSE,stringsAsFactors=F)
pop_order = sort(unique(popmap$V1))

### Load admixture files 
sfiles <- list.files(pattern = "*.Q", full.names=T)

### Basic usage
slist <- readQ(files=sfiles, filetype = "basic")

### Qlist attributes
attributes(slist)

### Dataframe attributes
head(attributes(slist[[1]]))

# read labels for STRUCTURE runs
labset <- read.table("444ind_pop.txt",header=FALSE,stringsAsFactors=F)

### Length of labels equal to number of individuals?
nrow(labset)

### Check if labels are a character data type
sapply(labset, is.character)
class(labset)

### Give a palett color
col2 <- c('blue',"red")
col5 <- c("purple","orange1","red1","chartreuse4","deeppink2")
col3 <- c('darkblue',"gray","blue",'cyan2')
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col11 <- brewer.pal(n=11, name="Spectral")

################# SPECIES ####################
### Create a qplot for K = 2 considering two species
slist1 <- alignK(slist[1]) 
plotQ(slist1,  clustercol= colorBlindBlack8  ,
      ,showsp=FALSE,grplab = labset,ordergrp=T,imgtype="pdf",
      showlegend=T, legendpos="right", legendkeysize = 6, legendtextsize = 6)

### Create a qplot from K = 2 to K =8
slist8 <- alignK(slist[10]) 
plotQ(slist8, 
      clustercol= col11, grplab = labset,ordergrp=T,
      grplabcol = "black",  grplabpos = 0.6, grplabspacer = 0.5, font="Helvetica",
      grplabangle=90, showsp=T,
      showtitle = F, titlelab = "", titlesize = 4,
      showsubtitle = F, subtitlelab = F, subtitlesize = 3,
      showlegend=T, legendpos="right", legendkeysize = 6, legendtextsize = 6,
      basesize=11,
      grplabsize=2,
      panelratio=c(1,2),
      panelspacer = 0.2, dpi = 600,imgtype="pdf")

################ Mentella
slist1 <- alignK(slist[1]) 
plotQ(slist1,clustercol= col4,grplab = labset,ordergrp=T,showlegend=T, legendpos="right", legendkeysize = 3, legendtextsize = 3,
      basesize=11,
      grplabsize=1.5, showsp=FALSE,imgtype="pdf")


################ Fasciatus
slist1 <- alignK(slist[1:3]) 
slist1 <- alignK(slist[1]) 
plotQ(slist1, 
      clustercol= col4, showsp=FALSE,grplab = labset,ordergrp=T,imgtype="pdf",showlegend=T, legendpos="right",)

