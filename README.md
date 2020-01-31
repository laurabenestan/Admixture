Fishery genomics of Sebastes: investigating population structure using ADMIXTURE
--------

Using ADMIXTURE, a program customed for SNP datasets only, we aimed to define the genetic units present in two fish species: *Sebastes faciatus* and *S. mentella*. ADMIXTURE estimates individual ancestries by efficiently computing maximum likelihood estimates in a parametric model. To better understand this program, read the manual page [Alexander 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146885/).
Those two species are curently fished in some regions of the Atlantic Canada while in othe regions the stocks are endangered. 
Using a wide dataset of SNPs markers may help us to better delineate the stock structure than with the use of microsatellites.

### Run Admixture in bash
Go into the folder where your bed file is. Add the path of this file in your terminal by typing:
```{r, engine = 'bash', eval = FALSE}
for K in 1 2 3 4 5; do admixture --cv=10 -B2000 -j8 nameofyourfile.bed $K | tee log${K}.out; done
```

### Collect the cross validation information obtained from the log files
```{r, engine = 'bash', eval = FALSE}
grep -h CV log*.out>cross_validation.txt
done
```

### Take the right order for individual id using the tfam file information
```{r, engine = 'bash', eval = FALSE}
cut -f 1 nameofyourfile.tfam > id_admixture.txt
done
```

# 1. Select the optimal nuber of clusters
#### Remove last features
```{r}
rm(list=ls())
library(stringr)
library(ggplot2)
library(dplyr)
```
#### Download the **cross-validation** results you have previously created via bash command.
```{r}
cv <- read.table("cross_validation.txt")
```

#### Analyze the **cross-validation** results
Then, add a K-cluster column indicating the number of K you test and select only two columns of interest, CV and K.
```{r}
cv$K <- c(1,2,3,4,5)  
CV <- select(cv, V4,K)
```
Rename your two columns CV and K-cluster
```{r}
colnames(CV) <- c("CV","K")
```

create a summary of your **cross-validation** results
```{r}
colnames(CV) <- c("CV","K")
```

### Visualising the results
Do a graph showing the cross validation results (the optimal number of clusters has the lowest cross validation error)
```{r}
graph_title="Cross-Validation plot"
x_title="K"
y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV))
graph_1+geom_line()+scale_x_continuous(breaks=c(1,2,3,4,5))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))
```
Save the graph
```{r}
ggsave("Admixture_cross-validation.pdf",width=7,height=5,dpi=600)
dev.off()
```
![Admixture cross-validation results.](Cross-validation.png)


# 2. Analyse Q estimates results

### Remove last features
```{r}
rm(list=ls())
ls()
```

### Download libraries
```{r}
library(reshape2)
library(plyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
```

### Read file.Q with K values the most likely (the smallest CV value)
```{r}
admixture <- read.table("24603snps_860ind.2.Q")
```

### Add one column with the individuals names and the population they belong to
```{r}
id <- read.table("444ind_admixture.txt",header=FALSE)
admixture <- cbind(id,admixture)
admixture$POP <- substr(admixture$V1, 1,5)
```

### Rename columns
```{r}
colnames(admixture) <- c("IND", "K1","K2","POP")
```

### Gather the colum "K1" and "K2" into the column "ANCESTRY"
```{r}
admixture_long <- melt(admixture,id.vars=c("IND","POP"),variable.name="ANCESTRY",value.name="PERC")
names(admixture_long)
class(admixture_long$ANCESTRY)
levels(admixture_long$ANCESTRY)
```

### Subset only the individuals showng more than 50% of ancestry with one genetic cluster
```{r}
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
```

### Save the map
```{r}
ggsave("Sebastes_cercles.pdf",width=15,height=10,dpi=600,units="cm",useDingbats=F)
```

### Check percent of ancestry per sampling locations
```{r}
aggregate(admixture[, 2:6], list(admixture$POP), mean)
aggregate(admixture[, 2:6], list(admixture$IND), max)
```

### Check the variation of % in coancestry for each genetic group
```{r}
admixture[, "max"] <- apply(admixture[, 2:6], 1, max)
summary(admixture)
```

### Check the indivdiuals that could not be clearly attributed to one genetic cluster
```{r}
admixture_subset <- subset(admixture, subset=admixture$max >= 0.5)
```

### Report how many individuals per cluster
```{r}
admixture_cluster <- select(admixture, K1,K2,K3,K4,K5)
admixture %>% 
  gather(POP,cnt, K1:K5) %>% 
  group_by(POP) %>% 
  slice(which.max(cnt)) 
```

### Cretae a pop map regarding thh cluster found
```{r}
admixture_subset$CLUSTER <- colnames(admixture_subset)[apply(admixture_subset,1,which.max)]
admixture_results <- select(admixture_subset, IND, CLUSTER)
```

### Save the files
```{r}
write.table(admixture_results, 'Admixture_results_fas_K4.txt',quote=FALSE, row.names=FALSE, sep="\t", dec=".")
table(admixture_subset$CLUSTER, admixture_subset$POP)
```

### Subset individuals belonging to GSL
```{r}
admixture$GSL <- ifelse(admixture$K4==admixture$max, 'GSL','no GSL')
gsl_subset <- subset(admixture, subset=admixture$GSL=='GSL')
table(gsl_subset$POP)
aggregate(gsl_subset[, 5], list(gsl_subset$POP), mean)
```

### Number of individuals bellow 80%
```{r}
quantile(admixture$max)
group_unknown <- subset(admixture, subset=admixture$max<0.5)
group_unknown$POP <- substr(group_unknown$IND,1,5)
table(group_unknown$POP)
```

############### ANALYSE RESULTS IN POP HELPER

### Install the dependency packages and library
```{r}
install.packages(c("Cairo","devtools","ggplot2","gridExtra","gtable","tidyr"),dependencies=T)
library(pophelper)
```

### Check version
```{r}
packageDescription("pophelper", fields="Version")
```

### Install the current version of pophelper
```{r}
devtools::install_github('royfrancis/pophelper', force=TRUE)
```

### Load population map
```{r}
popmap =  read.delim("444ind_admixture.txt",header=FALSE,stringsAsFactors=F)
pop_order = sort(unique(popmap$V1))
```

### Load admixture files
```{r}
sfiles <- list.files(pattern = "*.Q", full.names=T)
```

### Basic usage
```{r}
slist <- readQ(files=sfiles, filetype = "basic")
```

### Qlist attributes
```{r}
attributes(slist)
```

### Dataframe attributes
```{r}
head(attributes(slist[[1]]))
```

# Read labels for STRUCTURE runs
```{r}
labset <- read.table("860ind_pop.txt",header=TRUE,stringsAsFactors=F)
```

### Length of labels equal to number of individuals?
```{r}
nrow(labset)
```

### Check if labels are a character data type
```{r}
sapply(labset, is.character)
class(labset)
```

### Give a palett color
```{r}
col2 <- c('blue',"red")
```

### Create a qplot for K = 2 considering two species
```{r}
slist1 <- alignK(slist[1]) 
plotQ(slist1,  clustercol= col2,
      ,showsp=FALSE,grplab = labset,ordergrp=T,imgtype="pdf",
      showlegend=T, legendpos="right", legendkeysize = 6, legendtextsize = 6)
```

