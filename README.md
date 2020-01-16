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

### Swith to R to analyze the ADMIXTURE results 
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
![Admixture cross-validation results.](Admixture.png)






