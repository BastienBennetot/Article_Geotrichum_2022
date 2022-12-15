#Script by Bastien BENNETOT bastien.bennetot@u-psud.fr

# Setup of the environnment -----------------------------------------------
#Set your working directory
setwd("~/ownCloud/Thèse/gcandidum/cnv")
#A function to display thousands in plots using a k notation
ks <- function (x) { number_format(accuracy = 1,
                                   scale = 1/1000,
                                   suffix = "k",
                                   big.mark = ",")(x) }
# Function for %not in%
`%nin%` = Negate(`%in%`)
# Need more digits otherwise the e^5 notation fucked up everything
options("scipen"=100, "digits"=4)
#Needed library
`%nin%`=Negate(`%in%`)
library(ggplot2)
library(gdata)
library(tidyr)
library(plyr,include.only = "rbind.fill")
library(dplyr)
library(reshape2)
library("ape") 
library("phangorn")
library("ggtree")
library("adegenet")
library("ggplot2")
library("grid")
library("gridExtra")
library("ggrepel")
library("svglite")
library("ggnetworx")
library(scales)
library(stringr)
library(gdata)
library("FactoMineR")
library("factoextra")
library(caTools)
library(egg)
library(MASS)

palette_geo=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

# Import and formating ---------------------------------
temp = list.files(path = "data/ref_LMA244",full.names = T,pattern="*CNVs.p.value.txt")
#Keep name of the ref for later 
ref=unlist(str_split(string = temp[1],pattern = "/"))[2]

myfiles = lapply(temp, read.table,header = T,sep = "\t")
names(myfiles)<-temp
df<-do.call(gdata::combine,myfiles)

#df$chr<-as.numeric(str_remove_all(df$chr,pattern = "CCBN01000|\\.1")) #Keep only numbers in scaffold column
df$chr<-as.numeric(str_remove_all(df$chr,pattern = "QQZM010|\\.1")) #Keep only numbers in scaffold column

df$length<- df$end-df$start # Calculate length of each CNV

df$strains<-str_remove_all(df$source,pattern = "data/ref_LMA244/|_sorted_score_pos_markdup_q1.bam_CNVs.p.value.txt")

#Remove cnv with non significant p-value
#df<-df[df$WilcoxonRankSumTestPvalue<0.05,]
#df<-df[df$KolmogorovSmirnovPvalue<0.05,]
 
df[df$copy.number>6,]$copy.number<-6# Reducing cnv number to a maximum 

## Information about strains -----------------------------------------------
pop<-read.csv2("../../pop_info/Gcandidum summary - Feuille 1.tsv",sep='\t') # Import a csv with strains information
pop<-pop %>% add_count(population,name = "popsize") #Add a count of strains per species
df<-merge(df, pop, by.x = "strains") # Merge both dataframe

strains_order=scan("../../pop_info/order_strains_geotrichum", character())
strains_order=rev(str_remove(strains_order,pattern = "CLIB918"))

df=df%>%filter(df$strains%in%strains_order)
##Store csv file of cnv---------------------------
temp=df 
temp$source=NULL
temp$WilcoxonRankSumTestPvalue=NULL
temp$KolmogorovSmirnovPvalue=NULL
write.csv2(row.names = F,x =temp,file = "CNV_dataframe.csv")#Write the dataframe to a file

#Get a summary stat of cnv-----------------
summary_stat<-df%>%
  group_by(population,status)%>%
  summarise(number_cnv=n(),mean_length=mean(length),sd_length=sd(length))
summary_stat<-df%>%
  group_by(population,strains,status)%>%
  summarise(number_cnv=n())%>%
  group_by(population,status)%>%
  summarise(mean_number_cnv=mean(number_cnv))

# Split : one line per window ----------------------------------------------------------------------
# For each cnv create a column called window containing all windows contained in this cnv span
df$window<-apply(df,1,function(x) paste(seq(from=x['start'],to=x['end'],by=500),collapse = "|") )
df$window<-gsub("\\|*\\w*$", "", df$window)# This line remove the last number which is the end of the last window
# Create a column "window" for each cnv : contains start position of all 500pb windows that concern this CNV. Using pipe as a separator
# For exemple a CNV that start at 500 and end at 2000 give a column window with "500|1000|1500"
df<- df%>%
  mutate(window=strsplit(as.character(window),'\\|'))%>%
  unnest(window)
# Split all CNV line to have one line per window per CNV
df$window_start<-as.numeric(df$window)
df$window_end<-as.numeric(df$window)+500
# At this point we have for each line the border of the window and still have the information about border of the CNV region
df$strain_cluster<-paste(df$strains,df$population,sep = "_")# Add a column that gives information of strains and population

# Clustering analysis -------------------
df_clust<-reshape2::dcast(df, chr+window_start+window_end~strain_cluster,value.var = "copy.number")
#Create a matrix with columns (after the 3rd) corresponding to strains and each line a 500pb window
df_clust[is.na(df_clust)] <- 1 #Each NA was a void in previous dataframe which means in fact no deletion and no duplications. Number variation is 1 
df_clust<-df_clust[,-c(1:3)] # Remove the first columns that are not strains CNV
d<-dist(t(df_clust),method = "euclidean") # Calculating euclidean distance of the matrix. Each strains has an ordered vector of CNV number that reflects its identity
# Making a NJ tree from the euclidean distance
nj_clust<-NJ(d)
#Prepare tip label of the tree
vsttemp<-sapply(str_split(nj_clust$tip.label,pattern = "_"), `[`, 2)
nj_clust$vsttemp<-vsttemp
ggtree(nj_clust)+geom_tiplab(size=3) ##Tree plot

tree<-hclust(d,method="ward.D2") #Clustering on the distance matrix calculated previously
plot(tree, hang = -1, cex = 0.6)
rect.hclust(tree, k =7, border = 2:5)# Add rectangle for 4 cluster
new_order<-unique(df$strains)[tree$order] #Store the order of strains to order heatmap later
##Use ML order----------------
new_order=strains_order
#Heatmap clustered in both CNV and strains
#vsttemp<-heatmap(as.matrix(df_clust))


# PCA -----------------------
df_pca<-as.data.frame(t(df_clust))
df_pca$population=sub(".*_", "", rownames(df_pca))
df_pca$strains<-sub("_.*", "", rownames(df_pca))
pca1<-PCA(df_pca[,-which(colnames(df_pca)%in% c(  "population","strains" ))],scale.unit = TRUE, ncp = 3, graph = F)
df_pca$pc1 <- pca1$ind$coord[, 1] # indexing the first column
df_pca$pc2 <- pca1$ind$coord[, 2]  # indexing the second column
p<-ggplot(data = df_pca, aes(x = pc1, y = pc2,color=population,label=strains))+
  geom_text(aes(label=strains))+geom_point()+
  xlab(label = paste("PC1",round(pca1$eig[1,"percentage of variance"],digits = 2),"%"))+
  ylab(label =paste("PC2",round(pca1$eig[2,"percentage of variance"],digits = 2),"%"))
ggsave(plot = p,filename = paste0(ref,"_PCA_CNV.png"),units="cm", width=30, height=15, dpi = 300)
#Reordering strains
df$strains<-as.factor(df$strains)
df<-droplevels(df)
missing_strains=as.character(unique(df$strains[df$strains%nin%new_order]))
new_order=c(new_order,missing_strains)
df$strains<-reorder(df$strains,new.order = new_order)

#Remove cnv loss that concerns almost all strains, they are just reads that goes away in genome and not in FM013 reference
number_strains<-length(unique(df$strains))
df<-df%>%group_by(chr,window_start,window_end,copy.number)%>%mutate(strains_concerned=n()/number_strains)
df<-df%>%filter(!(strains_concerned>0.98&copy.number==0))

# Import scaffold borders for plots later
#borders<-read.table(file = "563_pacbio_tgs_gapfiller.fasta.fai",sep = "\t",header =F,col.names = c("chr","length","cumul","NA1","NA2"))
borders<-read.table(file = "GCA_013365045.1_LMA-244_clib_genomic.fasta.fai",sep = "\t",header =F,col.names = c("chr","length","cumul","NA1","NA2"))
borders$chr<-as.numeric(str_remove_all(borders$chr,pattern = "QQZM010|\\.1"))

# HEATMAP plot for each scaffold------------------------------------------
name<-c()
df$population=as.factor(df$population)
for (i in sort(unique(df$chr))){
#Dataframe restricted to this scaffold  
df_scaffold<-df[df$chr==i,]
borders_scaffold<-borders[borders$chr==i,]
## Big CNV discovery script --------------------------
sliding_window_size<-1#Multiply this number by the window length (it was 500pb)
minimum_interet_pop<-0.60#Minimum frequency of individual in the most concerned population that we want to conserve 
minimum_size_cnv<-1 #Multiply this number by the window length (it was 500pb)


#Prepare a temporary dataframe
temp<-df_scaffold %>%
  group_by(chr,window_start,population,popsize) %>%
  summarise(count=n(),missing=popsize-n(),seq=paste0(copy.number, collapse = ","))%>%
  distinct()
# Add a column with a sequence of copy number. This sequence is not complete because copy number of 1 are missing.
# That is why we add a count of missing copy number of 1 (no cnv)
temp$seq2<-apply(matrix(temp$missing),1,function(x) paste0(replicate(x, 1), collapse = ","))
#Create a column with a sequence of 1 coresponding to missing rows of copy equal to 1 (no cnv)
temp<-temp%>%mutate(seq3=paste(seq,seq2,sep = ","))
#Get both sequence together
temp$cnv_freq<-apply(matrix(temp$seq3),1,function(x)   mean(as.numeric(strsplit(x,split = ",")[[1]])!=1))
#Calculate the frequence (per pop) of cnv in this window. Number of copy number different from 1
temp$median_cnv_variation<-apply(matrix(temp$seq3),1,function(x) median(as.numeric(strsplit(x,split = ",")[[1]])))
# Get the median copy number 
temp<-temp%>%ungroup()%>%complete(window_start,population,fill = list(median_cnv_variation=1,count=0,cnv_freq=0))
#Fill missing species for some windows, the median cnv variation is 1 because no cnv and cnvfreq equal to 0
temp$population<-as.character(temp$population)

## categorizing cnv --------------------------------------------------------
#CNV that differ between geoC and other
temp<-temp%>%
  group_by(window_start)%>%
  mutate(geoc_other_cnv=1*(  median_cnv_variation[population%in%"GeoC"] %nin% median_cnv_variation[population %in% c("Cheese_1","Cheese_2","Cheese_3")]))
#CNV that differ between GeoB and GeoA
temp<-temp%>%
  group_by(window_start)%>%
  mutate(geob_cheese_cnv=1*(median_cnv_variation[population%in%"GeoB"] %nin% median_cnv_variation[population %in% c("Cheese_1","Cheese_2","Cheese_3")]))
#CNV that differ between cheese
temp<-temp%>%
  group_by(window_start)%>%
  mutate(cheese_diff_cnv=case_when((median_cnv_variation[population%in%"Cheese_1"] %nin% median_cnv_variation[population %in% c("Cheese_2","Cheese_3")])~1,
         (median_cnv_variation[population%in%"Cheese_2"] %nin% median_cnv_variation[population %in% c("Cheese_1","Cheese_3")])~1,
        (median_cnv_variation[population%in%"Cheese_3"] %nin% median_cnv_variation[population %in% c("Cheese_1","Cheese_2")])~1
))

# Don't need anymore the species information. Keep the maximum frequency (per pop) of cnv for each window
temp<-temp%>%group_by(window_start)%>%
  summarise(geoc_other_cnv=max(geoc_other_cnv),geob_cheese_cnv=max(geob_cheese_cnv),cnv_freq=max(cnv_freq),cheese_diff_cnv=max(cheese_diff_cnv))

#Fill missing windows that correspond to windows of no interest. We will need them to smooth windows of interest
# We will smooth gaps using a sliding max
temp<-temp%>%ungroup()%>%
  complete(window_start=full_seq(window_start,500),fill=list(geoc_other_cnv=0,geob_cheese_cnv=0,cheese_diff_cnv=0))
#We create a threshold, we will consider an interesting window if at least a certain percentage of individuals in a population have a cnv
temp$beyondthreshold<-(temp$cnv_freq>minimum_interet_pop)*1
#Sliding maximum to smooth gaps between region of interest. The sliding window size parameter can be modified
temp<-temp%>%mutate(sliding_beyondthreshold=runmax(beyondthreshold,k=sliding_window_size),
                    sliding_geoc_other=runmax(geoc_other_cnv,k = sliding_window_size),
                    sliding_cheese_diff_cnv=runmax(cheese_diff_cnv,k = sliding_window_size),
                    sliding_geob_cheese=runmax(geob_cheese_cnv,k = sliding_window_size))
#Add number of scaffold currently running
temp$chr<-i

# Now that we have boolean of membership for each window to different list, what would be perfect
# would be to merge windows to have only regions from beginning to their end. That is what
# we will do in next function

#We want a list of CNV that matter because they are not carried by a low percentage of individuals in at least one population
temp$grp<-with(rle(temp$sliding_beyondthreshold), rep(seq_along(lengths), lengths))#Create an index of regions
list_cnv_important<-temp %>%
  group_by(grp) %>%#Group by regions
  mutate(Counter = n(),begin=min(window_start),end=max(window_start)+500)%>%#Keep a count of windows in each region, start and end of this region
  ungroup() %>%
  filter(sliding_beyondthreshold>0&Counter>=minimum_size_cnv)%>% # Remove regions that are not positive for our criteria
  #We use a threshold to remove some regions that are shorter than a minimum_size_cnv
  distinct(chr,begin,end)# Get a list or region of hotspot CNV


# Now we want a list of cnv that differ between geoCvsother
temp$first_dom<-temp$sliding_beyondthreshold*temp$sliding_geoc_other
temp$grp<-with(rle(temp$first_dom), rep(seq_along(lengths), lengths))
list_cnv_geoc_other<-temp %>%
  group_by(grp) %>%
  mutate(Counter = n(),begin=min(window_start),end=max(window_start)+500)%>%
  ungroup() %>%
  filter(first_dom>0&Counter>=minimum_size_cnv)%>%
  distinct(chr,begin,end)# Get a list or region of hotspot CNV

# Now we want a list of cnv that differ between biforme VS camemberti.sl
temp$second_dom<-temp$sliding_beyondthreshold*temp$sliding_geob_cheese
temp$grp<-with(rle(temp$second_dom), rep(seq_along(lengths), lengths))
list_cnv_geob_cheese<-temp %>%
  group_by(grp) %>%
  mutate(Counter = n(),begin=min(window_start),end=max(window_start)+500)%>%
  ungroup() %>%
  filter(second_dom>0&Counter>=minimum_size_cnv)%>%
  distinct(chr,begin,end)# Get a list or region of hotspot CNV

# Now we want a list of cnv that differ between biforme VS camemberti.sl
temp$cheese_diff_cnv<-temp$sliding_beyondthreshold*temp$sliding_cheese_diff_cnv
temp$grp<-with(rle(temp$cheese_diff_cnv), rep(seq_along(lengths), lengths))
list_cnv_cheese_diff<-temp %>%
  group_by(grp) %>%
  mutate(Counter = n(),begin=min(window_start),end=max(window_start)+500)%>%
  ungroup() %>%
  filter(cheese_diff_cnv>0&Counter>=minimum_size_cnv)%>%
  distinct(chr,begin,end)# Get a list or region of hotspot CNV

## Prepare the plot --------------------------------------------------------
#First part : cnv that are important
if (nrow(list_cnv_important)==0){
  p1<-ggplot()+scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    ggtitle(paste ("Scaffold N°",i),subtitle ="CNV hotspot")
} else {
  p1<-ggplot()+
    geom_rect(data=list_cnv_important,mapping =  aes(xmin=begin, xmax=end, ymin = 0, ymax = 1),fill="grey20")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name ="" )+
    ggtitle(paste ("Scaffold N°",i),subtitle ="CNV hotspot")
}
#Second part : cnv that are different between GeoC and other
if (nrow(list_cnv_geoc_other)==0){
  p2<-ggplot()+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    ggtitle(label="GeoC vs Cheese")
} else {
  p2<-ggplot()+
    geom_rect(data=list_cnv_geoc_other,mapping =  aes(xmin=begin, xmax=end, ymin = 0, ymax = 1),fill="grey20")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    ggtitle(label="GeoC vs Cheese")
}
#Third part : cnv that are different between GeoB and GeoA
if (nrow(list_cnv_geob_cheese)==0){
  p3<-ggplot()+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    ggtitle(label = "GeoB vs Cheese")
} else {
  p3<-ggplot()+
    geom_rect(data=list_cnv_geob_cheese,mapping =  aes(xmin=begin, xmax=end, ymin = 0, ymax = 1),fill="grey20")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    ggtitle(label="GeoB vs Cheese")
}

#4th part
if (nrow(list_cnv_cheese_diff)==0){
  p4<-ggplot()+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    ggtitle(label="Cheese diff")
} else {
  p4<-ggplot()+
    geom_rect(data=list_cnv_cheese_diff,mapping =  aes(xmin=begin, xmax=end, ymin = 0, ymax = 1),fill="grey20")+
    theme(axis.title.y=element_text(angle=0, vjust = 0.5))+
    scale_x_continuous(labels = NULL,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
    scale_y_continuous(labels=NULL,breaks=NULL,name = "")+
    ggtitle(label="Cheese diff")
}

#Fifth part is the heatmap
#To keep all strains we need a line in dataframe with each strains
temp_plot=df_scaffold%>% #We just need a window
  ungroup()%>%complete(window,strains,fill = list(window_start=0,copy.number=1))%>% #We add missing combination of strains with copy number equal to 1
  filter(copy.number==1)%>% # we will merge with the previous data so we don't want to have duplicated lines
  dplyr::select(strains,window,window_start,copy.number)
temp_plot=merge(temp_plot,pop,by.x = "strains")# Re add population info which was NA after complete function
temp_plot=rbind.fill(temp_plot,df_scaffold)# Merge both dataframe


p5=ggplot(data = temp_plot)+
geom_tile(data = temp_plot,aes(y = strains,x=window_start,fill=copy.number))+
 facet_grid(tree_order_cluster~.,scales = "free",space = "free_y")+
labs(fill = "Number variation")+scale_x_continuous(labels = ks,expand = c(0,0),limits = c(0,max(borders_scaffold$length)))+
  scale_fill_gradientn(colors = c("red","white" ,"darkblue"),
                       values = rescale(c(0,1,6)),
                       breaks=c(0,1,2,3,4,5,6),
                       labels=c("  0","  1","  2","  3","  4","  5",">5"))+
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(color = "black"),
        panel.spacing = unit(0, "lines"),
        axis.title.y=element_blank(),
        plot.margin=unit(c(0,0,0,0), "mm"))

p6=ggplot(data=temp_plot)+
  geom_tile(aes(y = strains,x=1,fill=population))+
  facet_grid(tree_order_cluster~.,scales = "free",space = "free_y")+
  scale_fill_manual(values = palette_geo,na.value = "white")+
  theme_void()+
  theme( strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_rect(color = "black"),
        panel.spacing = unit(0, "lines"))+
  guides(fill="none")
empty=ggplot()+theme_void()
test=read.tree("../gcandidum_snp/MLtree_IQtree_treefile.tre.txt")
test=midpoint(test)
test$tip.label= rep('',length(test$tip.label))
p7=ggtree(test)+geom_tiplab(align = T)+ theme(plot.margin=unit(c(0,0,0,0), "mm"))
#Save the whole plot
ggsave(plot = ggarrange(empty,p1,empty,empty,p2,empty,empty,p3,empty,empty,p4,empty,p7,p5,p6,
                        ncol=3,nrow=5,widths =c(5,50,1),heights = c(rep(1,4),25)),
       filename = paste0(ref,"_scaffold",i,".png"),units="cm", width=60, height=30, dpi = 300)
##Store cnv lists of interest for further analysis ---------------
temp_cnv<-gdata::combine(list_cnv_important,list_cnv_geoc_other,list_cnv_geob_cheese,list_cnv_cheese_diff)
assign(paste0("list_cnv_",i,sep = ""),temp_cnv)#Save this scaffold list with a specific name
name<-c(name,paste("list_cnv_",i,sep = ""))#Save the name associated with this scaffold


}


# Save a CNV list of interest ---------------------------------------------
all_cnv_interet<-bind_rows(mget(name),.id="scaffold")#Merge all CNV of interest lists
all_cnv_interet$scaffold<-gsub(all_cnv_interet$scaffold,pattern = "list_cnv_",replacement = "")#Make a more readable name for each scaffold
all_cnv_interet$scaffold<-as.numeric(all_cnv_interet$scaffold)
#Preparing a column with scaffold name identical to ref just in case
all_cnv_interet$scaffold_name<- "QQZM010"
all_cnv_interet[all_cnv_interet$scaffold<10,]$scaffold_name<-paste(all_cnv_interet[all_cnv_interet$scaffold<10,]$scaffold_name,"0",sep = "")
all_cnv_interet$scaffold_name<-paste(sep = "",all_cnv_interet$scaffold_name,all_cnv_interet$scaffold )
all_cnv_interet$scaffold_name=paste0(all_cnv_interet$scaffold_name,".1")
write.csv2(row.names = F,x = all_cnv_interet,file = paste0(ref,"_CNV_list_interet.csv"))


# VST --------------------------------
save<-df
#Making pairwise combination of species
pair<-combn(unique(df$population)[!is.na(unique(df$population))], 2)
all_pair<-apply(format(pair),2,paste0,sep="_",collapse="")
all_pair<-str_replace_all(all_pair,pattern = " ",replacement = "")
all_pair<-str_sub(all_pair,end = -2)
################Computing FST pairwise for each pair of population
for (i in 1:ncol(pair)){
combi<-pair[,i]# For a specific pairwise combination
df2<-df[df$population%in%combi,]#Reduce df to the species of the desired combination

### We need a dataframe with windows that contains CNV. For each windows we want number of samples for each copy number for each population
#### Prepare a dataframe that will contains rows related to copy number equal to 1 for windows that have cnv
##### Add population size and sampling size of non cnv individuals (sampling size of cnv known individuals - population size) for each copy number and population
      vsttemp<-df2%>%
        group_by(chr,window_start,population)%>%
        dplyr::summarize(effectif=popsize-n(),popsize=popsize)%>%
        distinct()
      vsttemp$copy.number<-1#Their copy number is equal to 1
#### Prepare the dataframe for individuals that contains cnv. We want sample size for each copy number and population size
      vsttemp2<-df2%>%
        group_by(chr,window_start,population,copy.number)%>%
      dplyr::summarize(effectif=n(),popsize=popsize)%>%
        distinct()
#### Bind the two dataframe. Now we have all rows needed for windows that contains cnv
vsttemp3<-rbind.data.frame(vsttemp,vsttemp2)
vsttemp3<-droplevels(vsttemp3)# Just to be safe
popu_size<-vsttemp3%>%group_by(population)%>%dplyr::summarise(popsize=popsize)%>%distinct()#Store the population size for each species for later use
vsttemp3$popsize<-NULL #Don't need population size anymore
## Add missing combination of copy.number, population for each scaffold and windows start. We want a grouping by windows start to avoid having high copy number for windows that have few 
vsttemp4<-vsttemp3%>%group_by(chr,window_start)%>%complete(copy.number,population,fill = list(effectif=0))# Fill effectif by zeros. No individuals are concerned
vsttemp4<-merge(vsttemp4,popu_size,by = "population")# Add population size
vsttemp4$freq<-vsttemp4$effectif/vsttemp4$popsize#Add frequency of this copy number for this window in the population

#Calculation of Ht
ht<-vsttemp4%>%
group_by(chr,window_start,copy.number)%>%
  summarize(ti=sum(effectif)/sum(popsize))%>%
  group_by(chr,window_start)%>%
  summarize(Ht=1-sum(ti^2))

#Calculation of Hs
hs<-vsttemp4%>%
  group_by(chr,window_start,population)%>%
  summarise(varfreq=(1-sum(freq^2))*popsize,popsize=popsize)%>%
  group_by(chr,window_start)%>%
  summarise(Hs=sum(varfreq)/sum(popsize))

final<-merge(hs,ht,by.x = c("chr","window_start")) #Merge both dataframe ht and hs
final$fst<-(final$Ht-final$Hs)/final$Ht #Calculation of fst
final$Hs<-NULL
final$Ht<-NULL
assign(paste0("final",all_pair[i],sep = "_"),final)#Store the fst of this combination of population in a specific dataframe

}

name<-paste("final",all_pair,"_",sep = "")#Create a list of name for each dataframe related to combination
all_fst<-bind_rows(mget(name),.id="pop")#Add an id column to each dataframe
all_fst$pop<-str_sub(all_fst$pop,start = 6,end=-2)#Remove the "final_"

library(ggsci)
for (i in sort(unique(all_fst$chr))){
  #Restrict to chr dataframe
  all_fst_chr<-all_fst[all_fst$chr==i,]
  borders_chr<-borders[borders$chr==i,]
  
  ## Add window of no difference between population that have fst equal to zero
  all_fst_chr<-all_fst_chr%>%group_by(pop)%>% complete(window_start=fullseq(range=c(0,max(borders_chr$length)),size=500),fill=list(fst=0,chr=i))
  
  
  window_length=10000
  all_fst_chr$window=findInterval(x =all_fst_chr$window_start ,seq(1,max(all_fst_chr$window_start),window_length))
  all_fst_chr<-all_fst_chr%>%group_by(pop,window)%>%summarise(meanfst=mean(fst))
  all_fst_chr$window=all_fst_chr$window*window_length
  
p2<-ggplot(data =all_fst_chr,aes(y = meanfst,x=window,col=pop))+
 # geom_smooth(span=0.1,alpha=0.1,method = "loess",se = FALSE,formula = "y ~ x")+
  geom_point(size=1)+
  geom_line()+
  theme_bw()+scale_color_npg()+#xlim(0,max(df_plot_chr$dxy$pos.end)+1)+
  xlab("Position (kb)")+ ylab("Fst")+ theme(legend.key = element_rect(fill = "grey80"),legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),legend.title=element_blank())

ggsave(plot = p2,filename = paste0(ref,"_fst_chr",i,".png"),units="cm", width=60, height=30, dpi = 300)
}




