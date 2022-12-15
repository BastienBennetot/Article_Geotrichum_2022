# Library import ----------------------------------------------------------
library(PopGenome)
library(tidyverse)
library(stats)
library(ggplot2)
library(stringr)
library(svglite)
library(reshape2)
library(gdata)
library(xtable)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggsci)
library(MASS)
library(scales)
setwd("~/ownCloud/Thèse/gcandidum/genomic_scan")
# Definition of each population -------------------
#Strains info import
strains_info<- read.csv(file = "../../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t")
#Population and individual names
cheese1=strains_info[strains_info$population%in%"Cheese_1","strains"]
cheese2=strains_info[strains_info$population%in%"Cheese_2","strains"]
cheese3=strains_info[strains_info$population%in%"Cheese_3","strains"]
GeoB=strains_info[strains_info$population%in%"Mixed origin","strains"]
GeoC=strains_info[strains_info$population%in%"Wild","strains"]
#GeoA=c(cheese1,cheese2,cheese3)
GeoA=strains_info[strains_info$Cluster%in%"GeoA","strains"]

#WARNING : YOU HAVE TO DO THIS THE FIRST TIME
#Need to split vcf into scaffold-vcf file in scaffoldtest folder
#YOU HAVE TO LAUNCH THOSE THREE COMMAND ONCE. It will create a subfolder with a vcf for each scaffold. 
#By doing get.sum.data, you will be able to see which scaffolds have no biallelic.sites.
#They have to be removed for the analysis to work because whole.data = FALSE in sliding.window.transform function (see manual)
#VCF_split_into_scaffolds(VCF.file = "snps_pass.vcf","scaffolclib_918")
#df_genome<-readData(path = "scaffolclib_918",format = "VCF")
#get.sum.data(df_genome)

# Summary statistics mean over the genome ---------------------------------
#Once vcf file per scaffold is done using VCF_split_into_scaffold function you can read all scaffold vcf file in the right folder
df_genome<-readData(path = "scaffolclib_918",format = "VCF",include.unknown = T)
#set populations 
df_genome <- set.populations(df_genome,list(cheese1,cheese2,cheese3,GeoB,GeoC,GeoA))

##Compute population statistics------------------
df_genome<-diversity.stats(df_genome,pi=T) 
#Divide by number of sites to have a statistics per site and make a weighted mean over scaffolds to have whole genome
pi=apply(df_genome@nuc.diversity.within/df_genome@n.sites, 2, weighted.mean,df_genome@n.sites)
#Same with neutrality stat
df_genome<-neutrality.stats(df_genome,FAST = T)
theta_wat=apply(df_genome@theta_Watterson/df_genome@n.sites, 2, weighted.mean,df_genome@n.sites)
##Preparing the table of Nei's Pi and Watterson theta
diversity=rbind(pi,theta_wat)
name=colnames(diversity)#Get colnames to modify them
name=str_replace_all(string = name,pattern = "pop 1","Cheese_1")
name=str_replace_all(string = name,pattern = "pop 2","Cheese_2")
name=str_replace_all(string = name,pattern = "pop 3","Cheese_3")
name=str_replace_all(string = name,pattern = "pop 4","Mixed-origin")
name=str_replace_all(string = name,pattern = "pop 5","Wild")
name=str_replace_all(string = name,pattern = "pop 6","Cheese clade")
#Replace colnames by the modified vector
colnames(diversity)=name
diversity=formatC(diversity,digits = 3,format = "e")#Use a format with less digits and power of ten like "1.00e10"
write.table(diversity,"diversity.csv",sep="\t")#Save diversity stat^

df_genome<-F_ST.stats(df_genome,mode = "nucleotide")
fst=apply(t(df_genome@nuc.F_ST.pairwise), 2, weighted.mean,df_genome@n.sites)
dxy=apply(t(df_genome@nuc.diversity.between)/df_genome@n.sites, 2, weighted.mean,df_genome@n.sites)

pairwise=rbind(fst,dxy)
name=colnames(pairwise)
name=str_replace_all(string = name,pattern = "pop1","Cheese_1")
name=str_replace_all(string = name,pattern = "pop2","Cheese_2")
name=str_replace_all(string = name,pattern = "pop3","Cheese_3")
name=str_replace_all(string = name,pattern = "pop4","Mixed-origin")
name=str_replace_all(string = name,pattern = "pop5","Wild")
name=str_replace_all(string = name,pattern = "pop6","Cheese clade")

colnames(pairwise)=name
pairwise=as.data.frame(pairwise)
pairwise$stat=rownames(pairwise)

data_long <- gather(pairwise, population,value,colnames(pairwise)[colnames(pairwise)!="stat"], factor_key=TRUE)
data_long=data_long%>%separate(col = population,into = c("pop1","pop2"),sep = "/")

pfst=ggplot(subset(data_long,stat=="fst"),aes(pop1,pop2,fill=value,label=formatC(value,digits = 3 )))+
  geom_tile()+geom_text(color="white")+theme_bw()+ylab("")+xlab("Fst")
pdxy=ggplot(subset(data_long,stat=="dxy"),aes(pop1,pop2,fill=value,label=formatC(value,digits = 2,format ="e" )))+
  geom_tile()+geom_text(color="white")+theme_bw()+ylab("")+xlab("Nucleotive diversity between (Dxy)")

ggsave(filename ="pairwise.png" ,plot = ggpubr::ggarrange(pdxy,pfst),width = 12,height = 4)



# By sliding window -------------------------------------------------------

df_genome<-readData(path = "scaffolclib_918",format = "VCF",include.unknown = T)
#set populations 
df_genome <- set.populations(df_genome,list(cheese1,cheese2,cheese3,GeoB,GeoC,GeoA))

#get.sum.data(df_genome)
#Using sliding window
win<-7500 #Window size.
df_genome <- sliding.window.transform(df_genome,    width = win, 
                                    jump = 5000,
                                    type = 2,
                                    whole.data = FALSE)
#Number of window in the analysis
length(df_genome@region.names)

##### Computing statistical#####
df_genome<-neutrality.stats(df_genome,FAST = T)
df_genome<-F_ST.stats(df_genome,mode = "nucleotide")
df_genome<-diversity.stats(df_genome,pi=T) 
df_genome<-diversity.stats.between(df_genome)


######  Preparing data #####
#https://evolutionarygenetics.github.io/Chapter8.html
# extract nucleotide diversity and correct for window size
nd <- df_genome@nuc.diversity.within/win
#Thta watterson
theta<-df_genome@theta_Watterson/win
# make population name vector
pops<-c("Cheese_1","Cheese_2","Cheese_3","Mixed-origin","Wild","Cheese clade")
# set population names
colnames(nd) <- paste0(pops, "_pi")
colnames(theta) <- paste0(pops, "_wat")

# extract fst values
fst <- t(df_genome@nuc.F_ST.pairwise)
# extract dxy - pairwise absolute nucleotide diversity
dxy <- df_genome@nuc.diversity.between/win
#Get column names to replace it
x <- colnames(fst)
# replace all occurrences of each population with right name 
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)
x <- sub("pop3", pops[3], x)
x <- sub("pop4", pops[4], x)
x <- sub("pop5", pops[5], x)
x <- sub("pop6", pops[6], x)
# replace forward slash
x <- sub("/", "_", x)
#Replace old colnames with new one 
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

#Create a dataframe with all statistics
cam_data<-cbind.data.frame(nd,theta,fst,dxy)

#Each window is named using scaffold and position, so it's easy to retrieve these informations
#Get scaffold name
cam_data$pos.scaffold<-word(df_genome@region.names, 1)
cam_data$pos.scaffold<-as.numeric(str_replace_all(cam_data$pos.scaffold,pattern = "CCBN0100|\\.1",replacement = ""))
#Get scaffold beginning
cam_data$pos.begin<-as.numeric(word(df_genome@region.names, 2))
#Get scaffold end
cam_data$pos.end<-as.numeric(word(df_genome@region.names, 4))
#Get all position in kb
cam_data$pos.begin<-cam_data$pos.begin/1000 #in kb
cam_data$pos.end<-cam_data$pos.end/1000 #in kb
#Calculate window center for plotting
cam_data$pos.mean<- (cam_data$pos.begin+cam_data$pos.end)/2
#Get a number identifier for each window
cam_data$window<-rownames(cam_data)
#Get a statistic summary of this dataframe
print.xtable(xtable(summary(cam_data),digits = -2),"html",file="summary.html")#Get summary of this dataframe
#__________________________________________________________________________________________________________________________
# Get dataframe ready for plot
df_plot<- melt(cam_data,id.vars = c("pos.mean","pos.begin","pos.end","pos.scaffold","window"),value.name = "value")
df_plot$analysis<-str_sub(df_plot$variable,start = -3)
df_plot$analysis<-str_replace(df_plot$analysis,pattern = "_",replacement = "")
df_plot$analysis<-str_replace(df_plot$analysis,pattern = "wat",replacement = "watterson_theta")
df_plot$variable<-str_replace(df_plot$variable,pattern = "_pi|_dxy|_fst|_wat",replacement = "")
df_plot$variable<-str_replace_all(df_plot$variable,pattern = "camemberti_var_",replacement = "cam.")
df_plot$pos.scaffold<-as.numeric(str_sub(df_plot$pos.scaffold,start = -2))
df_plot<-df_plot[order(df_plot$variable),]
df_plot$variable<-as.factor(df_plot$variable)
df_plot=df_plot%>%filter(analysis%in%c("pi","watterson_theta")|
                           variable%in%c("Cheese_1_Cheese_3","Cheese_1_Mixed-origin","Cheese_1_Wild",
                                "Cheese_1_Wild","Cheese_3_Mixed-origin","Cheese_3_Wild",
                                "Mixed-origin_Cheese clade","Mixed-origin_Wild","Wild_Cheese clade"))

colpal=c("darkred","#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")
colpal2=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF")

outlier_pi<-df_plot[0,]
outlier_dxy<-df_plot[0,]


gene<-read.csv("../Geotrichum_candidum.annotations_ref_CLIB918.txt",sep = "\t",header = T)
gene=gene%>%dplyr::select(-c(gDNA,mRNA,CDS.transcript,Translation))
gene$scaffold=as.numeric(str_remove_all(gene$Contig,"CCBN01000|\\.1"))

temp = list.files(path ="../mktest/results_MKtest/to take/" ,pattern="*.tsv",full.names =T )
myfiles = lapply(temp, read.table,sep = "\t",header = T)
mktest=do.call(rbind,myfiles)
gene$MK_test="Not positively selected"
gene[gene$GeneID%in%mktest$GeneID,"MK_test"]="Positively selected"
###Small plot to see the distribution of pi and dxy

plotpi=ggplot(df_plot%>%filter(analysis=="pi"),aes(x=value,fill=variable))+
  geom_density()+scale_fill_manual(values = colpal)+
  geom_vline(data = df_plot%>%filter(analysis=="pi")%>%group_by(variable)%>%summarise(threshold=quantile(value, 0.05,na.rm = T)),mapping =aes(xintercept=threshold))+
  facet_wrap(.~variable,scales = "free")+
  ylab("Density")+
  xlab(expression(paste("Values of nucleotide diversity (",pi,") for each population")))+
  guides(fill="none")+
  theme(text = element_text(size=20))

plotdxy=ggplot(df_plot%>%filter(analysis=="dxy"),aes(x=value,fill=variable))+
  geom_density()+scale_fill_manual(values = colpal2)+
  geom_vline(data =df_plot%>%filter(analysis=="dxy")%>%group_by(variable) %>% summarise(threshold=quantile(value, 0.99,na.rm = T)),mapping =aes(xintercept=threshold))+
  facet_wrap(.~variable,scales = "free")+
  ylab("Density")+
  xlab(expression(paste("Values of absolute divergence (",D[xy],") for each pairwise comparison")))+
  guides(fill="none")+
  theme(text = element_text(size=20))
ggsave("distribution_pi.png",plot = plotpi,width = 10,height = 8)
ggsave("distribution_dxy.png",plot = plotdxy,width = 10,height = 8)



######## Overall plot by scaffold ####
for (i in unique(df_plot$pos.scaffold)){#setdiff(1:26,18)
# pi_cam_scaffold<-pi_cam[pi_cam$pos.scaffold==i,]  
# fst_cam_scaffold<-fst_cam[fst_cam$pos.scaffold==i,]
df_plot_scaffold<-df_plot[df_plot$pos.scaffold==i,]
df_plot_scaffold<-split(df_plot_scaffold, df_plot_scaffold$analysis,drop = T)
df_plot_scaffold$fst[df_plot_scaffold$fst$value<0&(!is.na(df_plot_scaffold$fst$value)),"value"]=0
df_plot_scaffold$pi<-droplevels.data.frame(df_plot_scaffold$pi)

df_pi_outlier<-df_plot_scaffold$pi %>% group_by(variable) %>% filter(quantile(value, 0.05,na.rm = T)>value)
df_pi_outlier=df_pi_outlier %>%
  group_by(pos.scaffold,window)%>%
  mutate(any_geoc=any(variable%in%"Wild"))%>%
  filter(any_geoc==F)%>%filter(variable%in%c("Cheese_1","Cheese_2","Cheese_3"))
df_pi_outlier$any_geoc=NULL

df_dxy_outlier<-df_plot_scaffold$dxy %>% group_by(variable) %>% filter(quantile(value, 0.99,na.rm = T)<value)



####pi
p1<-ggplot(data =df_plot_scaffold$pi,
           aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1.5,colour="black",
            aes(group=variable)) +
  geom_line(size=1)+
  scale_colour_manual(values=colpal)+
  geom_point(data = df_pi_outlier ,
             mapping = aes(y = value,x=pos.mean,color=variable),
             size=2,color="black")+
  theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  ggtitle(paste0("Scaffold N°", i))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),na.value=0,
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l",outside = T)+
  coord_cartesian(clip = "off")+
  xlab("Position (kb)")+
  ylab("Nucleotide diversity per site")+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=20),
        legend.key = element_rect(fill = "grey80"),
        legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),
        legend.title=element_blank())
###THETHA WATTERSON
p2<-ggplot(data =df_plot_scaffold$watterson_theta,
           aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1.5,colour="black",
            aes(group=variable)) +
  geom_line(size=1)+
  scale_colour_manual(values = colpal)+
  theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),na.value=0,
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l",outside = T)+
  coord_cartesian(clip = "off")+
  xlab("Position (kb)")+
  ylab("Watterson theta per site")+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=20),
        legend.key = element_rect(fill = "grey80"),
        legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),
        legend.title=element_blank())


#Dxy
p3<-ggplot(data =df_plot_scaffold$dxy,
           aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1)+
  geom_point(data = df_dxy_outlier ,
             mapping = aes(y = value,x=pos.mean,color=variable),
             size=2,
             color="black")+
  scale_color_npg()+
  theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  xlab("Position (kb)")+
  ylab("Dxy")+ 
  theme(legend.key = element_rect(fill = "grey80"),
        text = element_text(size=20),
        legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),
        legend.title=element_blank())

p4<-ggplot(data =df_plot_scaffold$fst,
           aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1)+
  scale_x_continuous(expand = c(0,0))+theme_bw()+
  scale_color_npg()+
  xlab("Position (kb)")+ 
  ylab("Fst")+
  theme(legend.key = element_rect(fill = "grey80"),
        text = element_text(size=20),
        legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),
        legend.title=element_blank())

p5=ggplot()+
  geom_rect(data = gene%>%filter(scaffold==i), 
            aes(xmin = Start/1000, xmax = Stop/1000, ymin = 0, ymax = 1,fill = MK_test)
                )+
  scale_fill_manual(values = c("black","red"))+
  scale_x_continuous(expand = c(0,0),
                     limits = range(df_plot_scaffold$fst$pos.mean))+
  theme_bw()+
  ylab("Genic region")+
  xlab("Position (kb)")+
  theme(axis.text.y=element_blank(),
        text = element_text(size=20),
        axis.ticks.y=element_blank())

ggsave(plot =egg::ggarrange(p1,p2,p3,p4,p5,ncol =  1,heights = c(rep(4,4),1)),filename = paste0("scaffold",i,".jpg"),width=25, height=15,dpi = 300)

#Save outlier windows for pi and fst
df_pi_outlier=df_pi_outlier%>%group_by(pos.begin,pos.scaffold)%>%mutate(variable=paste0(variable, collapse = ";"),value=max(value))%>%distinct(pos.begin,pos.scaffold,variable, .keep_all = T)
outlier_pi<-rbind.data.frame(outlier_pi,df_pi_outlier)


outlier_dxy<-rbind.data.frame(outlier_dxy,df_dxy_outlier)
}
#__________________________________________________________________________________________________________________________
outlier_genomic_scan=rbind.data.frame(outlier_pi,outlier_dxy)
outlier_genomic_scan$scaffold_name<- "CCBN0100000"
outlier_genomic_scan[outlier_genomic_scan$pos.scaffold<10,]$scaffold_name<-paste(outlier_genomic_scan[outlier_genomic_scan$pos.scaffold<10,]$scaffold_name,"0",sep = "")
outlier_genomic_scan$scaffold_name<-paste(sep = "",outlier_genomic_scan$scaffold_name,outlier_genomic_scan$pos.scaffold )
outlier_genomic_scan$scaffold_name=paste0(outlier_genomic_scan$scaffold_name,".1")
outlier_genomic_scan$begin=outlier_genomic_scan$pos.begin*1000
outlier_genomic_scan$end=outlier_genomic_scan$pos.end*1000
write_csv2(x = outlier_genomic_scan,file =   "outlier_genomic_scan.csv")


window_position<-data.frame(tot=df_genome@region.names)
window_position$scaffold<-word(window_position$tot,1)
window_position$begin<-word(window_position$tot,2)
window_position$end<-word(window_position$tot,4)
window_position$scaffold<-as.numeric(str_sub(window_position$scaffold,start = -2))
window_position$tot<-NULL
window_position$window<-rownames(window_position)
write_csv2(x = window_position,file  = "window_position.csv")


#Test if window of low pi are duplicated more than excpected
#Single sampling
sum(sample(4936,size =494 , replace = T, prob = NULL)%>%duplicated())

vector_sampling=c()
#1000 repeats
for (i in 1:10000) {
  new=sum(sample(4936,size =494 , replace = T, prob = NULL)%>%duplicated())
  vector_sampling=c(vector_sampling,new)
}
hist(vector_sampling)
#P.value of 23 being higher than expected by random
sum(vector_sampling>23)/10000










