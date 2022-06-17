# set up of the environment -----------------------------------------------

library(ape) 
library(phangorn)
library(ggsci)
library(ggtree)
library(ggplot2)
library(svglite)
library(ggnetworx)
library(rgl)
library(ggtreeExtra)
library(ggnewscale)
setwd("~/ownCloud/Thèse/gcandidum/gcandidum_snp")

#Load strains informations
info_strains=read.csv("../../pop_info/Gcandidum summary - Feuille 1.tsv",sep = "\t",header = T)
#Prepare the palette
palette=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")
#Read the tre file obtained by running IQtree
test=read.tree("98seq_snp.fasta.treefile.tre")

test=midpoint(test)#Midpoint root the tree
#Plot the ML tree
d <- data.frame(node=c(132,156,101,114,107), type=c( "Cheese_1","Cheese_2","Cheese_3","GeoB","GeoC"))
p=ggtree(test)%<+% info_strains+
  geom_hilight(data=d, aes(node=node, fill=type),extend=1 , gradient = TRUE,gradient.length.out=5) +
  geom_nodelab(hjust=1.2,vjust=-0.2)+
  geom_tiplab(align = T)+
  geom_tiplab(aes(subset=(to_keep_Sophie==1|other_cheese==1)),fontface = "bold",align = T)+
   scale_color_npg(na.value="black")+ geom_treescale(offset = 0.1)+
    scale_fill_manual(values=palette)

ggsave(filename = "ML_IQ_tree.png",p,width=20,height = 20)

ggsave(filename = "ML_IQ_tree.svg",p,width=20,height = 20)

# Get the order of the tree -----------------------------------------------
#It is usefull for other analysis to order strains accordingly
#test=read.tree("MLtree_IQtree_treefile.tre.txt")
d=fortify(test)
dd = subset(d, isTip)
order_tip=dd$label[order(dd$y, decreasing=TRUE)]
write.table(x = order_tip,row.names = F,col.names = F,file = "../../admixture/order_of_strains.csv")

# Try of splitstree -------------------------------------------------
#splitree=read.nexus.networx(file = "splitree.nex")

DATA <- read.dna(file = "snp.fasta", format = "fasta")
#Convert to the format type that we need :
DATA_phyDat <- phyDat(DATA, type = "DNA", levels = NULL)
DATA_dnabin <- fasta2DNAbin(file = "snp.fasta")

dist_DATA <- dist.dna(DATA, model = "JC69")
DATA_nj <- njs(dist_DATA)

net=neighborNet(dist_DATA)
write.nexus.networx(net,"splitree.nex")

#The splittree is right but angles are messy. To clean it :
#Read on splitree4 and redo equal angle
#Then save it again
splitree=read.nexus.networx(file = "splitree.nex")
#Now we can plot it
split=ggtree(splitree,layout = "equal_angle")
split$data$strains=split$data$label
split$data=merge(split$data,info_strains,by.x = "strains",all.x = T)

split=split+geom_splitnet() +geom_treescale(offset = 0.001)+
geom_tippoint(aes(fill=population),pch=21,size=5)+
  scale_fill_manual(values = palette,na.value = "white")# geom_tiplab(aes(label=label))

ggsave(file="NJnetwork_gcand.svg", plot=split,width = 10,height = 15)


# Heatmap next to ther tree ----------------------------------------------------
df=read.csv2(file = "admixture5.csv") #Plot the admixture for K=5 populations infered
df=df%>%dplyr::select(strains,V1:V5)%>%pivot_longer(!strains,names_to = "variable",values_to = "values")
palette_geo=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")
test=read.tree("98seq_snp.fasta.treefile.tre")

test=midpoint(test)#Midpoint the tree
#Plot the tree
d <- data.frame(node=c(132,156,101,114,107), type=c( "Cheese_1","Cheese_2","Cheese_3","Mixed-origin","Wild"))
p=ggtree(test)%<+% info_strains+
  geom_nodelab(hjust=1.2,vjust=-0.2)+
  geom_tiplab(align = T)+
  geom_tiplab(aes(subset=(to_keep_Sophie==1|other_cheese==1)),fontface = "bold",align = T)+
  scale_color_npg(na.value="black")+ geom_treescale(offset = 0.1)+
 # geom_text2(aes(subset=!isTip, label=node), hjust=-.3)+
  geom_hilight(data=d, aes(node=node, fill=type) ) +
  scale_fill_manual(values=palette)
p2=p+ 
  new_scale_fill() + 
  geom_fruit(data = df,geom=geom_bar,mapping=aes(y=strains, x=values,fill=variable),pwidth=0.05,
             stat="identity",offset=0.05,
             orientation="y", # the orientation of axis.
             axis.params=list(
               axis="x", # add axis text of the layer.
               text.angle=-45, # the text size of axis.
           
               hjust=0  # adjust the horizontal position of text of axis.
             ),
             grid.params=list() # add the grid line of the external bar plot.
  )+ 
  scale_fill_manual(values = palette_geo[c(2,3,1,4,5)])+guides(fill='none') +
  theme(#legend.position=c(0.96, 0.5), # the position of legend.
    legend.background=element_rect(fill=NA), # the background of legend.
    legend.title=element_text(size=7), # the title size of legend.
    legend.text=element_text(size=6), # the text size of legend.
    legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
  ) 

p2
ggsave(filename = "ML_admixture5.png",p2,width=25,height = 20)

ggsave(filename = "ML_admixture5.svg",p2,width=25,height = 20)
################################################################################




#Another tree formated plot
test=read.tree("98seq_snp.fasta.treefile.tre")
test=midpoint(test)
p=ggtree(test)%<+% info_strains+
  geom_nodelab(hjust=1.2,vjust=-0.2)+
  geom_tiplab(aes(color=population),align = T)+
  geom_tiplab(aes(subset=(to_keep_Sophie==1|other_cheese==1)),fontface = "bold",align = T)+
  scale_color_npg(na.value="black")+ geom_treescale(offset = 0.1)+
  #geom_hilight(data=d, aes(node=node, fill=type) ) +
  scale_fill_manual(values=palette)

p+ 
  new_scale_fill() + 
  geom_fruit(data = df,geom=geom_bar,mapping=aes(y=strains, x=values,fill=variable),pwidth=0.05,
             stat="identity",offset=0.05,
             orientation="y", # the orientation of axis.
             axis.params=list(
               axis="x", # add axis text of the layer.
               text.angle=-45, # the text size of axis.
               
               hjust=0  # adjust the horizontal position of text of axis.
             ),
             grid.params=list() # add the grid line of the external bar plot.
  )+ 
  scale_fill_manual(values = palette_geo[c(2,3,1,4,5)])+guides(fill='none') +
  theme(#legend.position=c(0.96, 0.5), # the position of legend.
    legend.background=element_rect(fill=NA), # the background of legend.
    legend.title=element_text(size=7), # the title size of legend.
    legend.text=element_text(size=6), # the text size of legend.
    legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
  ) 
