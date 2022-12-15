#Prepare and plot histogram of population inference for different K------------------------------------
## Initialization-----------------------
#Set working directory
setwd("~/ownCloud/Thèse/admixture_gcandidum/")
#Uploading of libraries
`%nin%` = Negate(`%in%`)
library(gdata)
library(reshape2)
library(ggplot2)
library(ggsci)
library(stringr)
library(dplyr)
library(zoo)
library(tidyr)

#Import informations about strains
strains_info<-read.csv(file = "../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t")
# Column names needed : strains , population ,

#Import strains name used in NGS admix and rename them if needed
pop<-read.table(header = F,"strains_order")# This list strains contains strains used in the analysis
colnames(pop)<-c("strains")
pop$strains=str_replace(pop$strains,pattern =  "NT12_ESE00548",replacement = "NT12")
pop$strains=str_replace(pop$strains,pattern =  "VTTC4559_ESE00549",replacement = "VTTC4559")


#### Define your range of K that you used, Kmin and Kmax
kmin<-2
kmax<-6
# Read inferred admixture proportions file
temp2<-read.table("results_qopt/NGSadmix_K2_repeat_1.qopt")
temp3<-read.table("results_qopt/NGSadmix_K3_repeat_1.qopt")
temp4<-read.table("results_qopt/NGSadmix_K4_repeat_1.qopt")
temp5<-read.table("results_qopt/NGSadmix_K5_repeat_1.qopt")
temp6<-read.table("results_qopt/NGSadmix_K6_repeat_1.qopt")
temp7<-read.table("results_qopt/NGSadmix_K7_repeat_1.qopt")
# temp8<-read.table("results_qopt/NGSadmix_K8_repeat_1.qopt")
# temp9<-read.table("results_qopt/NGSadmix_K9_repeat_1.qopt")
# temp10<-read.table("results_qopt/NGSadmix_K10_repeat_1.qopt")


#It would be highly nice and beautiful if you had a strains order of an ML tree. This helps everything to look better
#This can be done swiftly by doing those commands :
#test=read.tree("MLtree_from_IQtree_treefile.tre.txt")
#d=fortify(test)
#dd = subset(d, isTip)
#order_tip=dd$label[order(dd$y, decreasing=TRUE)]
#Then you can just save this order somewhere
#write.table(x = order_tip,row.names = F,col.names = F,file = "life/is/great/order_of_strains.csv")




# First step : Find the order of strains that order the better the plot -------------------------

#Create a big dataframe with all data for each K and population ID
df<-cbind.data.frame(mget(paste("temp",kmin:kmax,sep = "")))
#Get a column name vector
colonne_name<-paste("K",rep(kmin:kmax, times=kmin:kmax),"_",colnames(df),sep = "")
#Rename the dataframe with NGSadmix analysis accordingly
colnames(df)<-colonne_name
# Create the distance matrix
d <- dist(df, method = "euclidean")
# Hierarchical clustering
hc1 <-hclust(d, method = "ward.D2" )
# Plot the obtained dendrogram to have an idea of how it looks
plot(hc1, cex = 0.6, hang = -1)

#You can choose either to follow the hclust order of strains or follow an order defined by you (for example an ML tree)

##1st choice Use the hclust order ----------------------------------------------------
#Unquote next paragraph if you want to use this one

# pop$strains <- factor(pop$strains, levels = pop[hc1$order,"strains"])
# pop$strains<-pop[order(pop$strains),] #reordering levels of strains to follow the MLtree_order

##2nd choice Use your own order ----------------------------------------------------

ml_order=read.table("order_of_strains.csv")$V1#Import the order
pop$strains <- factor(pop$strains, levels = ml_order)#Reorder the factor of pop$strains



# Second step : Identify which colors goes where and in which order----------------------------

## First attribute correspondance between different K ----------------------

# The difficulty is that for different K you have a different number of colors and a different order because columns are independantly ordered
df2<-as.data.frame(t(df))#We need the transposed of the dataframe
# Compute dissimilarity matrix
d <- as.matrix(dist(df2, method = "euclidean"))
save<-c()#This will be a vector where we save correspondance between higher and lower K

#For each K populations we gradually see which one correspond the most to the previous ones
for (i in kmin:(kmax-1)) {
#We define coordinate of interest inside the distance matrix "d" for K and K+1
n<-i
beginr<-n*(n+1)/2-kmin+2
endr<-beginr+n
beginc<-n*(n-1)/2-kmin+2
endc<-beginc+n-1
test<-d[beginc:endc,beginr:endr]

#Now we have to select for each K+1 which K is the closest
save_duplicates<-list()
for (j in 1:n) {
  save_duplicates[[j]]<-apply(test, 1,function(x){which(x==sort(x,partial=j)[j])})
}
save_duplicates<-do.call(rbind,save_duplicates)#This give the order that matches
#In case two new populations (at K) are both having the same K+1 population as the closest population we have to make a choice
save_duplicates[1,]<-replace(save_duplicates[1,], duplicated(save_duplicates[1,]), NA)#If that is the case, we replace by a NA the duplicated assignation
save_duplicates<-na.locf(save_duplicates,fromLast=T)#And we move to the next closest

#Just in case there is multiple duplicated possibilities we redo this from 2 to n
for (j in 2:n) {
  save_duplicates[1:j,]<-t(apply(save_duplicates[1:j,], 1, function(x) replace(x, duplicated(x), NA)))
  save_duplicates<-na.locf(save_duplicates,fromLast=T)
}

numeric_vector<-save_duplicates[1,]
#If there is still duplicated that were impossible to resolve, we just do random attributions
while (any( duplicated(numeric_vector))) {
  numeric_vector[duplicated(numeric_vector)]<-sample(which(1:(n+1)%nin%numeric_vector),1)
}
save<-c(save,paste(numeric_vector,collapse = ",",sep = ""))
}


## Attribute colors to the different K based on the attribution of correspondance between K -------------------------------------------------------------------------

### Choose your palette -----------------------------------------------------

#From Paul Tol: https://personal.sron.nl/~pault/
#Tol_muted <- c("blue", '#117733', "gray60", "#00CC99", '#999933',"gray20",  "Maroon 1", '#AA4499', '#DDDDDD','#E41A1C','black','cyan','darkred','darkgoldenrod1','#44AA99')
palette_geo=c("Maroon 1","#00CC99" ,"darkred" , "gray20","blue", "gray60")
#scales::show_col(palette_geo)#If you want to see how your colors look

### Get palette in the right order for each K -----------------------------------------------------
#Define the used palette
palette=palette_geo[1:kmax]
save_palette<-paste(palette_geo[1:kmax],sep = "",collapse = ",")#Save the palette in a string way

#For each K from K to K-1 attribute colors from the previous one
#The trick is that each number that define a population over all K is defined from a K to another K+1. I mean that the "1" for K=2 may not be the same as "1" in K=4 because it just say the coordinate between K
for (j in (kmax-1):kmin) {
  palette<- palette[  as.numeric(str_split(save[j-1],pattern = ",",simplify = T)[1,])  ]#Attribute color for the known K
  save_palette<-c(save_palette,paste(palette,collapse = ",",sep = ""))#Save the new color palette name with one less color
}
save_palette<-rev(save_palette)#It was in decreasing order so we have to reverse it

### Make a vector with each new color when you increase K  -----------------------------------------------------
newcolor<-c()
for (i in 1:(kmax-kmin)) {
  temp<-str_split(save_palette[i+1],pattern = ',',simplify = T)[1,]
  templag<-str_split(save_palette[i],pattern = ',',simplify = T)[1,]
  newcolor<-c(newcolor,which(temp%nin%templag))
}


# Plot histograms of different K ----------------------------------
## Single K plot function definition------------------
plot_histo <- function(pop,temp,palette) {
  #Bind NGSadmix results with the strains name
  df<-cbind.data.frame(pop,temp)
  #Melt different K populations assignations to different rows instead of columns
  df<- melt(df, id.vars = c("strains"))
  #How much K do we have ?
  k_number<-length(unique(df$variable))
  #Grab the right color for your plot using the K number
  mycol<-str_split(save_palette[k_number-kmin+1],pattern = ",",simplify = T)
  #Arrange by values then number the order of membership for each population by strains
  df<-df%>%arrange(value)%>%group_by(strains)%>%mutate(ordre= row_number())

  #Just a specific case if that is NOT your graph for your lowest K
    if(k_number-kmin!=0){
    p<-ggplot(data=df, aes(x=strains, y=value,fill=variable,group=ordre))+
      geom_bar(stat="identity")+
      theme_void()+ 
      theme(axis.text.y = element_text( size=12))+
      scale_fill_manual(label ="new color",name=paste("K = ",k_number,sep = ""),
                        values = mycol,
                        breaks = levels(df$variable)[newcolor[k_number-kmin]])+
      theme(legend.title.align = 0.5,
            plot.margin = unit(c(1, 5.5, 1, 5.5), "points"),
            text=element_text(size=21))
    }
  #Just a specific case if that is your graph for your lowest K. Because we want all color to be said to be newcolor
  if(k_number==kmin){
    p<-ggplot(data=df, aes(x=strains, y=value,fill=variable,group=ordre)) +
      geom_bar(stat="identity")+
      theme_void()+
      theme(axis.text.y = element_text( size=12),
            plot.title = element_text(hjust = 0.5,size=21))+
      ggtitle("Membership of each strains to K clusters")+
      scale_fill_manual(label =rep("new color",length(levels(df$variable))),
                        name=paste("K = ",k_number,sep = ""),
                        values = mycol)+
      theme(legend.title.align = 0.5,
            plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
            text=element_text(size=21))
  }
  
  return(p) #Return the plot for the K. Good job !  
}
## Single K as a whole function definition------------------

#A quick function if you wanted to plot a single K because you are a jackass. You want strains to be labelled on the bottom
plot_histo_alone <- function(pop,temp,palette) {
  #Bind NGSadmix results with the strains name
  df<-cbind.data.frame(pop,temp)
  #Melt different K populations assignations to different rows instead of columns
  df<- melt(df, id.vars = c("strains"))
  #How much K do we have ?
  k_number<-length(unique(df$variable))
  #Grab the right color for your plot using the K number
  mycol<-str_split(save_palette[k_number-kmin+1],pattern = ",",simplify = T)
  #Arrange by values then number the order of membership for each population by strains
  df<-df%>%arrange(value)%>%group_by(strains)%>%mutate(ordre= row_number())
  #Reverse the order of strains factor
  df$strains=forcats::fct_rev(df$strains)

    p<-ggplot(data=df, aes(y=strains, x=value,fill=variable,group=ordre)) +
      geom_bar(stat="identity")+
      theme(axis.text.y = element_text( size=12),
            plot.title = element_text(hjust = 0.5,size=21))+
      ggtitle("Membership of each strains to K clusters")+
      scale_fill_manual("K populations",values = palette_geo)+
      theme(legend.title.align = 0.5,
            plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
            text=element_text(size=21))
  return(p)  
}

# Time for real shit ----------------

# Get ID and pop info for each individual tab delimited 
strains_info<-read.csv(file = "../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t")
# 
# ##If you want a csv to remember this great analysis in addition to the plot ----------------------------
# #Here an exemple for K equal to 5 populations
# #Bind NGSadmix results with the strains name
# df5<-cbind.data.frame(pop,temp5)
# #Melt different K populations assignations to different rows instead of columns
# df5<- melt(df5, id.vars = c("strains"))
# #How much K do we have ?
# k_number<-length(unique(df5$variable))
# #Grab the right color for your plot using the K number
# mycol<-str_split(save_palette[k_number-kmin+1],pattern = ",",simplify = T)
# #Arrange by values then number the order of membership for each population by strains
# df5<-df5%>%arrange(value)%>%group_by(strains)%>%mutate(ordre= row_number())
# #Reverse the order of strains factor
# df5$strains=forcats::fct_rev(df5$strains)
# #Long to wide
# df5=df5%>%select(strains,variable,value)%>%spread(key = variable,value = value)
# #Add strains informations
# df5=merge(df5,strains_info,by = "strains",all =  T)
# 
# #A specific part because I had some clonal strains and did not kept them in the NGSadmix analysis but can still spread the results afterward within a clonal group
# df5=df5%>%select(strains,V1:V5,Clonal.group)%>%filter(!is.na(Clonal.group))#Remove not well annotated clonal group, in fact polyploid strains that were in strains_info
# df5=df5%>%group_by(Clonal.group)%>%fill(V1:V5,.direction = "downup")%>%ungroup()#Spread NGSadmix results within clonal group

#Save it as a csv
#write.csv2("../gcandidum/gcandidum_snp/admixture5.csv",x = df5)


## Plot for each K ---------------------------------------------------------
p2<-plot_histo(pop,temp2,save_palette)
p3<-plot_histo(pop,temp3,save_palette)
p4<-plot_histo(pop,temp4,save_palette)
p5<-plot_histo(pop,temp5,save_palette)
p6<-plot_histo(pop,temp6,save_palette)
p7<-plot_histo(pop,temp7,save_palette)
#p8<-plot_histo(pop,temp8,save_palette)
#p9<-plot_histo(pop,temp9,save_palette)
#p10<-plot_histo(pop,temp10,save_palette)


## additional plot for informations about strains --------------------------
### Strains and real population assignation -------------------------------------------------------------------------

# Get ID and pop info for each individual tab delimited 
pop2=merge(pop,strains_info)
pop3=pop2[order(pop2$strains),]

#In fact the palette geo has to be in a different order between K assignation and and here because we order in a different way K pop and real populations
palette_geo=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

#A plot with real population informations and strains name
leg_strains<-ggplot(data = pop3)+
  geom_tile(aes(fill=population,x=strains,y=1),color="black")+
  scale_fill_manual(name="Population infered on NJ tree with SNP data",values =palette_geo,na.value = "white")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"),
        text=element_text(size=21))+
  scale_y_continuous( expand = c(0, 0)) +
  scale_x_discrete(labels = pop3$strains)

### Adding informations about milk from which cheese was made and st --------
#A palette for each origin of milk
pal_milk=c("goldenrod1","purple","firebrick","darkgreen","tan4","darkblue","black")
#A plot to label origin of milk for each strains
leg5=ggplot(data = pop3)+
  geom_tile(aes(fill=type_fromage,x=strains,y=1),color="black")+theme_void()+
  geom_text(aes(label=Environment,y=1,x=strains),angle=90)+
  scale_fill_manual(name="Milk type",values = pal_milk)+
  theme(plot.margin = unit(c(1, 5.5, 1, 5.5), "points"),
        text=element_text(size=21),
        legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1,override.aes = list(size = 1)))

## Make the whole plot and save it-----------------
#Get the full list of plots
plot_list<-mget(c(paste0("p",kmin:kmax),"leg5","leg_strains"))
#Save it
ggsave(plot =egg::ggarrange(plots=plot_list,ncol = 1,  heights = c(rep(3,kmax-kmin+1),5,1)  ),
       filename = "NGSadmix_results_ward.D2_euclidean.png",units="cm", width=35, height=25, dpi = 300)



# Checking likeliness of different run for same k : multiplot ---------------------------------------------------------------


for (f in kmin:kmax) {# Look for each K population
##Import inputs -------------
K<-f
#Get all filenames for this K
temp = list.files(path = "results_qopt",pattern=paste("NGSadmix_K",K,"_repeat_[0-9]*.qopt",sep = ""),full.names = T)
#How much repeats do you have ?
repeats=length(temp)
#Read them
myfiles = lapply(temp, read.table)
#Bind them in a dataframe, each column being a repeat
test<-do.call(cbind,myfiles)
#Rename columns
colnames(test)<-rep(1:repeats,each=K)

## Hierarchical clustering dendogram plot----------------
# The goal is to understand which population n is similar to another one between repeats
#Compute the distance matrix
d <- dist(t(test), method = "euclidean")
#Hierarchical clustering
hc1 <- hclust(d, method = "ward.D2" )
# Plot the obtained dendrogram
#This gives you an idea of how similar your repeats are. You expect to have K branch in this dendogram that are K populations with as many leafs as repeats
jpeg(filename = paste("tree_K",K,sep = ""),units="cm", width=70, height=10, res = 300)
plot(hc1, cex = 0.6, hang = -1,cex=0.3)
rect.hclust(hc1,k = K,border = "blue")
rect.hclust(hc1,k = K+1,border = "red")
dev.off()

## Hierarchical clustering admixture plot----------------
#reordering levels of strains to follow the ML order so that is readable
ml_order=read.table("order_of_strains.csv")$V1
#Reorder factors accordingly
pop$strains <- factor(pop$strains, levels = ml_order)


### get correspondance of populations for each repetition of histogram--------------

#Transpose the dataframe
df2<-as.data.frame(t(test))
#Rename rows accordingly
rownames(df2)<-paste("rep_",rep(1:repeats,each=K),"_V_",rep(1:K,times=100),sep="")
#Compute the dissimilarity matrix
d <- as.matrix(dist(df2, method = "euclidean"))
save<-c()
for (i in 1:(repeats-1)) {#Start from a repeat and compare it to the next one
  n<-i
  beginc<-n*K +1
  endc<-n*K+K
  beginr<-1
  endr<-K
  temp<-d[beginc:endc,beginr:endr]#Get the right subset of the distance matrix for this comparison

  #Get the correspondance using the closest one
  save_duplicates<-list()
  for (j in 1:K) {
    save_duplicates[[j]]<-apply(temp, 1,function(x){which(x==sort(x,partial=j)[j])})
  }
  #In case there is duplicates we can attribute to the second closest or third or so on
  save_duplicates<-do.call(rbind,save_duplicates)
  save_duplicates<-t(apply(save_duplicates,1,function(x){replace(x,duplicated(x),NA)}))
  save_duplicates<-na.locf(save_duplicates,fromLast=T)

  numeric_vector<-save_duplicates[1,]
  
# If there is no perfect closest then we attribute it randomly
  while (any( is.na(numeric_vector))) {
    numeric_vector[is.na(numeric_vector)]<-sample(which(1:K%nin%numeric_vector),1)
  }
  save<-c(save,paste(numeric_vector,collapse = ",",sep = ""))
}

###Attribute color to the correspondance vector-------------------
#Big palette
Tol_muted<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#000000' , '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff','#fffac8')
palette<-Tol_muted[1:K]
save_palette<-paste(Tol_muted[1:K],sep = "",collapse = ",")

for (j in 1:(repeats-1)){#Replace numbers by the right color
  numeric_vector<-as.numeric(str_split(save[j],pattern = ",",simplify = T)[1,])
  palette_temp<-palette[numeric_vector]
save_palette<-c(save_palette,paste(palette_temp,collapse = ","))
}

##Get a better order between repeats so it looks nicer-----------------------
#First we will need to reorder each K population for each n repeats so we will be able to make the distance matrixes of each K populations
#This is the better K order for each repeats
better_order<-c(1:K,as.numeric(str_split(paste0(save,collapse = ","),pattern = ",",simplify = T)[1,]) )
#Add number of rows when this repeat happened
analysis_repeat<-rep(0:(repeats-1),
                     each=K)*K
#The new order for all rows
better_order_rows<-better_order+analysis_repeat
#Reorder population within repeats
df3<-df2[better_order_rows,]
#We will label rownames to make it nicer later
df3$rep<-rep(1:100,each=K)#What repeat ?
df3$var<-rep(1:K,times=100)# What K population
rownames(df3)<-paste("rep",df3$rep,"var",df3$var)
# For each K population we will make a distance matrix of between repeats and normalize it by its max value
list_distance_matrix<-list()
for (i in 1:K) {
  matrix_temp<-dist(df3[df3$var==i,colnames(df3)%nin%c("var","rep")], method = "euclidean")
 matrix_temp<- matrix_temp/max(matrix_temp)
list_distance_matrix[[i]]<-matrix_temp
}
#Now that we have a distance matrix of each population, we want to add them to get an overall picture
new_matrix<-Reduce('+', list_distance_matrix)
#Now we can cluster between repeats of the analysis
hc2 <- hclust(new_matrix, method = "ward.D2" )
#If you want to plot the obtained dendrogram
#plot(hc2, cex = 0.6, hang = -1)

#Here is the plot order to make it nicer. Great job!
order_for_plot<-hc2$order


## Time to plot ------------------------------------------------------------

#We have what we needed
#- Correspondance between colors of different repeats
#- An order of repeats that optimize the visualization

#We have to transpose the whole matrix. Operations on rows are easier
new<-as.data.frame(t(test))
#Rename by K population infered and repeats
new$rep<-rep(1:100,each=K)
new$var<-rep(1:K,times=100)
#Wide to long
new<-new%>%gather(individual,value,-c(rep,var))
#A bit of reformating of strains names that were columns
new$individual<-as.numeric(str_sub(start=2,new$individual))
new$individual<-as.factor(new$individual)
levels(new$individual)<-pop$strains
#Reorder strains factor to have a nice order depending on the ml_order. This will make the plot looks better
new$individual <- factor(new$individual, levels = ml_order)
#Arrang by the value and add a number for which the membership of the strains is higher
new<-new%>%arrange(value)%>%group_by(individual,rep)%>%mutate(ordre= row_number())
new$var<-as.factor(new$var)

## Prepare the list of each repetition histogram plot-----------------
myplot_list<-list()
for (i in 1:repeats) {
#Get the color for this repeat
  my_col<-str_split(save_palette[i],pattern = ",",simplify = T)
#Make the plot for this repeat and store it in the list at its index number
myplot_list[[i]]<-ggplot(data=new[new$rep==i,], aes(x=individual, y=value,fill=var,group=ordre)) +#,group=ordre
  geom_bar(stat="identity")+
  scale_fill_manual(values=my_col)+
  theme_void()+
  theme(legend.title.align = 0.5,
        axis.title.y = element_text(),
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
        text=element_text(size=21))+guides(fill="none")+
  ylab(paste("rep",i))
}


##Facultative part : to add informations about strains ------------------
pop$strains <- factor(pop$strains, levels = ml_order)#Reorder accordingly to the ML
pop2=merge(pop,strains_info)
pop3=pop2[order(pop2$strains),]
#Strains name and real population informations
leg_strains<-ggplot(data = pop3)+
  geom_tile(aes(fill=population,x=strains,y=1),color="black")+
  scale_fill_manual(name="Population infered on NJ tree with SNP data",values =pal_d3("category10")(10))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"),
        text=element_text(size=21))+
  scale_y_continuous( expand = c(0, 0)) +
  scale_x_discrete(labels = pop3$strains)
##Final plot ---------------------------
myplot_list_reordered<-myplot_list[order_for_plot]
#If you want to add strains name and real population as a legend 
myplot_list_reordered[[101]]<-leg_strains
#Save the plot of your whole repeats
ggsave(filename = paste("Plot_rep_","K",K,sep = ""),
       device = "png",
       plot = egg::ggarrange(plots = myplot_list_reordered,ncol = 1),
       width=10,
       height=50,
       dpi = 300,
       limitsize = F,
       bg="white")
}

