
# Library and env configuration -------------------------------------------
setwd("~/ownCloud/Thèse/gcandidum/Sophie")#Need to change accordingly to the location of script and data
library(ggplot2)
library(dplyr)
library(stringr)
library(rstatix)
library(tidyr)
library(ggpubr)
library(PupillometryR)
library(Hmisc)
library(ggtext)
library(grid)
library(gtable)
library(lemon)
library(kableExtra)
library(ggrepel)

# Import DATA -------------------------------------------------------------
#Informations about strains
info_strains=read.csv("../../pop_info/Gcandidum summary - Feuille 1.tsv",sep = "\t",header = T)

# Growth_media ------------------------------------------------------------
df=read.csv2("croissance_mesure.csv",sep = "\t")#Import data

#reformat data group names
df$Date=as.factor(str_replace_all(string = df$Date,pattern = "jours",replacement = "days"))
df$Date=factor(df$Date,levels =c("7 days","11 days" ,"14 days" ) )
df$Temperature=as.factor(paste(df$Temperature,"°C",sep = ""))
#Add strains informations
df=merge(df,info_strains,by.x = "NUM.ESE")
#Change the order of media factor
df$Media=factor(df$Media,levels = c("Cheese","YPD","Minimal media"))
#Remove strains that does not belong to a population and have NA in measures
df=df%>%filter(!is.na(population)&!is.na(Measure))
df=droplevels(df)
#Make a dataframe that count number of strains
df_count=df%>%group_by(population,Media,Temperature,Date)%>%tally()
#Set a palette of color
palette=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

## Prepare Stat for radial growth on media -----------------------------------------
model=lm(data = df,formula = Measure~
           Date+
           population+
           Media+
           Temperature+
           population*Date+
           population*Media+
           population*Temperature)
summary(model)
anova(model)

#Save the anova to html
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the effect of media, temperature on radial growth of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "Anova_growth.html")

#Compute post-hoc tukey test
stat=df%>%
  group_by(Media,Date,Temperature)%>%
  tukey_hsd(Measure~population)
#Save the post hoc to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the effect of media, temperature on radial growth of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "Anova_growth2.html")

#Prepare the right position of p value on graphs
stat=stat%>%
  add_xy_position()

stat2=stat %>%filter(p.adj.signif=="ns")#We keep ns in a temp file for storage
max_area=df%>%group_by(Media,Date,Temperature)%>%
  summarise(max=max(Measure)+0.01)%>%
  dplyr::select(Media,Date,Temperature,max)#We need to know the maximum value of the graph to put brackets
max_area=max_area%>%ungroup()%>%complete(Media,Date,Temperature)
stat=merge(stat,max_area,by =c("Date","Temperature","Media")) %>%
  filter(p.adj.signif!="ns")%>%#We remove unsignificant comparison
  group_by(Media,Date,Temperature)%>%
  mutate(y.position =max+row_number()*0.08*max)#Correct the automatic spacing between bracket once the unsignificative ones are removed
stat$max=NULL#Remove the max column
stat=rbind.data.frame(stat2,stat)#For the code to work we need even the non significant even if they are not shown

#Write to a file the post hoc test 
write.csv("stat_media.csv",row.names = F,x = stat%>%select(Temperature,Media,Date,group1,group2,estimate,conf.low,conf.high,p.adj,p.adj.signif))

## Plot of radial growth on media with statistical test --------------------
p=ggplot(df,aes(x =population,y = Measure))+
  #Display points in a jitter way
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  #Plot a vertical line to show standard deviation
  stat_summary(data = subset(df,other_cheese!=1),
               fun = mean,
                fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  #Plot an horizontal dashed line to show the mean
  stat_summary(data = subset(df,other_cheese!=1),
               aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+  
  theme_bw()+
  #Add p value to the graph
  stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,size=3,inherit.aes = FALSE,hide.ns =T)+
  #Add bottom text that synthetizes how much sample for each condition
  geom_text(data =  df_count,aes(x = population,y = -Inf,label=paste("n =",n)),vjust=-1)+
  #aesthetics tweaking
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  #Split data in a grid of plots between temperature condition and media tested
  facet_grid(Media+Temperature~Date,scales = "free_y")+
  #Specific color palette
  scale_fill_manual(values =palette)+
  ylab("Colony radial growth (cm)")+
  guides(fill="none")
# Some plots were empty because we did not have all combinations of media and temperature so we replace those graph by white graphs
grob <- ggplotGrob(p)
idx <- which(grob$layout$name %in% c("panel-2-3", "panel-3-2", "panel-3-3","panel-5-2","panel-5-3","panel-7-2","panel-7-3"));
for (i in idx){ grob$grobs[[i]] <- nullGrob()}  

ggsave("croissance_milieu_temperature_stat.png",plot = grob,width = 10,height = 20)


## Subplot for 10 degree radial growth on media -----------------------------
df2=subset(df,Date=="14 days"&Temperature=="10°C")
df2_count=subset(df_count,Date=="14 days"&Temperature=="10°C")

### Statistical test --------------------------------------------------------

model=lm(data = df2,formula = Measure~
           population+
           Media+
           population*Media)
anova(model)

#Compute post-hoc test
stat=df2%>%
  group_by(Media)%>%
  tukey_hsd(Measure~population)%>%
  add_xy_position()

#Compute stat position on plot
stat2=stat %>%filter(p.adj.signif=="ns")#We keep ns in a temp file for torage
max_area=df2%>%group_by(Media)%>%
  summarise(max=max(Measure)+0.01)%>%
  dplyr::select(Media,max)#We need to know the maximum value of the graph to put brackets
max_area=max_area%>%ungroup()%>%complete(Media)
stat=merge(stat,max_area,by ="Media") %>%
  filter(p.adj.signif!="ns")%>%#We remove unsignificant comparison
  group_by(Media)%>%
  mutate(y.position =max+row_number()*0.08*max)#Correct the automatic spacing between bracket once the unsignificative ones are removed
stat$max=NULL#Remove the max column
stat=rbind.data.frame(stat2,stat)#Remove the max column


### Plot --------------------------------------------------------------------
p=ggplot(df2,aes(x =population,y = Measure))+
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  stat_summary(data = subset(df2,other_cheese!=1),
               fun = mean,
                fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(data = subset(df2,other_cheese!=1),
               aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+ theme_bw()+
  geom_text(data =df2_count ,aes(x = population,y = -Inf,label=paste("n =",n)),vjust=-1)+
  stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  facet_wrap(~Media)+
  scale_fill_manual(values =palette)+
  
  ylab("Colony radial growth (cm)")+
  guides(fill="none")
ggsave("croissance_milieu_10_stat.png",plot = p,width = 15,height = 5)


## Subplot for 25 degree radial growth on media -----------------------------
df2=subset(df,Date=="7 days"&Temperature=="25°C")
df2_count=subset(df_count,Date=="7 days"&Temperature=="25°C")

### Statistical test for 25 --------------------------------------------------------

summary(model)
model=lm(data = df2,formula = Measure~
           population+
           Media+
           population*Media)
#Save the anova to hmtl
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the effect of media on radial growth of G. candidum at  25°C")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "Anova_growth25.html")

#Compute post-hoc test
stat=df2%>%
  group_by(Media)%>%
  tukey_hsd(Measure~population)

#Save post-hoc to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the effect of media on radial growth of G. candidum at 25°C")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "Anova_growth25_2.html")

#Compute p-values position
stat=stat%>%
  add_xy_position()

stat2=stat %>%filter(p.adj.signif=="ns")#We keep ns in a temp file for storage
max_area=df2%>%group_by(Media)%>%
  summarise(max=max(Measure)+0.01)%>%
  dplyr::select(Media,max)#We need to know the maximum value of the graph to put brackets
max_area=max_area%>%ungroup()%>%complete(Media)
stat=merge(stat,max_area,by ="Media") %>%
  filter(p.adj.signif!="ns")%>%#We need to know the maximum value of the graph to put brackets
  group_by(Media)%>%
  mutate(y.position =max+row_number()*0.08*max)#Correct the automatic spacing between bracket once the unsignificative ones are removed
stat$max=NULL#Remove the max column
stat=rbind.data.frame(stat2,stat)#Remove the max column

### Plot --------------------------------------------------------------------
p=ggplot(df2,aes(x =population,y = Measure))+
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  stat_summary(data = subset(df2,other_cheese!=1),
               fun = mean,
                fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(data = subset(df2,other_cheese!=1),
               aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+ theme_bw()+
  geom_text(data =df2_count ,aes(x = population,y = -Inf,label=paste("n =",n)),vjust=-1)+
  stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  facet_wrap(~Media)+
  scale_fill_manual(values =palette)+
  ylab("Colony radial growth (cm)")+
  guides(fill="none")
ggsave("croissance_milieu_25_stat.png",plot = p,width = 20,height = 7)


# Salt tolerance ----------------------------------------------------------------
df=read.csv("salinite.csv",header = T,sep = "\t",dec=",")
df=merge(df,info_strains,by.x = "NUM.ESE")
df=droplevels(subset(df,!is.na(df$population)))
df$Salt=paste(df$Salt,"% Salt",sep = "")

## Salt statistical test -------------------------------------------------
model=lm(data = df,Measure~population+Salt+Salt*population)
summary(model)
anova(model)
#Save the linear model to html
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the effect of salt content on radial growth of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "Salt.html")

#Comput post-hoc test
stat=df%>%
  group_by(Salt)%>%
  wilcox_test(Measure~population)
  
#Save post-hoc to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the effect of salt content on radial growth of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "Salt2.html")  

#Compute p-value postion on plots
stat=stat%>%add_xy_position()

stat2=stat %>%filter(p.adj.signif=="ns")
max_area=df%>%group_by(Salt)%>%
  summarise(max=max(Measure)+0.01)%>%
  dplyr::select(Salt,max)
max_area=max_area%>%ungroup()%>%complete(Salt)
stat=merge(stat,max_area,by ="Salt") %>%
  filter(p.adj.signif!="ns")%>%
  group_by(Salt)%>%
  mutate(y.position =max+row_number()*0.08*max)
stat$max=NULL
stat=rbind.data.frame(stat2,stat)
stat=stat%>%ungroup()%>%mutate(xmin=xmin,xmax=xmax)


## Salt plot----------------------------------------------
p=ggplot(df, aes(x = population, y = Measure)) + 
  scale_fill_manual(values =palette) +
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  stat_summary(fun = mean,
                fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+
  theme_bw()+
   facet_wrap(~Salt,nrow=1)  +
  ylab(label = "Colony radial growth (cm)")+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  geom_text(data = df %>%filter(!is.na(Measure))%>% group_by(Salt, population) %>% tally(), 
            aes(x = population, y = -Inf, 
                label = paste("n =",n)), 
            position = position_dodge(width = 0.9),vjust=-1, color = "black")+
  stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  guides(fill="none")
ggsave("Salt_stat.png",plot = p,width = 20,height = 7)


# Growth by milk type ---------------------------------------------
df = read.table(file = "croissance_typedelait.csv", header = T, sep = ";", dec = ",")
df = merge.data.frame(df, info_strains, by = "NUM.ESE")
df=subset(df,!is.na(df$population))
df=df%>%rename(Measure=Mean)
df$Media=str_replace(string = df$Media,pattern = "Chevre",replacement ="Goat" )
df$Media=str_replace(string = df$Media,pattern = "Brebis",replacement ="Sheep" )
df$Media=str_replace(string = df$Media,pattern = "Vache",replacement ="Cow" )
df2=subset(df,other_cheese!=1)

#Statistical linear model
model=lm(data = df,Measure~population*Media)
summary(model)
anova(model)
#We save the linear model to html
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the effect of milk species on radial growth of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "milk_type.html")

#We compute post-hoc test
stat=df%>%
  group_by(Media)%>%
  tukey_hsd(Measure~population)

#Save post-hoc test to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the effect of milk species on radial growth of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "milk_type2.html")  

#Compute p-value position
stat=stat%>%add_xy_position()
stat2=stat %>%filter(p.adj.signif=="ns")
max_area=df%>%group_by(Media)%>%
  summarise(max=max(Measure)+0.01)%>%
  dplyr::select(Media,max)
max_area=max_area%>%ungroup()%>%complete(Media)
stat=merge(stat,max_area,by ="Media") %>%
  filter(p.adj.signif!="ns")%>%
  group_by(Media)%>%
  mutate(y.position =max+row_number()*0.08*max)
stat$max=NULL
stat=rbind.data.frame(stat2,stat)
stat=stat%>%ungroup()%>%mutate(xmin=xmin,xmax=xmax)

#Do the plot
p=ggplot(df, aes(x = population, y = Measure)) + 
  scale_fill_manual(values =palette) +
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  stat_summary(data = df2,
               fun = mean,
               fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(data = df2,
               aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+ theme_bw()+
  facet_wrap(~Media,nrow=1)  +xlab(label ="Population") + 
  ylab(label = "Colony radial growth (cm)")+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  geom_text(data = df %>%filter(!is.na(Measure))%>% group_by(Media, population) %>% tally(), 
            aes(x = population, y = -Inf, 
                label = paste("n =",n)), 
            position = position_dodge(width = 0.9), color = "black",vjust=-1)+
 stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  guides(fill="none")
ggsave("Type_of_milk.png",plot = p,width = 20,height = 7)

#IRIS data for opacity measure --------------------------------
df=read.csv("iris_new_data3.csv",sep = "\t")
df=merge(df,info_strains,by.x = "NUM.ESE")
df=df %>%filter(!is.na(population))
#Linear model
model=lm(data = df,opacity~population)

#Save the linear model to html
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the effect of population on opacity of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "opacity_new.html")

#Compute post-hoc test
stat=df%>%
  tukey_hsd(opacity~population)
#Save post-hoc to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the effect of population on opacity of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "opacity_new2.html") 

#Compute pvalue position on plot
maxi_opacity=max(df$opacity)
stat=stat%>%add_xy_position()
stat=stat%>% mutate(max=maxi_opacity+0.01,y.position=0)
stat2=stat %>%filter(p.adj.signif=="ns")
stat=stat %>%
  filter(p.adj.signif!="ns")%>%
  mutate(y.position =max+row_number()*0.08*max)
stat=stat%>%ungroup()%>%mutate(xmin=xmin,xmax=xmax)

#Plot
p=ggplot(df, aes(x = population, y = opacity)) + 
  scale_fill_manual(values =palette) +
  geom_jitter(aes(fill=population),position = position_jitter(seed = 1,width=0.1),size=3,pch=21,color="black")+
  stat_summary(data = subset(df,other_cheese!=1),
               fun = mean,
               fun.min = function(x)  if((mean(x) - sd(x))<0){0}else{mean(x) - sd(x)}, 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(data = subset(df,other_cheese!=1),
               aes(ymax = ..y.., ymin = ..y..),fun = mean,geom="errorbar",linetype="dashed")+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
xlab(label ="Population") + 
  ylab(label = "Integral opacity score")+
  geom_text(data = df %>%filter(!is.na(opacity))%>% group_by(population) %>% tally(), 
            aes(x = population, y = -Inf, 
                label = paste("n =",n)), 
            position = position_dodge(width = 0.9) ,vjust=-1,color = "black")+
  stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  guides(fill="none")

ggsave("new_iris.png",plot = p,width = 10,height = 7)


# Competiton with mat of G. candidum and a dot of competitor --------------------------------------------------------------

df=read.csv("competition.csv",sep = "\t",dec = ",")
df$NUM.ESE[df$NUM.ESE=="Nabis"]<-NA
df[is.na(df$NUM.ESE),"NUM.ESE"]="control"
df_solo_competitor=df%>%filter(NUM.ESE=="control"&Competitor_species!="Debaryomyces hansenii")
df_solo_competitor=df_solo_competitor%>%group_by(Competitor_species)%>%summarise(Max=round(mean(Mean),digits = 2))

soucheok_compet=info_strains[info_strains$competition_tapis_point==1&info_strains$to_keep_Sophie==1,"NUM.ESE"]
df=droplevels(subset(df,NUM.ESE%in%c(soucheok_compet,"control")))

df=df%>%ungroup()%>%complete(NUM.ESE,Competitor,fill=list(Mean=0))
df=df%>%group_by(Competitor)%>%fill(Competitor_species,Temp,.direction = "downup")

df=merge(df,info_strains,all.x = T,by.x = "NUM.ESE")
df=merge(df,df_solo_competitor,by="Competitor_species")
df[df$NUM.ESE=="control","population"]="No G. candidum"

## Statistical linear model for competition ----------------------------------------
model=lm(df,formula = Mean~population+Competitor_species+population*Competitor_species)

#Save the linear model to html
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the competition abilities of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "competition.html")

#Compute post-hoc test
stat=df%>%
  group_by(Competitor_species)%>%
  tukey_hsd(Mean~population)

#Save post-hoc test to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the competition abilities of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "competition2.html")  


#Compute pvalue position on the plot
stat=stat%>%add_xy_position()
stat2=stat %>%filter(p.adj.signif=="ns")
max_area=df%>%group_by(Competitor_species)%>%
  summarise(max=max(Mean)+0.01)%>%
  dplyr::select(Competitor_species,max)
max_area=max_area%>%ungroup()%>%complete(Competitor_species)
stat=merge(stat,max_area,by ="Competitor_species") %>%
  filter(p.adj.signif!="ns")%>%
  group_by(Competitor_species)%>%
  mutate(y.position =max+row_number()*0.08*max)
stat$max=NULL
stat=rbind.data.frame(stat2,stat)


## Plot competition with stat test -----------------------------------------
p=ggplot(data = df,mapping = aes(y=Mean,x=population))+
  facet_wrap(.~Competitor_species)+
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(aes(ymax = ..y.., ymin = ..y..),fun = mean,geom="errorbar",linetype="dashed")+
  stat_pvalue_manual(stat%>%filter(Competitor_species%nin%"Penicillium roqueforti"), label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  guides(fill="none")+theme_bw()+
  ylab("Radial growth of the competitor (mm)")+
  geom_text(data = df %>%filter(!is.na(Mean)) %>%group_by(Competitor_species, population) %>% tally(), 
            aes(x = population, y = -Inf, 
                label = paste("n =",n)), 
            position = position_dodge(width = 0.9),vjust=-1, color = "black")+
  scale_fill_manual(values = palette,na.value="black")+ 
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))

ggsave("competition_direct.png",plot = p,width = 15,height = 8)




# Competition in splitted petri dish: only volatiles can influence ---------------------------------------------------

df=read.csv(file = "Competition_volatiles.csv",sep = "\t",header = T,dec = "," )
df[is.na(df$NUM.ESE),"NUM.ESE"]="Control"
df=merge(x = df,y = info_strains,all.x =T,by= "NUM.ESE")
palette=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")
df[df$NUM.ESE=="Control","population"]="No G. candidum"
df=df%>%filter(Competitor!=645&!is.na(population))
# Remove rows of competitor equal to 645 because it was contaminated
df_control=df%>%
  filter(population=="No G. candidum")%>%
  group_by(Competitor,Competitor_species,population)%>%summarise(Mean=mean(Mean,na.rm = T))

#Linear model
model=lm(data = df,formula = Mean~Competitor_species+population+Competitor_species*population)
#Save anova of this model to html
tab=anova(model)%>%as.data.frame()%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  rename("P.value"=`Pr(>F)`)%>%
  kable(caption = "Anova table of the linear model about the competition abilities of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "competition_VC.html")

#Compute post-hoc test
stat=df%>%
  group_by(Competitor_species)%>%
  tukey_hsd(Mean~population)

#Save post-hoc to html
tab=stat%>%
  mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
  kable(caption = "Post-hoc test of the linear model about the competition abilities of G. candidum")%>%
  kable_classic(html_font = "arial")
save_kable(x = tab,file = "competition_VC2.html")  

#Compute pvalue position on plot
stat=stat%>%add_xy_position()
stat2=stat %>%filter(p.adj.signif=="ns")
max_area=df%>%group_by(Competitor_species)%>%
  summarise(max=max(Mean,na.rm = T)+0.01)%>%
  dplyr::select(Competitor_species,max)
stat=merge(stat,max_area,by ="Competitor_species") %>%
  filter(p.adj.signif!="ns")%>%
  group_by(Competitor_species)%>%
  mutate(y.position =max+row_number()*0.08*max)
stat$max=NULL
stat=rbind.data.frame(stat2,stat)

p=ggplot(data = df,mapping = aes(y=Mean,x=population))+
  facet_wrap(.~Competitor_species)+
  geom_jitter(aes(fill=population),size=3,width=0.1,pch=21,color="black")+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(aes(ymax = ..y.., ymin = ..y..),fun = mean,geom="errorbar",linetype="dashed")+
  stat_pvalue_manual(stat, label = "p.adj", tip.length = 0.01,inherit.aes = FALSE,hide.ns =T)+
  guides(fill="none")+theme_bw()+
  ylab("Radial growth of the competitor (mm)")+
  geom_text(data = df %>%filter(!is.na(Mean)) %>%group_by(Competitor_species, population) %>% tally(), 
            aes(x = population, y = -Inf, 
                label = paste("n =",n)), 
            position = position_dodge(width = 0.9),vjust=-1, color = "black")+
  
  scale_fill_manual(values = palette,na.value="black")+ 
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))

ggsave(p,filename = "competition_volatiles.png",width = 15,height = 8)

#A command to combine All html -----------------------
#pandoc -f html -t html -o statistics_phenotype.html Anova_growth.html Anova_growth2.html Anova_growth25.html Anova_growth25_2.html Salt.html Salt2.html milk_type.html milk_type2.html opacity.html opacity2.html competition.html competition2.html competition_VC.html competition_VC2.html
#
#

