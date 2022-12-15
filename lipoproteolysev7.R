
# Environment setup ---------------------------------------------------------------


setwd("~/ownCloud/Thèse/gcandidum/lipoproteolyse")
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpp)
library(rstatix)
library(kableExtra)
library(ggpubr)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(multcomp)

palette=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

# Info import -------------------------------------------------------------
strains=read.csv("strains.csv")
strains[strains$LMA.ID=="","LMA.ID" ]=strains[strains$LMA.ID=="","ESE.ID" ]
strains[strains$ESE.ID=="","ESE.ID" ]=strains[strains$ESE.ID=="","LMA.ID" ]
strains$strains=strains$LMA.ID

gcand_summary=read.csv2("../../pop_info/Gcandidum summary - Feuille 1.tsv",sep = "\t")
gcand_summary$strains=str_remove(gcand_summary$strains,pattern = "-")
gcand_summary$population=str_replace(gcand_summary$population,pattern = "-",replacement = "_")


# Lypolysis ---------------------------------------------------------------
df_lipo=read.csv("lypolysis.csv")

#Reformat the data and add strains informations
df_lipo=df_lipo%>%pivot_longer(ESE00516:LMA1849, names_to = "strains", values_to = "Measure")
df_lipo=merge(df_lipo,strains,by = "strains")
df_lipo$strains=df_lipo$ESE.ID
df_lipo$strains=str_replace_all(df_lipo$strains,"LMA[0]*",replacement = "LMA")
df_lipo=df_lipo%>%group_by(T.,Day,X,strains)%>%summarise(Measure=mean(Measure,na.rm=T))


df_lipo=merge(df_lipo,gcand_summary,by="strains")
df_lipo=df_lipo%>%filter(!is.na(population))

df_lipo_count=df_lipo%>%group_by(T.,Day,X,population)%>%tally()
df_lipo=df_lipo%>%pivot_wider(names_from=X,values_from = Measure)

df_lipo=df_lipo%>%  rename(Measure=Lipolysis)
df_lipo=df_lipo%>%
  rename(Temperature=T.)


## Statistical analysis --------------------------------------------------------------------
#Discretize data
df_lipo=df_lipo%>%mutate(discrete_measure=(Measure!=0)*1)
#Make a binomial model
model <- glmer(discrete_measure ~population+(1 |Day )+(1 | Temperature), data = df_lipo, family = binomial)
model
summary(model)
#Save model to html
tab_model(model,file = 'Anova_lipo.html')

# summary(glht(reg, mcp(population="Tukey")))
## Plot --------------------------------------------------------------------

p=ggplot(df_lipo,aes(x =population,y = Measure))+
  geom_jitter(aes(fill=population),size=3,width=0.1,height = 0,pch=21,color="black")+
  stat_summary(fun = mean,
               fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+  theme_bw()+
  geom_text(data =df_lipo_count ,aes(x = population,y = -1,label=paste("n =",n)),size=4)+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  facet_grid(Temperature~Day, labeller = label_both,scales="free")+
  scale_fill_manual(values =palette,na.value = "white")+

  ylab("Length of lysis (mm)")+
  guides(fill="none")

ggsave("lipolyse_lypolysis.png",plot = p,width = 10,height = 7)


# Proteolysis -------------------------------------------------------------


df_proteo=read.csv("proteolysis.csv")

#Reformat data and add strains informations
df_proteo$X=df_proteo$Measure
df_proteo$Measure=NULL
df_proteo=df_proteo%>%pivot_longer(ESE00516:LMA1849, names_to = "strains", values_to = "Measure")
df_proteo=merge(df_proteo,strains,by = "strains")
df_proteo$strains=df_proteo$ESE.ID
df_proteo$strains=str_replace_all(df_proteo$strains,"LMA[0]*",replacement = "LMA")
df_proteo=df_proteo%>%group_by(T.,Day,X,strains)%>%mutate(var=paste(Measure,collapse = "|"))
df_proteo=df_proteo%>%group_by(T.,Day,X,strains,var)%>%summarise(Measure=mean(Measure,na.rm=T))

df_proteo=merge(df_proteo,gcand_summary,by="strains")
df_proteo=df_proteo%>%filter(!is.na(population))

df_proteo_count=df_proteo%>%group_by(T.,Day,X,population)%>%tally()
 df_proteo=df_proteo%>%pivot_wider(names_from=X,values_from = Measure)%>%rename(Measure=Proteolysis)


df_proteo=df_proteo%>%
  rename(Temperature=T.)

## Statistical analysis ------------------
#Discretize data
df_proteo=df_proteo%>%mutate(discrete_measure=(Measure!=0)*1)
#Binomial model
model <- glmer(discrete_measure ~population+(1 |Day )+(1 | Temperature), data = df_proteo, family = binomial)
model
summary(model)
#Save model to html
tab_model(model,file = 'Anova_proteo.html')
#Do post hoc test
post_test=summary(glht(model, mcp(population="Tukey")))
sink("post_hoc_proteo.txt")
print(post_test)
sink()  # Save post hoc test

## Plot --------------------------------------------------------------------
p=ggplot(df_proteo,aes(x =population,y = Measure))+
  geom_jitter(aes(fill=population),size=3,width=0.1,height = 0,pch=21,color="black")+
  stat_summary(fun = mean,
                fun.min = function(x) max(0,mean(x) - sd(x)), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "linerange")+
  stat_summary(aes(ymax = ..y.., ymin = ..y..),
               fun = mean,geom="errorbar",
               linetype="dashed")+  theme_bw()+

  geom_text(data =df_proteo_count ,aes(x = population,y = -1,label=paste("n =",n)),size=4)+
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,face = "bold",color="black"))+
  scale_fill_manual(values =palette,na.value = "white")+
  facet_grid(Temperature~Day, labeller = label_both,scales="free")+
  ylab("Length of lysis (mm)")+
  guides(fill="none")
ggsave("proteolyse_proteolysis.png",plot = p,width = 10,height = 7)

