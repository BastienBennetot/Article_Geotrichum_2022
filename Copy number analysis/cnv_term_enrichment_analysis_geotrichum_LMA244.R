library(tidyr)
library(dplyr)
require(data.table) 
library(stats)
setwd("~/ownCloud/Thèse/gcandidum/cnv")

df<-read.csv2("ref_LMA244_CNV_list_interet.csv")
ref="ref_LMA_244"
df<-df[df$source!="list_cnv_important",]


gene<-read.csv("../Geotrichum_candidum_LMA_244.annotations.txt",sep = "\t",header = T)
gene=gene%>%dplyr::select(-c(gDNA,mRNA,CDS.transcript,Translation))
gene$chr<-as.numeric(str_remove_all(gene$Contig,pattern = "QQZM010|\\.1")) #Keep only numbers in scaffold column



#gene$scaffold_name=gene$Contig
gene$begin=gene$Start
gene$end=gene$Stop
df$chr<-as.character(df$chr)
gene$chr<-as.character(gene$chr)

setDT(gene)  ## convert to data.table without copy
setDT(df)##Convert to data.table type
setkey(gene,chr, begin, end)#Set key, it will play by scaffold_name on begin and end
cnv_gene<- foverlaps(df,gene ,type="any")# Godly-tier function that check if some interval in one dataframe overlap interval in another dataframe and merge if it is the case. Moreover you can do that by scaffold. One word A W E S O M E
write.csv(cnv_gene,file = paste0(ref,"_cnv_annotated.csv"),row.names = F)

df$length=df$end-df$begin
mean(df$length)

df%>%group_by(source)%>%tally()#CNV regions
cnv_gene%>%filter(is.na(GeneID))%>%group_by(source)%>%tally()#Non genic regions
temp=cnv_gene%>%filter(!is.na(GeneID))
temp%>%group_by(source)%>%tally()# Number of genes within CNV regions
temp%>%filter(source=="list_cnv_geoc_other")%>%group_by(InterPro)%>%tally()#hypothetical protein count


compute_enrichment=function(cnv_gene,gene,list_to_filter,func_choice,sep){
  cnv_gene<- cnv_gene%>%
    mutate(!!func_choice := strsplit(as.character(!!sym(func_choice)),sep))%>%
    unnest(!!func_choice)
  gene<- gene%>%
    mutate(!!func_choice:=strsplit(as.character(!!sym(func_choice)),sep))%>%
    unnest(!!func_choice)

cnv_gene2=cnv_gene%>%filter(source==list_to_filter)
gene_interest=na.omit(unique(cnv_gene2$GeneID))
gene=gene%>%mutate(selected=GeneID%in%gene_interest)
tested_term=gene%>%group_by_at(func_choice)%>%mutate(n_CNV_GO=sum(selected),n_GO=n())
tested_term=tested_term%>%ungroup()%>%distinct(!!!syms(func_choice),n_CNV_GO,n_GO)%>%filter(n_CNV_GO!=0)
tested_term$n_CNV<-nrow(cnv_gene2)# Need the number of gene concerned of interest
tested_term$n_gene_total<-nrow(gene) #Total number of gene in the genome
tested_term=tested_term%>%rowwise()%>%mutate(pvalue=fisher.test(matrix(c(n_CNV_GO,n_CNV-n_CNV_GO,n_GO-n_CNV_GO,n_gene_total-n_GO-n_CNV+n_CNV_GO),nrow=2),alternative='greater')$p.value)
tested_term=tested_term%>%ungroup()%>%mutate(significant=pvalue<(0.05/n()))
tested_term$list=list_to_filter
return(tested_term)
}

library(foreach)
test=foreach(i=unique(cnv_gene$source),.combine = rbind.data.frame)%do%
compute_enrichment(cnv_gene = cnv_gene,gene=gene,list_to_filter = i,func_choice = "InterPro",sep=";")
write.table("enrichment_CNV_LMA244.csv",x = test,row.names = F,sep = "\t")
