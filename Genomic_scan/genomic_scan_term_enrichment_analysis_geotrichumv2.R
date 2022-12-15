
# Env setup ---------------------------------------------------------------
library(tidyr)
library(dplyr)
require(data.table) 
library(stringr)
library(stats)
library(foreach)
setwd("~/ownCloud/Thèse/gcandidum/genomic_scan")


# Data import -------------------------------------------------------------
df<-read.csv2("outlier_genomic_scan.csv")
genomic_annotation<-read.csv("../Geotrichum_candidum.annotations_ref_CLIB918.txt",sep = "\t",header = T)
genomic_annotation=genomic_annotation%>%dplyr::select(-c(gDNA,mRNA,CDS.transcript,Translation))

genomic_annotation$scaffold_name=genomic_annotation$Contig
genomic_annotation$begin=genomic_annotation$Start
genomic_annotation$end=genomic_annotation$Stop
df$scaffold_name<-as.character(df$scaffold_name)
genomic_annotation$scaffold_name<-as.character(genomic_annotation$scaffold_name)

#Remove Cheese_2 from pi enrichment analysis because it did not have enough diversity to be analysed
df=df%>%filter(variable!="cheese2")

#Unnest the pi which can be attributed to multiple populations 
df=df%>%mutate(variable=strsplit(as.character(variable),';'))%>%
  unnest(variable)

# Foverlaps of position and annotated genes -------------------------------

setDT(genomic_annotation)  ## convert to data.table without copy
setDT(df)##Convert to data.table type
setkey(genomic_annotation,scaffold_name, begin, end)#Set key, it will play by scaffold_name on begin and end
genomic_scan_gene<- foverlaps(df,genomic_annotation ,type="any")# Godly-tier function that check if some interval in one dataframe overlap interval in another dataframe and merge if it is the case. Moreover you can do that by scaffold. One word A W E S O M E
write.csv(genomic_scan_gene,file = "genomic_scan_annotated.csv",row.names = F)

temp=genomic_scan_gene%>%filter(!is.na(GeneID))
temp%>%group_by(variable,analysis)%>%tally()
df%>%group_by(variable,analysis)%>%tally()


# Get number of genes involved in analysis for the fisher exact test --------
number_of_pi=genomic_scan_gene%>%filter(analysis=="pi")%>%filter(!is.na(GeneID))%>%
  group_by(variable)%>% tally()
number_of_dxy=genomic_scan_gene%>%filter(analysis=="dxy")%>%filter(!is.na(GeneID))%>%
  group_by(variable)%>% tally()

genomic_scan_gene=genomic_scan_gene%>%filter(!is.na(GeneID))



#genomic_scan_gene=genomic_scan_gene[grep(x = genomic_scan_gene$Product,pattern = "lysophospholipase|galactonate|protease",ignore.case = T),]

#Function to compute enrichment test
compute_enrichment=function(potentially_enriched_gene,genomic_annotation,list_to_filter,func_choice,sep,number_genes_picked){
  potentially_enriched_gene<- potentially_enriched_gene%>%
    mutate(!!func_choice := strsplit(as.character(!!sym(func_choice)),sep))%>%
    unnest(!!func_choice)
  genomic_annotation<- genomic_annotation%>%
    mutate(!!func_choice:=strsplit(as.character(!!sym(func_choice)),sep))%>%
    unnest(!!func_choice)

potentially_enriched_gene2=potentially_enriched_gene%>%filter(variable==list_to_filter)
gene_interest=na.omit(unique(potentially_enriched_gene2$GeneID))
genomic_annotation=genomic_annotation%>%mutate(selected=GeneID%in%gene_interest)
tested_term=genomic_annotation%>%group_by_at(func_choice)%>%mutate(n_target_GO=sum(selected),n_GO=n())
tested_term=tested_term%>%ungroup()%>%distinct(!!!syms(func_choice),n_target_GO,n_GO)%>%filter(n_target_GO!=0)


#tested_term$n_target<-nrow(potentially_enriched_gene2)# Need the number of gene concerned of interest
tested_term$n_target<-as.numeric(number_genes_picked[number_genes_picked$variable==list_to_filter,"n"])# Need the number of gene concerned of interest
tested_term$n_gene_total<-nrow(genomic_annotation) #Total number of gene in the genome



tested_term=tested_term%>%rowwise()%>%mutate(pvalue=fisher.test(matrix(c(n_target_GO,n_target-n_target_GO,n_GO-n_target_GO,n_gene_total-n_GO-n_target+n_target_GO),nrow=2),alternative='greater')$p.value)
tested_term=tested_term%>%ungroup()%>%mutate(significant=pvalue<(0.05/n()))
tested_term$list=list_to_filter
return(tested_term)
}


# Pi ----------------------------------------------------------------------

genomic_scan_gene_pi=genomic_scan_gene%>%filter(analysis=="pi")
test=foreach(i=unique(genomic_scan_gene_pi$variable),.combine = rbind.data.frame)%do%
compute_enrichment(potentially_enriched_gene = genomic_scan_gene_pi,
                   genomic_annotation=genomic_annotation,
                   list_to_filter = i,
                   func_choice = "Product",
                   sep=";",
                   number_genes_picked=number_of_pi)

#Filter only Product inside
test_reduced=test[grep(x = test$Product,pattern = "lact|lipase|protease|lipid"),]
test_reduced=test_reduced%>%ungroup()%>%mutate(p_value_adjusted=pvalue*n() ,significant=pvalue<(0.05/n()))

write.table(test_reduced,"enriched_interest_term.csv",sep = "\t",row.names = F)

# dxy ---------------------------------------------------------------------
genomic_scan_gene_dxy=genomic_scan_gene%>%filter(analysis%in%"dxy"&variable%in%c("Cheese_1_Wild","Cheese_2_Wild","Cheese_3_Wild"))

number_of_dxy=number_of_dxy%>%group_by(variable)%>%summarise(n=sum(n))


test2=foreach(i=unique(genomic_scan_gene_dxy$variable),.combine = rbind.data.frame)%do%
  compute_enrichment(potentially_enriched_gene = genomic_scan_gene_dxy,
                     genomic_annotation=genomic_annotation,
                     list_to_filter = i,
                     func_choice = "Product",
                     sep=";",
                     number_genes_picked=number_of_dxy)

#Filter only Product inside
test_reduced2=test2[grep(x = test2$Product,pattern = "lact|lipase|protease|lipid"),]
test_reduced2=test_reduced2%>%ungroup()%>%mutate(p_value_adjusted=pvalue*n() ,significant=pvalue<(0.05/n()))


write.csv2(test_reduced2,"fisher_test_dxy.csv",row.names = F)
write.csv2(test_reduced,"fisher_test_pi.csv",row.names = F)

