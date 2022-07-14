library(tidyverse)
library(RColorBrewer)
setwd("/media/logen/sdc/h.pylori/snps_nanopore")
strains=dir("/media/logen/sdc/h.pylori/snps_nanopore/snp16S")
snp16S=data.frame()
snp23S=data.frame()

# read in phased long SNVs (including SNPs and indels)
for(i in 1:length(strains)){
  tmp=read.delim(paste0("snp16S/",strains[i],"/phased_merge_output.vcf.gz"),header=F,comment.char = "#")
  tmp$strain=strains[i]
  tmp2=read.delim(paste0("snp23S/",strains[i],"/phased_merge_output.vcf.gz"),header=F,comment.char = "#")
  tmp2$strain=strains[i]
  snp16S=rbind(snp16S,tmp)
  snp23S=rbind(snp23S,tmp2)
}
snpslong=rbind(snp16S,snp23S)
colnames(snpslong)=c("CHROM" , "POS" ,   "ID"   ,  "REF"   , "ALT"   , "QUAL"  , "FILTER", "INFO" ,  "FORMAT", "SAMPLE","strain")
snpslong=snpslong[which(snpslong$FILTER =="PASS"),]

# The strains ID may be different in long reads and short reads, i.e. Hpfe0001, hp0001; change (ids of short reads) to it.
samples=read.csv("/media/logen/sdc/h.pylori/samples.csv",header=T)
head(samples)
samples_sub=samples[which(samples$strain.1 %in% strains),]
sam_il=unique(samples_sub$strain)
sam_il
samples_sub= samples_sub %>% distinct(strain.1, .keep_all = TRUE)

# read in snps and indels of short reads respectively
setwd("/media/logen/sdc/h.pylori/snps_illumina/")
snpdf=data.frame();indeldf=data.frame()
for(i in 1:nrow(samples_sub)){
  df=read.delim(paste0("16S/pile2cns/",samples_sub$strain[i],"_16S.snp"),header = T)
  df2=read.delim(paste0("23S/pile2cns/",samples_sub$strain[i],"_23S.snp"),header = T)
  df$strain=samples_sub$strain.1[i]
  df2$strain=samples_sub$strain.1[i]
  snpdf=rbind(snpdf,df,df2)
  
  df=read.delim(paste0("16S/pile2cns/",samples_sub$strain[i],"_16S.indel"),header = T)
  df$strain=samples_sub$strain.1[i]
  df2=read.delim(paste0("23S/pile2cns/",samples_sub$strain[i],"_23S.indel"),header = T)
  df2$strain=samples_sub$strain.1[i]
  indeldf=rbind(indeldf,df,df2)
}

# filter SNVs foud by short reads with Variance frequency of more than 5%
snpshort=rbind(snpdf, indeldf)
snpshort$VarFreq=as.numeric( sub("%","",snpshort$VarFreq))
snpshort=snpshort[which(snpshort$VarFreq > 5),]

# reference of strain 26695 site 702 755 756 789 1320 are note sure, discard SNVs in this sites
snpshort16S=snpshort[which(snpshort$Chrom == "16S" & snpshort$Position > 100 & snpshort$Position <1602),]
snpshort16S=snpshort16S[which(! snpshort16S$Position %in% c(702,755,756,789,1320)),]
snpshort23S=snpshort[which(snpshort$Chrom == "23S_5S" & snpshort$Position > 100 & snpshort$Position < 3431),]
snpslong16S=snpslong[which(snpslong$CHROM== "16S" & snpslong$POS > 100 & snpslong$POS < 1602),]
snpslong16S=snpslong16S[which(! snpslong16S$POS %in% c(702,755,756,789,1320)),]
snpslong23S=snpslong[which(snpslong$CHROM== "23S_5S" & snpslong$POS > 100 & snpslong$POS < 3431),]
snpshort=rbind(snpshort16S,snpshort23S)
snpslong=rbind(snpslong16S,snpslong23S)
setwd("/media/logen/sdc/h.pylori/snp_compare/")
write.table(snpshort16S,file="snpshort16S.tab",sep = "\t",row.names = F,col.names = T)
write.table(snpshort23S,file="snpshort23S.tab",sep = "\t",row.names = F,col.names = T)
write.table(snpshort,file="snpshort.tab",sep = "\t",row.names = F,col.names = T)
write.table(snpslong,file="snpslong.tab",sep="\t",row.names = F,col.names = T)
write.table(snpslong16S,file="snpslong16S.tab",sep="\t",row.names = F,col.names = T)
write.table(snpslong23S,file="snpslong23S.tab",sep="\t",row.names = F,col.names = T)

########################################################################################################
# draw picture of comparison of short and long of 16S
head(snpslong16S)
head(snpshort16S)
snpshort16S$method="illumina"
snpslong16S$method="nanopore"
colsneed=c("chr","pos","strain","method")
colnames(snpshort16S)[c(1,2,20,21)]=colsneed
colnames(snpslong16S)[c(1,2,11,12)]=colsneed

# combine SNVs of short and long to a single dataframe
mix=data.frame()
strains= unique(snpslong$strain)
for(i in 1:length(strains)){
  sublong=snpslong16S[which(snpslong16S$strain == strains[i]),colsneed]
  subshort=snpshort16S[which(snpshort16S$strain == strains[i]),colsneed]
  submix=rbind(sublong,subshort)
  submix$dup=ave(submix$pos,submix$pos,FUN = length) > 1L
  submix$com="0"
  for(j in 1:nrow(submix)){if(submix$dup[j]==TRUE){submix$com[j]="Common"} else{submix$com[j]=submix$method[j]}}
  submix=submix[!duplicated(submix$pos),]
  mix=rbind(mix,submix)
}

# exlude two strains that are indicated to be contanminated
mix=mix[which(!mix$strain %in% c("Hpfe091","Hpfe110")),]
write.table(mix[,c(1,2,3,6)],file="snps_combine_16S.tab",quote=F,row.names = F,col.names = T)

# group strains
by_strain = mix %>% group_by(strain,com)  
by_strain2= by_strain %>% summarise(n=n())
# make a summary of snps found by long reads
by_strain_long=by_strain2[which(by_strain2$com %in% c("Common","nanopore")),]
sum_long=by_strain_long %>% summarise(sum=sum(n))
summary(sum_long)
# make a summary of snps found by short reads
by_strain_short=by_strain2[which(by_strain2$com %in% c("Common","illumina")),]
sum_short=by_strain_short %>% summarise((sum=sum(n)))
summary(sum_short)

# sort strains by the number of total SNPs found
by_strain3= by_strain2 %>% summarise(sum=sum(n))
by_strain3=by_strain3[order(by_strain3$sum),]
by_strain3$order=as.integer(rownames(by_strain3))
by_strain2$order=0
for(i in 1:nrow(by_strain2)){
 by_strain2$order[i]=as.integer(by_strain3$order[grepl(by_strain2$strain[i],by_strain3$strain)])
}
# arrange the order of snps that shown in the plot
by_strain2$com=factor(by_strain2$com,levels = c("illumina","nanopore","Common"))
df=data.frame(by_strain2)
head(df)

ggplot(df,aes(x=order,y=n,fill=com))+ geom_bar(stat="identity",position = "stack",alpha=0.5)+
  #scale_fill_manual(values = alpha(c(brewer.pal(6,"Set2")[c(1,2,3)]), .8))+
  #geom_line(position="stack",size=.25,color="black")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour="grey",fill = NA), strip.background = element_blank())

############################################################################################################################
# draw plot of 23S, stack_area plots
# extract data that are necessary for plot
nanopore_sub=data.frame(sites=snpslong23S$POS,strains=snpslong23S$strain)
head(nanopore_sub)
illumina_sub=data.frame(sites=snpshort23S$Position,strains=snpshort23S$strain)
head(illumina_sub)
strains=sort(unique(snpshort23S$strain))
strains=strains[which(! strains %in% c("Hpfe091","Hpfe110"))]

# calculate SNPs found by three methods
df=data.frame(strain=strains,comman=0,illu=0,nano=0,total=0)
df=df[which(! df$strain %in% c("Hpfe091","Hpfe110")),]
for(i in 1:length(strains)){
  nanopore_strain=nanopore_sub[which(nanopore_sub$strains==strains[i]),]
  illumina_strain=illumina_sub[which(illumina_sub$strains==strains[i]),]
  sites=unique(c(nanopore_strain$sites,illumina_strain$sites))
  common=nrow(nanopore_strain[which(nanopore_strain$sites %in% illumina_strain$sites),])
  nano=nrow(nanopore_strain[which(! nanopore_strain$sites %in% illumina_strain$sites),])
  illu=nrow(illumina_strain[which(! illumina_strain$sites %in% nanopore_strain$sites),])
  total=common + illu + nano
  df[i,2]=common
  df[i,3]=illu
  df[i,4]=nano
  df[i,5]=total
}
df=df[order(df$total),]
# summary of SNPs
comman_mut=sum(df$comman)
illu_mut=sum(df$illu)
nano_mut=sum(df$nano)
total_mut=sum(df$total)
df2=df
df2$strains=seq(1,98)
df2=df2[,c(2,3,4,6)]
colnames(df2)=c("common","illumina_only","nanopore_only","Hp_strains")
df_melt=reshape2::melt(df2,id.vars="Hp_strains",variable.name="sort",value.name="SNVs")
head(df_melt)
df_melt$sort=factor(df_melt$sort,levels = c("illumina_only","nanopore_only","common"))
ggplot(df_melt,aes(x=Hp_strains,y=SNVs,fill=sort))+  geom_area(position="stack")+
  scale_fill_manual(values = alpha(c(brewer.pal(6,"Set2")[c(1,2,3)]), .8))+
  #geom_line(position = "stack",size=0.5,color="grey",alpha=0.5)+
  theme(panel.background = element_blank(),panel.border = element_rect(colour="grey",fill = NA), strip.background = element_blank())

