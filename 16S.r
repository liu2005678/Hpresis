# finding SNPs of 16S
setwd("/media/logen/sd16b/Hp/16S/")
cleandatadir="/media/logen/sdc/h.pylori/2ndseq_cleandata/"
strains=dir(cleandatadir)
strainsr1=strains[grep("hp\\d+.clean.r1.fq.gz",strains)]
strainsr2=strains[grep("hp\\d+.clean.r2.fq.gz",strains)]
strainsr1;strainsr2
base=sub(".clean.r1.fq.gz","",strainsr1)
base
system("mkdir -p align")
system("mkdir -p pile2cns")

bowtie2_cmds=paste0("bowtie2 -p 24 -x index/16s -1 ", cleandatadir, strainsr1, " -2 ", cleandatadir, strainsr2, " -S align/",
                  base, "_16S.sam" )
samtools_flow2 =function(files_base,p=24,rm=TRUE){
  cmd6=c();cmd8=c()
  cmd6=paste0("samtools view  -@ ", p," -u ", files_base,".sam ",
              "| samtols sort -@ ",p, " -o ",files_base,"_sorted.bam -")
  cmd8=paste0("samtools index ",files_base, "_sorted.bam")
  cmd8=c(cmd6,cmd8)#,cmd9)
  return(cmd8)
}
samtools_cmds=samtools_flow2(paste0("align/",base,"_16S"))

pileup_cmds=paste0("samtools mpileup -f 16S.fa align/",base,"_16S_sorted.bam > pile2cns/",base,"_16S.mpileup") 
pile2snp_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2snp pile2cns/",base,"_16S.mpileup --min-coverage 100,  --min-reads2 10 >",
                     "pile2cns/",base,"_16S.snp")
pile2indel_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2indel pile2cns/",base,"_16S.mpileup --min-coverage 100,  --min-reads2 10 >",
                       "pile2cns/",base,"_16S.indel")
pile2cns_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2cns pile2cns/",base,"_16S.mpileup --min-coverage 100,  --min-reads2 10 >",
                     "pile2cns/",base,"_16S.cns")
cmds=c(bowtie2_cmds,samtools_cmds,pileup_cmds,pile2snp_cmds,pile2indel_cmds,pile2cns_cmds)
for(i in 1:length(cmds)){system(cmds[i])}

