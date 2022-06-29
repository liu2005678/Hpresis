# finding SNPs of 16S
setwd("/media/logen/sd16b/Hp/16S/")  # set work directory
cleandatadir="/media/logen/sdc/h.pylori/2ndseq_cleandata/"   #directory for clean illumina data 
strains=dir(cleandatadir)
strainsr1=strains[grep("hp\\d+.clean.r1.fq.gz",strains)]
strainsr2=strains[grep("hp\\d+.clean.r2.fq.gz",strains)]
strainsr1;strainsr2
base=sub(".clean.r1.fq.gz","",strainsr1)
base
system("mkdir -p align")
system("mkdir -p pile2cns")
system("mkdir -p index")

# to allow better aligning of reads that mapped to the two terminals of 16 S, 100 bp surrounding sequence (upstream 100 bp and dowstream 100 bp) were included in the 16S.fa
system("bowtie2-build -f 16S.fa index/16s")

# align reads to reference 23S sequence, keep aligned reads only (--no-unal) 
# the CPU threads used is 24 (-p 24) 
bowtie2_cmds=paste0("bowtie2 -p 24 --no-unal -x index/16s -1 ", cleandatadir, strainsr1, " -2 ", cleandatadir, strainsr2, " -S align/",
                  base, "_16S.sam" )
# sort the aligned reads
samtools_flow2 =function(files_base,p=24,rm=TRUE){
  cmd6=c();cmd8=c()
  cmd6=paste0("samtools view  -@ ", p," -u ", files_base,".sam ",
              "| samtols sort -@ ",p, " -o ",files_base,"_sorted.bam -")
  cmd8=paste0("samtools index ",files_base, "_sorted.bam")
  cmd8=c(cmd6,cmd8)#,cmd9)
  return(cmd8)
}
samtools_cmds=samtools_flow2(paste0("align/",base,"_16S"))

# make pileup file
pileup_cmds=paste0("samtools mpileup -f 16S.fa align/",base,"_16S_sorted.bam > pile2cns/",base,"_16S.mpileup") 
# snp and indel calling  A minimal coverage of 100 and variant supporting reads of more than 10 was set (--min-coverage 100,  --min-reads2 10)
pile2snp_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2snp pile2cns/",base,"_16S.mpileup --min-coverage 100,  --min-reads2 10 >",
                     "pile2cns/",base,"_16S.snp")
pile2indel_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2indel pile2cns/",base,"_16S.mpileup --min-coverage 100,  --min-reads2 10 >",
                       "pile2cns/",base,"_16S.indel")
pile2cns_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2cns pile2cns/",base,"_16S.mpileup --min-coverage 100,  --min-reads2 10 >",
                     "pile2cns/",base,"_16S.cns")
cmds=c(bowtie2_cmds,samtools_cmds,pileup_cmds,pile2snp_cmds,pile2indel_cmds,pile2cns_cmds)
# excute commands generated
for(i in 1:length(cmds)){system(cmds[i])}

