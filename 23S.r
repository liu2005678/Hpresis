# finding SNPs of 23S and 5S
setwd("/media/logen/sd16b/Hp/23S/")  # set work directory
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
system("bowtie2-build -f 23S_5S.fa index/23S")
# align reads to reference 23S sequence, keep aligned reads only (--no-unal) 
# the CPU threads used is 24 (-p 24) 
# to allow better align of reads that mapped to terminals, 100 bp upstream of 23S and 100 bp downstream of 5S
bowtie2_cmds=paste0("bowtie2 -p 24 --no-unal -x index/23S -1 ", cleandatadir, strainsr1, " -2 ", cleandatadir, strainsr2, " -S align/",
                    base, "_23S.sam" )
# sort the aligned reads
samtools_flow2 =function(files_base,p=24,rm=TRUE){
  cmd6=c();cmd8=c()
  cmd6=paste0("samtools view  -@ ", p," -u ", files_base,".sam ",
              "| samtools sort -@ ",p, " -o ",files_base,"_sorted.bam -")
  cmd8=paste0("samtools index ",files_base, "_sorted.bam")
  cmd8=c(cmd6,cmd8)#,cmd9)
  return(cmd8)
}
samtools_cmds=samtools_flow2(paste0("align/",base,"_23S"))

# make pileup file
pileup_cmds=paste0("samtools mpileup -f 23S_5S.fa align/",base,"_23S_sorted.bam > pile2cns/",base,"_23S.mpileup") 
# snp and indel calling  A minimal coverage of 100 and variant supporting reads of more than 10 was set (--min-coverage 100  --min-reads2 10)
pile2snp_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2snp pile2cns/",base,"_23S.mpileup --min-coverage 100  --min-reads2 10 >",
                     "pile2cns/",base,"_16S.snp")
pile2indel_cmds=paste0("java -jar /home/logen/programs/VarScan.v2.3.9.jar pileup2indel pile2cns/",base,"_23S.mpileup --min-coverage 100  --min-reads2 10 >",
                       "pile2cns/",base,"_23S.indel")
cmds=c(bowtie2_cmds,samtools_cmds,pileup_cmds,pile2snp_cmds,pile2indel_cmds)
for(i in 1:length(cmds)){system(cmds[i])}

