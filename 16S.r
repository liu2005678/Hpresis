# finding SNPs of 16S
setwd("/media/logen/sd16b/Hp/16S/")
cleandatadir="/media/logen/sdc/h.pylori/2ndseq_cleandata/"
strains=dir(cleandatadir)
strainsr1=strains[grep("hp\\d+.clean.r1.fq.gz",strains)]
strainsr2=strains[grep("hp\\d+.clean.r2.fq.gz",strains)]
strainsr1;strainsr2
base=sub(".clean.r1.fq.gz","",strainsr1)
base
bowtie2_map=paste0("bowtie2 -p 24 -x index/16s -1 ", cleandatadir, strainsr1, " -2 ", cleandatadir, strainsr2, " -S align/",
                  base, "_16S.sam" )
samtools_flow2 =function(files_base,p=24,rm=TRUE){
  cmd6=c();cmd8=c()
  cmd6=paste0("samtools view  -@ ", p," -u ", files_base,".sam ",
              "| samtols sort -@ ",p, " -o ",files_base,"_sorted.bam -")
  cmd8=paste0("samtools index ",files_base, "_sorted.bam")
  #cmd9=paste0("bamCoverage -bs 20 -p ",p, " -b ", files_base, "_sorted.bam ", " -o ", files_base, ".bw")
  #cmd9
  cmd8=c(cmd6,cmd8)#,cmd9)
  return(cmd8)
}
samtools_cmds=samtools_flow2(paste0("align/",base,"_16S"))



