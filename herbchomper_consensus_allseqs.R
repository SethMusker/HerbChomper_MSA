### HerbChomper Beta 0.3 - for review - E. Gardner - June 4, 2020

## Modified by S. Musker Aug 31, 2022
## Rather than specifying a reference and target sequence, the script calculates the majority rule consensus for the reference,
## and uses that to decide whether to chomp a bit of sequence.
## It does this for all sequences.
## alignments must be in fasta format

### Usage: Rscript herbchomper_consensus_allseqs.R -a [alignment in] -o [alignment out] -w [size of sliding window] -i [identity cutoff] -g [gap size necessary to restart trimming]

#get arguments
args = commandArgs(trailingOnly=TRUE)
#alignment in
seqfile<-args[grep("-a",args)+1]
#output alignment
outfile<-args[grep("-o",args)+1]
#sliding window size
slidingwindow<-as.numeric(args[grep("-w",args)+1])
#identity cutoff
identity<-as.numeric(args[grep("-i",args)+1])
#size of gap required to restart trimming
gapsize<-as.numeric(args[grep("-g",args)+1])

###### To use in the R environment, uncomment this section, set these five variables and then run from here down
##working directory
# setwd("~")
# ##alignment
# seqfile<-"seq.fasta"
# outfile<-"seq.chompedConsensusAll.fasta"
# ##sliding window size
# slidingwindow<-50
# #identity cutoff
# identity<-0.8
# #size of gap required to restart trimming
# gapsize<-10

#check dependency and load alignments
if ("seqinr" %in% rownames(installed.packages()) == FALSE) {install.packages("seqinr")}
library(seqinr)
read.fasta(paste(seqfile))->gene
w<-slidingwindow-1
g<-gapsize-1

geneAlCon<-consensus(read.alignment(paste(seqfile),format ="fasta"),method = "majority")

print(paste("Input file:",seqfile,collapse=" "))

for(target in 1:length(gene)){
  #count the number of non-gap characters in the target sequence:
  length1<-sum(gene[[target]][1:length(gene[[target]])]!="-")
  
  #find where to start - forward
  for (i in 1:length(gene[[target]])) {
    if (gene[[target]][i]=="-") {
      next
    }
    else {
      posF<-i
      break
    }
  }
  
  #find where to start - reverse
  for (i in length(gene[[target]]):1) {
    if (gene[[target]][i]=="-") {
      next
    }
    else {
      posR<-i
      break
    }
  }
  
  
  #trim forward end
  for (i in posF:posR) {
    if (sum(gene[[target]][i:(i+w)]==geneAlCon[i:(i+w)]) < (slidingwindow*identity)) {
      gene[[target]][i]<-"-"
    }
    else {
      break
    }
  }
  
  #trim reverse
  for (i in posR:posF) {
    if (sum(gene[[target]][i:(i-w)]==geneAlCon[i:(i-w)]) < (slidingwindow*identity)) {
      gene[[target]][i]<-"-"
    }
    else {
      break
    }
  }
  
  #fix up starting and ending positions so that we don't get errors by running up against the ends
  posF<-max(posF,gapsize)
  posR<-min(posR,length(gene[[target]])-gapsize-1)
  
  # print("CHOMPING... yum yum yum")
  
  #trim forward ends of gaps
  for (i in posF:posR) {
    if (sum(gene[[target]][i:(i+g)]=="-")==gapsize) {
      i+g->gapF	## position to start trimming after gap
      gapF<-min(gapF,posR)
      for (j in gapF:posR) {
        if((j+w) < length(geneAlCon)){
          goto<-(j+w)
        }else{
          goto<-length(geneAlCon)
        }
        if (sum(gene[[target]][j:goto]==geneAlCon[j:goto]) < (slidingwindow*identity)) {
          gene[[target]][j]<-"-"
        }
        else {
          break
        }
      }
    }
    else {
      next
    }
  }
  
  
  #trim reverse ends of gaps
  for (i in posR:posF) {
    if (sum(gene[[target]][i:(i-g)]=="-")==gapsize) {
      i-g->gapR	## position to start trimming after gap
      gapR<-max(gapR,posF)
      for (j in gapR:posF) {
        if(j>w){
          goto<-(j-w)
        }else{
          goto<-1
        }
        if (sum(gene[[target]][j:goto]==geneAlCon[j:goto]) < (slidingwindow*identity)) {
          gene[[target]][j]<-"-"
        }
        else {
          break
        }
      }
    }
    else {
      next
    }
  }
  
  length2<-sum(gene[[target]][1:length(gene[[target]])]!="-")
  cat(paste(length1-length2," non-gap characters out of ",length1," removed for sequence ",names(gene[target]),"\n",sep="",collapse=""))
  
}

write.fasta(gene,names=names(gene),file=paste(outfile))
