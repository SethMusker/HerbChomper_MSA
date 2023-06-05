### HerbChomper Beta 0.3 - for review - E. Gardner - June 4, 2020

## Modified by S. Musker 
## Aug 31, 2022 (v1)
## Dec 2022, (v2): Modified for better consensus-making and post-chomp filtering 
## Rather than specifying a reference and target sequence, the script calculates the majority rule consensus for the reference,
## and uses that to decide whether to chomp a bit of sequence.
## It does this for all sequences.
## alignments must be in fasta format

### Usage: Rscript herbchomper_consensus_allseqs.R -a [alignment in] -o [alignment out] -w [size of sliding window] -i [identity cutoff] -g [gap size necessary to restart trimming]
usage<-"Rscript herbchomper_consensus_allseqs.R -a [alignment in] -o [alignment out] -w [size of sliding window] -i [identity cutoff] -g [gap size necessary to restart trimming]"
help_usage <- "
REQUIRED:
-a [alignment in (fasta)] 
-o [alignment out (fasta; excludes extension)] 
-w [size of sliding window in bp. Don't make this too low. Try 50.] 
-i [identity cutoff below which windows are recoded as gaps. Don't make this too high. Try 0.8.] 
-g [gap size necessary to restart trimming; gaps below this size that are encountered while scanning right or left will cause chomping to terminate. Use herbchompernogaps_consensus_allseqs.R if you want to avoid this.]

OPTIONAL:
-c [default=1] [wide (1) or narrow (0) cut (for recoding windows to gaps)]

 Global fasta trimming options:
--min_n_sites_per_seq [default=100] [Before and after chomping, *sequences* with <= 'min_n_sites_per_seq' non-gap sites will be removed. If < 1 this is interpreted as a proportion.]
--min_prop_seqs_per_site [default=0.1] [Before and after chomping, *sites* <= 'min_prop_seqs_per_site' non-gap sites (proportion) will be removed]
 
 Consensus-building options:
--consensus_prop_gaps_allowed [default=0.9] [Sites with > 'consensus_prop_gaps_allowed' will be coded as missing in the consensus]
--consensus_minor_prop_ignore [default = 0.2] [Polymorphic sites with > 'consensus_minor_prop_ignore' minor allele frequency will be recoded as gaps in the consensus. Must be between 0 and 0.5]
"
## get arguments
args = commandArgs(trailingOnly=TRUE)
if ("-h" %in% args | "--help" %in% args | length(args) == 0) {
	cat("ARGUMENTS: \n-------",help_usage,"\n")
	cat("USAGE (minimal): \n-------\n",usage,"\n")
	q()
}
#alignment in
	seqfile<-args[which(args=="-a")+1]
	#output alignment
	outfile<-args[which(args=="-o")+1]
	#sliding window size
	slidingwindow<-as.numeric(args[which(args=="-w")+1])
	#identity cutoff
	identity<-as.numeric(args[which(args=="-i")+1])
	#size of gap required to restart trimming
	gapsize<-as.numeric(args[which(args=="-g")+1])

	## Optional arguments ##
	#after chomping, sequences with less than than this many non-gap sites will be removed
	if(any(grepl("-m",args) | grepl("--min_n_sites_per_seq",args))){
	  args[which(args=="--min_n_sites_per_seq")]<-"-m"
	  MNPS<-as.numeric(args[which(args=="-m")+1])
	}else{
	  MNPS<-10
	}

	# Sites with this many non-gap sites will be removed
	if(any(grepl("-p",args) | grepl("--min_prop_seqs_per_site",args))){
	  args[which(args=="--min_prop_seqs_per_site")]<-"-p"
	  MPPS<-as.numeric(args[which(args=="-p")+1])
	}else{
	  MPPS<-0.05
	}

	## Options for consensus-making
	if("--consensus_prop_gaps_allowed" %in% args){
	  consensus_prop_gaps_allowed<-as.numeric(args[which(args=="--consensus_prop_gaps_allowed")+1])
	}else{
	  consensus_prop_gaps_allowed<-0.9
	}

	if("--consensus_minor_prop_ignore" %in% args){
	  consensus_minor_prop_ignore<-as.numeric(args[which(args=="--consensus_minor_prop_ignore")+1])
	  if (consensus_minor_prop_ignore>0.5) stop("Consensus_minor_prop_ignore must be between 0 and 0.5!")
	  }else{
	  consensus_minor_prop_ignore<-0.2
	}

	cat("Input file:",seqfile,"\n")
	cat("Output file:",outfile,"\n")

	## Done parsing arguments ##

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
	gene<-read.fasta(paste(seqfile))

	countNonGaps<-function(x) sum(!x %in% c("-","N","n"))
	exceedsNumNonGap<-function(x,num_nonGap) countNonGaps(x)>num_nonGap
	exceedsPropNonGap<-function(x,prop_nonGap) (countNonGaps(x)/length(x))>prop_nonGap
	getFastaWidth<-function(fa) unique(unlist(lapply(fa,length)))

	trimFasta <- function(fa,min_n_sites_per_seq=10,min_prop_seqs_per_site=0.05) {
	# iteratively trims until there are
	## 1. no empty sequences AND
	## 2. no effectively empty sites
	fa_width <- getFastaWidth(fa)
	condition <- TRUE
	iter <- 0
	while (condition == TRUE) {
	iter <- iter+1
	cat("Trimming iteration\t",iter,"\n")
		# Find and remove gappy sites
		is_empty_site <- !apply(do.call(rbind,fa),2,meetsPropNonGap,prop_nonGap=min_prop_seqs_per_site)
		if (any(is_empty_site)) {
			cat(sum(is_empty_site),"out of",fa_width,"SITES have <",min_prop_seqs_per_site*100,"% sequences with data (not N, not gap). Removing them.\n")
			fa <- lapply(fa,function(x) x[-which(is_empty_site)])
			found_missing_sites<-TRUE
		} else {
			cat("No missing sites! Done trimming.\n")
			found_missing_sites<-FALSE
		}
		# Find and remove gappy sequences
		if (min_n_sites_per_seq>1) {
			is_empty_seq <- !unlist(lapply(fa,meetsNumNonGap,num_nonGap=min_n_sites_per_seq))
		} else {
			is_empty_seq <- !unlist(lapply(fa,meetsPropNonGap,prop_nonGap=min_n_sites_per_seq))
		}
		if (any(is_empty_seq)) {
			cat(sum(is_empty_seq),"out of",length(fa),"SEQUENCES have <",min_n_sites_per_seq,"sites with data (not N, not gap). Removing them.\n")
			cat("Sequences removed:",names(fa)[is_empty_seq],sep = "\n")
			found_missing_seqs<-TRUE
			fa <- fa[!is_empty_seq]
		} else {
			cat("No empty sequences found in this iteration.\n")
			found_missing_seqs<-FALSE
		}
		if (!found_missing_sites && !found_missing_seqs) condition <- FALSE
	}
	return(fa)
	}	

	##
	cat("Performing initial fasta filtering to remove sequences with <",MNPS,"non-gap/non-N sites.\n")
	gene<-trimFasta(gene,min_n_sites_per_seq=MNPS,min_prop_seqs_per_site=0)

	## smarter consensus

	fa_as_aln<-function(x) as.alignment(length(x),names(x),lapply(x,c2s),NA)

	consensus_ignoreGaps<-function(aln,prop_gaps_allowed,conflict_minor_prop_ignore) {
	  temp<-seqinr::consensus(aln,method="profile")
	  rnames<-rownames(temp)
	  gaprow<-which(rnames=="-")
	  ngaprow<-which(rnames=="n")
	  temp_nonGap<-temp[-c(gaprow,ngaprow),]
	  rnames_nonGap<-rownames(temp_nonGap)
	  ## 1. Get the majority sequence, ignoring gaps
	  ##    note that for now ties are resolved arbitrarily, but later will be recoded as gaps 
	  ##     due to the minor proportion ratio (which is bounded at [0,0.5])
	  majority_ignoreGaps<-unlist(apply(temp_nonGap,2,function(x)rnames_nonGap[which(x==max(x))][1]))
	  ## 2. Find SNPs (to recode as gaps eventually)
	  is_monomorphic<-apply(temp_nonGap,2,function(x)length(x[which(x!=0)])==1)
	  conflict_ratio<-apply(temp_nonGap,2,function(x)min(x[x!=0])[1]/sum(x[x!=0]))
	  passes_conflict_ratio<- (conflict_ratio <= conflict_minor_prop_ignore) 
	  ## 3. Get index of sites with < "prop_gaps_allowed" proportion gaps
	  prop_gaps<-apply(temp,2,function(x) (x[gaprow]+x[ngaprow])/aln$nb )
	  passes_prop_gaps<-prop_gaps <= prop_gaps_allowed
	  cat("Generating majority-rule consensus sequence while requiring:\n\ta). % of sequences with data:\t\t\t>",
		  100*(1-prop_gaps_allowed),"%\n\tAND\n\tb). If ambiguous, % conflicting sequences:\t<",
		  conflict_minor_prop_ignore*100,"%\n")
	  ngap<-sum(majority_ignoreGaps=="-")
	  ## change sites not passing above filters to gaps
	  majority_ignoreGaps[!(passes_prop_gaps)] <-"-"
	  ngap2<-sum(majority_ignoreGaps=="-")
	  cat("Results:\n\ta). Characters in the consensus recoded as gaps by missingness:\t",ngap2-ngap,"\n")
	  majority_ignoreGaps[!(passes_conflict_ratio) & !(is_monomorphic)]<-"-"
	  ngap3<-sum(majority_ignoreGaps=="-")
	  cat("\tAND\n\tb). Characters in the consensus recoded as gaps by ambiguity:\t",ngap3-ngap2,"\n")
	  cat("\nFinal consensus has ",length(majority_ignoreGaps)-ngap3,
		  " non-gaps out of a total of ",length(majority_ignoreGaps)," sites (",
		  signif(ngap3*100/length(majority_ignoreGaps),2),"% gaps).\n",sep="")
	  majority_ignoreGaps
	}

	geneAlCon<-consensus_ignoreGaps(fa_as_aln(gene),
									prop_gaps_allowed=consensus_prop_gaps_allowed,
									conflict_minor_prop_ignore=consensus_minor_prop_ignore)
	if("--write_consensus" %in% args){
	  write.fasta(geneAlCon,
				  names = paste0("HerbComper_consensus_gaps_",consensus_prop_gaps_allowed,
								 "_minor_conflict_",consensus_minor_prop_ignore),
				  file.out = paste0(outfile,
									".consensus.gaps_",consensus_prop_gaps_allowed,
									".minor_conflict_",consensus_minor_prop_ignore,
									".fasta"))
	}

	## begin Chomping
	cat("\n**************************\nBegin chomping.\n")
	w<-slidingwindow-1
	g<-gapsize-1

	for (target in 1:length(gene)) {
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
	  #trim forward ends of gaps
	  for (i in posF:posR) {
		if (sum(gene[[target]][i:(i+g)]=="-")==gapsize) {
		  gapF<-i+g	## position to start trimming after gap
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
		  gapR<-i-g	## position to start trimming after gap
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

	## remove any now-empty sequences and now-empty sites
	cat("\n**************************\n
		Chomping finished! Next step is to remove any now-empty sequences and now-empty sites.\n")

	gene<-trimFasta(gene,min_n_sites_per_seq=MNPS,min_prop_seqs_per_site=MPPS)

	cat("Output fasta",outfile,"has",getFastaWidth(gene),"SITES.\n")

	write.fasta(lapply(gene,toupper),names=names(gene),file=outfile)