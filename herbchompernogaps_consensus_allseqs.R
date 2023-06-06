### HerbChomper Beta 0.3 - for review - E. Gardner - June 4, 2020

## Modified by S. Musker 
## Aug 31, 2022 (v1)
## Dec 2022, (v2): Modified for better consensus-making and post-chomp filtering 
## Rather than specifying a reference and target sequence, the script calculates the majority rule consensus for the reference,
## and uses that to decide whether to chomp a bit of sequence.
## It does this for all sequences.
## alignments must be in fasta format

usage <- "Rscript herbchompernogaps_consensus_allseqs.R -a [alignment in] -o [alignment out] -w [size of sliding window] -i [identity cutoff]"
help_usage <- "
REQUIRED:
-a [alignment in (fasta)] 
-o [alignment out (fasta; excludes extension)] 
-w [size of sliding window in bp. Don't make this too low. Try 50.] 
-i [identity cutoff below which windows are recoded as gaps. Don't make this too high. Try 0.8.] 

OPTIONAL:
-c [default=1] [wide (1) or narrow (0) cut (for recoding windows to gaps)]

 Global fasta trimming options:
--trim_before [default=1] [Trim fasta before chomping? 1=yes, 0=no.]
--trim_after [default=1] [Trim fasta after chomping? 1=yes, 0=no.]
-m / --min_n_sites_per_seq [default=100] [Before AND after chomping, *sequences* with <= 'min_n_sites_per_seq' non-gap sites will be removed. If < 1 this is interpreted as a proportion.]
-p / --min_prop_seqs_per_site [default=0.1] [After chomping, *sites* <= 'min_prop_seqs_per_site' non-gap sites (proportion) will be removed]
 
 Consensus-building options:
--consensus_prop_gaps_allowed [default=0.9] [Sites with > 'consensus_prop_gaps_allowed' will be coded as missing in the consensus]
--consensus_minor_prop_ignore [default = 0.2] [Polymorphic sites with > 'consensus_minor_prop_ignore' minor allele frequency will be recoded as gaps in the consensus. Must be between 0 and 0.5]
--write_consensus [flag] [Write the consensus sequence as a separate fasta?]
--append_consensus [flag] [Append the consensus sequence to the output alignment?]
"
## get arguments
args <- commandArgs(trailingOnly=TRUE)
if ("-h" %in% args | "--help" %in% args | length(args) == 0) {
	cat("ARGUMENTS: \n-------",help_usage,"\n")
	cat("USAGE (minimal): \n-------\n",usage,"\n")
	q()
}
#alignment in
seqfile <- args[which(args == "-a")+1]
#output alignment
outfile <- args[which(args == "-o")+1]
#sliding window size
slidingwindow <- as.numeric(args[which(args == "-w")+1])
#identity cutoff
identity <- as.numeric(args[which(args == "-i")+1])

#wide or narrow cut
if ("-c" %in% args) {
	cutsize <- as.numeric(args[which(args == "-c")+1])
} else {
	cutsize <- 0
}

## Optional arguments ##
if ("--trim_before" %in% args) {
	trim_before <- as.numeric(args[which(args == "--trim_before")+1])
		if (!trim_before %in% c(0,1)) stop("--trim_before must be 0 or 1!")
} else {
	trim_before <- 1
}

if ("--trim_after" %in% args) {
	trim_after <- as.numeric(args[which(args == "--trim_after")+1])
		if (!trim_after %in% c(0,1)) stop("--trim_after must be 0 or 1!")
} else {
	trim_after <- 1
}

if ("--append_consensus" %in% args && (trim_after != 0)) {
	stop("Post-chomp trimming should be disabled if you want the consensus sequence appended to the output! Set --trim_after 0.")
}

#before and after chomping, sequences with less than than this many non-gap sites will be removed
if (any(grepl("-m",args) | grepl("--min_n_sites_per_seq",args))) {
	args[which(args == "--min_n_sites_per_seq")] <- "-m"
	MNPS <- as.numeric(args[which(args == "-m")+1])
} else {
	MNPS <- 100
}

# Sites with this many non-gap sites will be removed
if (any(grepl("-p",args) | grepl("--min_prop_seqs_per_site",args))) {
	args[which(args == "--min_prop_seqs_per_site")] <- "-p"
	MPPS <- as.numeric(args[which(args == "-p")+1])
} else {
	MPPS <- 0.1
}

## Options for consensus-making
if ("--consensus_prop_gaps_allowed" %in% args) {
	consensus_prop_gaps_allowed <- as.numeric(args[which(args == "--consensus_prop_gaps_allowed")+1])
} else {
	consensus_prop_gaps_allowed <- 0.9
}

if ("--consensus_minor_prop_ignore" %in% args) {
	consensus_minor_prop_ignore <- as.numeric(args[which(args == "--consensus_minor_prop_ignore")+1])
	if (consensus_minor_prop_ignore>0.5) stop("--consensus_minor_prop_ignore must be between 0 and 0.5!")
	} else {
	consensus_minor_prop_ignore <- 0.2
}

cat("Input file:",seqfile,"\n")
cat("Output file:",outfile,"\n")

## Done parsing arguments ##

###### To use in the R environment, uncomment this section, set these five variables and then run from here down
##working directory
# setwd("~")
# ##alignment
# seqfile <- "seq.fasta"
# outfile <- "seq.chompedConsensusAll.fasta"
# ##sliding window size
# slidingwindow <- 50
# #identity cutoff
# identity <- 0.8
# #size of gap required to restart trimming
# gapsize <- 10

#check dependency and load alignments
# if ("seqinr" %in% rownames(installed.packages()) == FALSE) install.packages("seqinr")
library(seqinr)
gene <- read.fasta(paste(seqfile))

## DEFINE functions ----
countNonGaps <- function(x) sum(!x %in% c("-","N","n"))
meetsNumNonGap <- function(x,num_nonGap) countNonGaps(x) > num_nonGap
meetsPropNonGap <- function(x,prop_nonGap) (countNonGaps(x)/length(x)) > prop_nonGap
getFastaWidth <- function(fa) unique(unlist(lapply(fa,length)))

if (length(getFastaWidth(gene)) > 1) stop("Fasta has varying sequence lengths -- is it aligned?")

trimFasta <- function(fa,min_n_sites_per_seq,min_prop_seqs_per_site) {
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

## smarter consensus
fa_as_aln <- function(x) as.alignment(length(x),names(x),lapply(x,c2s),NA)

consensus_ignoreGaps <- function(aln,prop_gaps_allowed,conflict_minor_prop_ignore) {
  temp <- seqinr::consensus(aln,method="profile") # produces a matrix of character counts per site
  rnames <- rownames(temp)
  gaprow <- which(rnames == "-")
  ngaprow <- which(rnames == "n")
  if ("-" %in% rnames || "n" %in% rnames) {
    temp_nonGap <- temp[-c(gaprow,ngaprow),]
  } else {
    temp_nonGap <- temp
  }
  rnames_nonGap <- rownames(temp_nonGap)
  ## 1. Get the majority sequence, ignoring gaps
  ##    note that for now ties are resolved arbitrarily, but later will be recoded as gaps
  ##    due to the minor proportion ratio (which is bounded at [0,0.5])
  majority_ignoreGaps <- unlist(apply(temp_nonGap,2,function(x) rnames_nonGap[which(x == max(x))][1]))
  ## 2. Find SNPs (to recode as gaps eventually)
  is_monomorphic <- apply(temp_nonGap,2,function(x) length(x[which(x!=0)]) == 1)
  conflict_ratio <- apply(temp_nonGap,2,function(x) min(x[x!=0])[1] / sum(x[x!=0]))
  passes_conflict_ratio <-  (conflict_ratio <= conflict_minor_prop_ignore)
  ## 3. Get index of sites with < "prop_gaps_allowed" proportion gaps
  prop_gaps <- apply(temp,2,function(x) (x[gaprow]+x[ngaprow]) / aln$nb)
  passes_prop_gaps <- prop_gaps <= prop_gaps_allowed
  cat("Generating majority-rule consensus sequence while requiring:\n\t(a). % of sequences with data:\t\t\t>",
      100*(1-prop_gaps_allowed),"%\n\tAND\n\t(b). If ambiguous, % conflicting sequences:\t<",
      conflict_minor_prop_ignore*100,"%\n")
  ngap <- sum(majority_ignoreGaps == "-")
  ## change sites not passing above filters to gaps
  majority_ignoreGaps[!(passes_prop_gaps)]  <- "-"
  ngap2 <- sum(majority_ignoreGaps == "-")
  cat("Results:\n\t(a). Characters in the consensus recoded as gaps by missingness:\t",ngap2-ngap,"\n")
  majority_ignoreGaps[!(passes_conflict_ratio) & !(is_monomorphic)] <- "-"
  ngap3 <- sum(majority_ignoreGaps == "-")
  cat("\tAND\n\t(b). Characters in the consensus recoded as gaps by ambiguity:\t",ngap3-ngap2,"\n")
  cat("\nFinal consensus has ",length(majority_ignoreGaps)-ngap3,
      " non-gaps out of a total of ",length(majority_ignoreGaps)," sites (",
      signif(ngap3*100/length(majority_ignoreGaps),2),"% gaps).\n",sep="")
  majority_ignoreGaps <- seqinr::as.SeqFastadna(majority_ignoreGaps)
  attr(majority_ignoreGaps,"Annot") <- paste0(">HerbComper_consensus_gaps_",prop_gaps_allowed,"_minor_conflict_",conflict_minor_prop_ignore)
  attr(majority_ignoreGaps,"name") <- paste0("HerbComper_consensus_gaps_",prop_gaps_allowed,"_minor_conflict_",conflict_minor_prop_ignore)
  majority_ignoreGaps <- list(majority_ignoreGaps)
  names(majority_ignoreGaps) <- paste0("HerbComper_consensus_gaps_",prop_gaps_allowed,"_minor_conflict_",conflict_minor_prop_ignore)
  return(majority_ignoreGaps) ## <-- a named list. Use [[1]] index to get the sequence itself.
}

calcIdentity <- function(sequence, reference, start, window) {
  s <-  sequence[start:(start+window)]
  r <- reference[start:(start+window)]
  s_is_r <- s == r
  s_is_missing <- s %in% c("-","n")
  r_is_missing <- r %in% c("-","n")
  s_and_r_are_not_missing <- !s_is_missing & !r_is_missing
  s_is_r_ignore_gap <- s_is_r & s_and_r_are_not_missing
  sequence_identity <- sum(s_is_r_ignore_gap)/sum(s_and_r_are_not_missing)
  if (is.na(sequence_identity)) sequence_identity <- 1 ## ignore empty windows. These will be trimmed by trimFasta()
  return(sequence_identity)
}

## Begin work! ----
if (trim_before==1) {
	cat("Performing initial fasta filtering to remove sequences with <",MNPS,"non-gap/non-N sites.\n")
	cat("Note that no sites (except completely empty sites) will be removed in this round!\n")
	gene <- trimFasta(gene,min_n_sites_per_seq=MNPS,min_prop_seqs_per_site=0)
} else {
	## just trim to remove all-gap sequences and sites
	cat("You have opted out of initial fasta filtering. Only empty sequences and sites will be trimmed.\n")
	gene <- trimFasta(gene,min_n_sites_per_seq=0,min_prop_seqs_per_site=0)
}

geneAlCon <- consensus_ignoreGaps(
	fa_as_aln(gene),
	prop_gaps_allowed=consensus_prop_gaps_allowed,
	conflict_minor_prop_ignore=consensus_minor_prop_ignore)

if ("--write_consensus" %in% args) {
	write.fasta(geneAlCon[[1]],
				names = paste0(outfile,"_HerbComper_consensus_gaps_",consensus_prop_gaps_allowed,
								"_minor_conflict_",consensus_minor_prop_ignore),
				file.out = paste0(outfile,
								".consensus.gaps_",consensus_prop_gaps_allowed,
								".minor_conflict_",consensus_minor_prop_ignore,
								".fasta"))
}

## begin Chomping
cat("\n**************************\nBegin chomping.\n")
w <- slidingwindow-1

for (target in 1:length(gene)) {
	#count the number of non-gap characters in the target sequence:
	length1 <- sum(gene[[target]][1:length(gene[[target]])]!="-")
	#find where to start - forward
	posF <- min(which(!gene[[target]] %in% c("-","n")))
	#find where to start - reverse
	posR <- max(which(!gene[[target]] %in% c("-","n")))

	## HERE starts the herbchompernogaps code  ----
	#scan forward
	# cutF <- 0
	# for (i in posF:min(posR,(length(gene[[target]])-slidingwindow))) {
	# 	if (sum(gene[[target]][i:(i+w)] == geneAlCon[i:(i+w)]) < (slidingwindow*identity)) {
	# 			cutF <- c(cutF,i)
	# 		} else {
	# 			next
	# 		}
	# }
	# #scan reverse
	# cutR <- -1
	# for (i in posR:max(posF,slidingwindow)) {
	# 	if (sum(gene[[target]][i:(i-w)] == geneAlCon[i:(i-w)]) < (slidingwindow*identity)) {
	# 			cutR <- c(cutR,i)
	# 		} else {
	# 			next
	# 		}
	# }

	#scan forward (While loop???)
	cutF <- 0
	for (i in seq(posF, posR, by = slidingwindow)) {
		if (calcIdentity(gene[[target]],geneAlCon[[1]],i,slidingwindow) < identity) {
			cutF <- c(cutF,c(i:(i+slidingwindow)))
		} else {
			next
		}
	}

	#scan reverse
	cutR <- -1
	for (i in seq(posR, posF, by = -1*slidingwindow)) {
		if (calcIdentity(gene[[target]],geneAlCon[[1]],i,slidingwindow) < identity) {
			cutR <- c(cutR,c((i+slidingwindow):i))
		} else {
			next
		}
	}
	# print(cutR);Sys.sleep(2)


	# cat("posF:",posF,"\n")
	# cat("cutF:",cutF,"\n")
	# cat("posR:",posR,"\n")
	# cat("cutR:",cutR,"\n")
	# cat("intersect(cutF,cutR): ",intersect(cutF,cutR),"\n")
	# cat("unique(c(cutF[2:length(cutF)],cutR[2:length(cutR)])):",unique(c(cutF[2:length(cutF)],cutR[2:length(cutR)])),"\n")
	#now cut the positions flagged in either the forward or reverse scans
	if (cutsize == 1) {
		## case when 0 sites trimmed
		if (length(cutF)>1 && length(cutR)>1) {
			gene[[target]][unique(c(cutF[2:length(cutF)],cutR[2:length(cutR)]))] <- "-"
		} else {
			gene[[target]][intersect(cutF,cutR)] <- "-"
		}
	} else {
		gene[[target]][intersect(cutF,cutR)] <- "-"
	}
	## HERE ends the herbchompernogaps code  ----
		length2 <- sum(gene[[target]][1:length(gene[[target]])]!="-")
		cat(paste(length1-length2," non-gap characters out of ",length1," removed for sequence ",names(gene[target]),"\n",sep="",collapse=""))
}

## remove any now-empty sequences and now-empty sites
cat("\n**************************\n
	Chomping finished! Next step is to remove any now-empty sequences and now-empty sites.\n")

if (trim_after==1) {
	gene <- trimFasta(gene,min_n_sites_per_seq=MNPS,min_prop_seqs_per_site=MPPS)
} else {
	## just trim to remove all-gap sequences and sites
	cat("You have opted out of fasta filtering after chomping. Only empty sequences and sites will be trimmed.\n")
	gene <- trimFasta(gene,min_n_sites_per_seq=0,min_prop_seqs_per_site=0)
}

if ("--append_consensus" %in% args) {
	gene <- append(gene,geneAlCon)
}

cat("Output fasta",outfile,"has",getFastaWidth(gene),"SITES.\n")
write.fasta(lapply(gene,toupper),names=names(gene),file=outfile)