# HerbChomper_MSA
A bioinformatic tool for trimming poorly-aligned ends from DNA sequences in multiple sequence alignments.

This is a fork of E. Gardner's repo https://github.com/artocarpus/HerbChomper.

### What's different?
I've added two scripts (`herbchomper_consensus_allseqs.R` and `herbchompernogaps_consensus_allseqs.R`) that chomp **all sequences** in an aligned multifasta based on the majority-rule consensus.
<!-- It's basically just a loop added to the original code in herbchomper.R with a few other minor changes to fix errors that came up with indexing base positions. -->

Updates:
1. Wrote my own consensus function (`consensus_ignoreGaps()`) that allows missingness and ambiguity (both tunable)
2. Added site and sequence filtering function (`trimFasta()`) based on missingness. This is performed both before and after chomping, with tunable parameters.

See E. Gardner's fork for more details.
