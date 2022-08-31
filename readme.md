# HerbChomper_MSA
A bioinformatic tool for trimming poorly-aligned ends from DNA sequences in multiple sequence alignments.

This is a fork of E. Gardner's repo https://github.com/artocarpus/HerbChomper.
I've added a script (herbchomper_consensus_allseqs.R) that chomps all sequences in an aligned multifasta based on the majority-rule consensus.
It's basically just a loop added to the original code in herbchomper.R with a few other minor changes to fix errors that came up with indexing base positions.

See E. Gardner's fork for more details.
