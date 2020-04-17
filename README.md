# pyramus

pyramus is exploration code for processing covid-19 data. 

composition is a simple program to give some stats on the composition of a sequence file.

They both use the seqio package, which can read sequences
transparently from any compressed or uncompressed fasta, fastq,
SAM/BAM/CRAM or vgp-tools .seq file, or from a custom binary sequence
file.  To read SAM/BAM/CRAM you must have installed htslib, and
likewise for VGP seq you must have installed vgp-tools.  See Makefile
for options on including or excluding these.
