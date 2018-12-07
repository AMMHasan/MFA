# MFA
R script for performing ***Marker Frequency Analysis (MFA)*** from the depth of coverage, calculated from NGS data, to accompany the manuscript "Genomic Analysis of DNA Double-Strand Break Repair in *Escherichia coli*" (January 2018. Methods in Enzymology; DOI:10.1016/bs.mie.2018.09.001)

This R script comes as a part of an umbrella project - Establishing a better genomic DNA extraction method for MFA in the [Leach lab](http://leach.bio.ed.ac.uk/index.html), Institute of Cell Biology, University of Edinburgh.

For the convenience of users, certain nomenclature patterns for the genome and the depth of coverage files are followed:
- The reference genome file is given in .fasta format and the name of the file should start with “genome”, like “genome_ref.fasta”.
- The depth of coverage data is read from the files starting with "depth" and separated by "_" like "depth_DL4184Aara.txt". The last part, here, is indicating any experimentally relevant details without any gap followed by ".txt".
- Names of the depth of coverage files from stationary phase culture is similar to that of logarithmic-phase data, plus, starting with "stn", like “stn_depth_DL4184Aara.txt".

### IMPORTANT: if you are missing the depth of coverage data from stationary phase culture, comment out from line 75 to line 120 in the "MFA_V2.1.R" R script.

