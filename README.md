<p align="center" >
    <img src="https://github.com/SAMtoBAM/TLHcrs_finder/blob/main/logo/TLHcrsFinder_logo.png" width=100%>
</p>


# TLHcrsFinder
TLHcrsFinder is a tools desgined to detect and analyse **T**elomere-**L**inked-**H**elicase **C**ontaining Region**S** (TLHcrs) repeats

TLHcrs repeats are, as the name describes, regions adjacent to telomeres that are conserved across chromosome ends and generally contain helicase genes <br/>
These repeats, such as the most perhaps most well known Y-prime region in _S. cerevisiae_, are conserved across diverse fungi <br/>
With the advent of Long-read sequencing we now have an increasing number of well assembled genomes with the subtelomeres intact <br/>
Therefore we are now primed for looking at the TLHcrs repeats and their evolution.

This tool was developed for XXXX (please cite this publication if you find this tool useful)

# Conda installation

    conda install samtobam::TLHcrsfinder

# How to use

    TLHcrsFinder.sh -a assembly.fa




## How does TLHcrsFinder work

TLHcrsFinder works in 10 main steps

1. Extract 50kb from the ends of contigs >100kb
2. Identify Telomeric sequences genomes wide and determine the number within the contig ends
3. Align all contig ends to one another
4. Calculate the average coverage of sliding windows within the 50kb ends
5. Extract contiguous regions with a coverage > (0.75*'the number of telomeric sequences')
6. Remove regions < 2kb
7. Cluster the nucleotide sequence of the remaining regions using the 80/80 principal
8. Consider the largest cluster as the TLHcrs repeat and take the largest version as the representative
9. Search the whole genome using BLASTn using the representative TLHcrs repeat
10. Plot the alignment and whole genome positions of TLHcrs repeats for manual verification/scrutiny



