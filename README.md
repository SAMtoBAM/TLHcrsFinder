<p align="center" >
    <img src="https://github.com/SAMtoBAM/TLHcrs_finder/blob/main/logo/TLHcrsFinder_logo.png" width=100%>
</p>


# TLHcrsFinder
TLHcrsFinder is a tools desgined to detect and analyse **T**elomere-**L**inked-**H**elicase **C**ontaining Region**S** (TLHcrs) repeats

TLHcrs repeats are, as the name describes, regions adjacent to telomeres that are conserved across chromosome ends and generally contain helicase genes <br/>
These repeats, such as the most perhaps most well known Y-prime region in _S. cerevisiae_, are conserved across diverse fungi <br/>
With the advent of Long-read sequencing we now have an increasing number of well assembled genomes with the subtelomeres intact <br/>
Therefore we are now primed for looking at the TLHcrs repeats and their evolution.

All that is required to run TLHcrsFinder is an assembly in fasta format; can be bgzip compressed

This tool was developed for [O'Donnell et al. 2025]() (please cite this publication if you find this tool useful)

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




## Manually identified a secondary repeat?
When looking a the alignments of the contig ends (and maybe manually reordering the ends etc), did you notice another region next to telomeres with several alignments? <br/>
Here is a few lines you can run in order to extract and detail another repeat <br/>
All you need to know is the coordinates for one of the repeats (will be considered the representative copy) <br/>
For Example: I extracted a secondary repeat in the assembly 'GCA030345115.fa' from contig 'GCA030345115_CP128282.1' between position 2 and 9881 (GCA030345115_CP128282.1:2-9881) <br/>
In my run of TLHcrsFinder on 'GCA030345115.fa', the prefix was GCA030345115 and therefore in the TLHcrsFinder output folder I ran the below to add the secondary repeat:

        ##my run of TLHcrsFinder commented out
        #TLHcrsFinder -a GCA030345115.fa -p GCA030345115 -o GCA030345115_TLHcrsFinder_output
        
        prefix="GCA030345115"
        assembly="GCA030345115.fa"
        coords="GCA030345115_CP128282.1:2-9881"

        ##move into the TLHcrsFinder output for 'GCA030345115.fa'
        cd GCA030345115_TLHcrsFinder_output/
        
        ##manually extract the second repeat for ${prefix} based on using the synteny alignment script (only the primary repeat was found previously)
        samtools faidx ${assembly} ${coords} > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.fa
        ##get the positions of this second repeat in the contig ends
        blastn  -subject contig_ends/${prefix}.${tipsize2}kb_ends.fa -query subtelomeric_repeats/${prefix}.manual_second.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.tsv
        ##create bed from nonredundant positions (used for plotting)
        echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.bed
        cat subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10  >> subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.bed
        ##and the same now for the whole genome
        blastn  -subject assemblies/${prefix}.fa -query subtelomeric_repeats/${prefix}.manual_second.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.tsv
        echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.bed
        cat subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10 >> subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.bed

        ##Now you can plot both side by side






