<p align="center" >
    <img src="https://github.com/SAMtoBAM/TLHcrs_finder/blob/main/logo/TLHcrsFinder_logo.png" width=100%>
</p>


# TLHcrsFinder
TLHcrsFinder is a tools desgined to detect and analyse **T**elomere-**L**inked-**H**elicase **C**ontaining Region**S** (TLHcrs) repeats

TLHcrs repeats are, as the name describes, regions adjacent to telomeres that are conserved across chromosome ends and generally contain helicase genes <br/>
These repeats, such as the most perhaps most well known Y-prime region in _S. cerevisiae_, are conserved across diverse fungi <br/>
With the advent of Long-read sequencing we now have an increasing number of well assembled genomes with the subtelomeres intact <br/>
Therefore we are now primed for looking at the TLHcrs repeats and their evolution.

All that is required to run TLHcrsFinder is an assembly in fasta format (can be bgzip compressed) or a tsv file containing a list of assemblies (first column: sample name; second column: assembly) <br/>
NOTE: Do not use hyphens '-' or other special characters in the sample/file names

This tool was developed for [O'Donnell et al. 2025]() (please cite this publication if you find this tool useful)

# Conda installation

    conda install samtobam::tlhcrsfinder

# How to use

    TLHcrsFinder.sh -a assembly.fa
    ##or
    TLHcrsFinder.sh -al list.tsv

    Required inputs:
    -a | --assembly     Genome assembly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz) with bgzip
    or
    -al | --assemblylist     A tsv file containing sample names in the first column and assembly paths in the second column

    Recommended inputs:
    -ts | --tipsize     Length of contig ends to be extracted for TLHcrs detection (Default: 50000)
    -tr | --telomererepeat       Telomeric repeat pattern (Default: TTAGGG)
    -ct | --covthreshold    The amount of coverage required for a region to be considered for TLHcrs clustering relative to the number fo telomeres (Default = 0.75)
    -sm | --sizemin     The minimum size of a region passing the coverage threshold to be considered as a potential TLHcrs region (Default: 2000)
    
    Optional parameters:
    -w | --window       Number of basepairs for window averaging for coverage (Default: 10)
    -s | --slide        Number of basepairs for the window to slide for coverage (Default: 5)
    -p | --prefix       Prefix for output (Default: TLHcrsFinder)
    -o | --output       Name of output folder for all results (default: TLHcrsFinder_output)
    -h | --help         Print this help message



## How does TLHcrsFinder work

TLHcrsFinder works in 10 main steps (each step is run on all assemblies provided -al)

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


## How to look for known functional domains
TLHcrsFinder looks for the conserved region however the actual expressed portion of the repeat (generally a helicase) is only contained within the repeat <br/>
Therefore you may be interested to find domains within the repeat as an indication of the function <br/>
Notably, in [O'Donnell et al. 2025]() using TLHcrsFinder we found several examples of repeats containing genes that were not helicases.

This step requires you to download a large database of conserved domains in order to search against <br/>

        ##first download the appropriate conserved domain library (https://pmc.ncbi.nlm.nih.gov/articles/PMC7378889/)
        wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz
        mkdir Cdd_LE
        mv Cdd_LE.tar.gz Cdd_LE
        cd Cdd_LE
        tar -zxvf Cdd_LE.tar.gz

        ##also download a library to read the cdd index numbers
        wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
        gunzip cddid.tbl.gz
        cd ../

Then use rpstblastn to search against these domains with your representative TLHcrs repeat OR your all TLHcrs repeats in your assembly <br/>
Here we are using the TLHcrs repeat representative for GCA030345115 <br/>
We also just grab the top 5 for the blast results for looking at the function using the cdd index table (this can definitely be increased if necessary)

        sample="GCA030345115"
        input="GCA030345115_TLHcrsFinder_output/subtelomeric_repeats/GCA030345115.repeat_rep.fa"
        tophits="5"

        ##get the outfmt 6 plus the actual sequence aigned in the query (our repeats)
        rpstblastn -query ${input} -db Cdd_LE/Cdd -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq" > ${sample}.repeat_rep.rpsblast_output.txt 
 
        ##now take the best 5 hits for each representative
        echo "repeat_representative;CDD;domain;title;description" | tr ';' '\t' > ${sample}.repeat_rep.rpsblast_output.besthit_description.tsv
        grep '>' ${input} | sed 's/>//g' | while read rep
        do
        grep "${rep}" ${sample}.repeat_rep.rpsblast_output.txt | head -n${tophits}
        done | while read line
        do
        cdd=$( echo "${line}" | cut -f2 | awk -F "|" '{print $3}' )
        cdd2=$( grep ^${cdd} Cdd_LE/cddid.tbl | awk -F "\t" '{print "CDD:"$1"\t"$2"\t"$3"\t"$4}' )
        echo "${line}" | awk -v cdd2="$cdd2" '{print $1"\t"cdd2}'
        done  >> ${sample}.repeat_rep.rpsblast_output.besthit_description.tsv



