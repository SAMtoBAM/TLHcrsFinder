#!/bin/bash
set -euo pipefail


##############################################################
##################### SETTING VARIABLES ######################
##############################################################

#default values, unless denoted when running PAQman
assembly=""
tipsize="50000"
telomererepeat="TTAGGG"
covthreshold="0.75"
sizemin="2000"
threads="1"
window="10"
slide="5"


prefix="TLHcrsFinder"
output="TLHcrsFinder_output"
help="nohelp"

## to clean up a bunch of output from the tools in order to reduce all the unnecessary output
cleanup="yes"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
    -a|--assembly)
    assembly="$2"
    shift
    shift
    ;;
    -ts|--tipsize)
    tipsize="$2"
    shift
    shift
    ;;
    -tr|--telomererepeat)
    telomererepeat="$2"
    shift
    shift
    ;;
    -ct|--covthreshold)
    covthreshold="$2"
    shift
    shift
    ;;
    -sm|--sizemin)
    sizemin="$2"
    shift
    shift
    ;;
    -t|--threads)
    threads="$2"
    shift
    shift
    ;;
    -w|--window)
    window="$2"
    shift
    shift
    ;;
    -s|--slide)
    slide="$2"
    shift
    shift
    ;;
    -p|--prefix)
    prefix="$2"
    shift
    shift
    ;;
    -o|--output)
    output="$2"
    shift
    shift
    ;;
    -h|--help)
    echo "
 
    TLHcrsFinder.sh -a assembly.fa
    
    Required inputs:
    -a | --assembly     Genome assembly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz)

    Recommended inputs:
    -ts | --tipsize     Length of contig ends to be extracted for TLHcrs detection (Default: 50000)
    -tr | --telomererepeat       Telomeric repeat pattern (default: TTAGGG)
    -ct | --covthreshold    The amount of coverage required for a region to be considered for TLHcrs clustering relative to the number fo telomeres (Default = 0.75)
    -sm | --sizemin     The minimum size of a region passing the coverage threshold to be considered as a potential TLHcrs region (Default: 2000)
    
    Optional parameters:
    -w | --window       Number of basepairs for window averaging for coverage (default: 10)
    -s | --slide        Number of basepairs for the window to slide for coverage (default: 5)
    -p | --prefix       Prefix for output (Default: TLHcrsFinder)
    -o | --output       Name of output folder for all results (default: TLHcrsFinder_output)
    -h | --help         Print this help message
    "
    exit
    ;;
    *)  # catch invalid args
    echo "ERROR: Unknown option: '$1'"
    echo "Run 'TLHcrsFinder.sh -h' to see valid options"
    exit 1
    ;;
    esac
done


##check the assembly is given and that I can find it
[[ $assembly == "" ]] && echo "ERROR: Path to assembly not found, assign using -a" && exit
assemblypath=$( realpath ${assembly} )
[ ! -f "${assemblypath}" ] && echo "ERROR: Cannot find path to assembly file provided by -a; check path is correct and file exists" && exit

##warning if the telomeric repeat and/or number of threads was not modified
[[ $telomererepeat == "TTAGGG" ]] && echo "WARNING: Using default telomeric repeat sequence (TTAGGG)"
[[ $threads == "1" ]] && echo "WARNING: Using 1 thread for analysis"

##get a kb relative tipsize for output file tags
tipsize2=$( echo ${tipsize} | awk '{print $1/1000}' )

##check if assembly is compressed or not
compressed=$( echo ${assembly} | awk -F "." '{if($NF == "gz" || $NF == "bgzip") {print "yes"} else {print "no"}}' )

#######################################################################################################
####################################BEGIN FINDING THE TLHCRS REPEAT#################################### 
#######################################################################################################


echo "######### Starting TLHcrsFinder"
echo "######### Output in ${output}"


##make an output directory and move into it
mkdir ${output}
cd ${output}


mkdir contig_ends
mkdir contig_ends_alignments
mkdir telomeric_repeats
mkdir contig_ends_coverage
mkdir subtelomeric_repeats
mkdir plotting_Rscripts

##create link to assembly in the output folder
ln -sf ${assemblypath} ./
assembly=$( echo ${assembly} | awk -F "/" '{print $NF}' )

##get the sametools index as a quick way to generate a file with contig sizes etc
samtools faidx ${assembly}

##getting the tips of contigs larger than 2 times the tip size
cat ${assembly}.fai  | awk -v tipsize="$tipsize" '{if($2 > (2*tipsize)) {print $1"\t"$2} }' | while read line
do
contig=$( echo "${line}" | awk -F "\t" '{print $1}' )
size=$( echo "${line}" | awk -F "\t" '{print $2}' )
if [ $size -gt ${tipsize} ]
then
end=$( grep ^${contig} ${assembly}.fai | awk -F "\t" -v contig="$contig" '{if($1 == contig) print $2}' | awk -v tipsize="$tipsize" '{if(($1-tipsize)> 0) {print ($1-tipsize)"-"$1} else {print "1-"$1}}' )
samtools faidx ${assembly} ${contig}:1-${tipsize}
samtools faidx ${assembly} ${contig}:${end}
fi
done > contig_ends/${prefix}.${tipsize2}kb_ends.fa

##can also label all regions with the canonical telomeric repeat
##use seqkit to locate the position of the canonical repeat then merge all those locations with a buffer of 7bp incase one repeat is off and export a bed file
echo "contig;start;end;sense" | tr ';' '\t' > telomeric_repeats/${prefix}.telomeres.bed
seqkit locate --ignore-case -p "TTAGGG" ${assembly} | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d 7 -c 4 -o distinct -i - | awk '{if($3-$2 > 11) print}' | awk -F "," '{print $1}' >> telomeric_repeats/${prefix}.telomeres.bed

echo "contig;start;end;sense" | tr ';' '\t' > telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed
seqkit locate --ignore-case -p "TTAGGG" contig_ends/${prefix}.${tipsize2}kb_ends.fa | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d 7 -c 4 -o distinct -i - | awk '{if($3-$2 > 11) print}' | awk -F "," '{print $1}' >> telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed

##subset the ends to only those with a good telomere sequence (>50bp)
#cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  | awk '{if($3-$2 > 50) print $1}' > contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.txt
#seqkit grep -f contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.txt contig_ends/${prefix}.${tipsize2}kb_ends.fa > contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.fa

##get simple bed file of just the ends to use for visualisation
samtools faidx contig_ends/${prefix}.${tipsize2}kb_ends.fa
echo "contig;start;end" | tr ';' '\t' > contig_ends/${prefix}.${tipsize2}kb_ends.bed
cat contig_ends/${prefix}.${tipsize2}kb_ends.fa.fai | awk -F "\t" '{print $1"\t1\t"$2}' >> contig_ends/${prefix}.${tipsize2}kb_ends.bed

##align the whole regions to themselves (once with all ends and another with just telomeric ends)
nucmer --maxmatch --threads ${threads} --delta contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.delta contig_ends/${prefix}.${tipsize2}kb_ends.fa contig_ends/${prefix}.${tipsize2}kb_ends.fa
##convert to paf format and remove self matches and the reciprical alignments
paftools.js delta2paf contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.delta | awk '{
    q = $1; r = $6;
    qs = $3; qe = $4;
    rs = $8; re = $9;

    if (q == r) next;  # Skip self-alignments (same query and reference)

    # Canonical key: always sort (query, reference) alphabetically
    if (q < r) {
        key = q "|" qs "|" qe "|" r "|" rs "|" re;
    } else {
        key = r "|" rs "|" re "|" q "|" qs "|" qe;
    }

    if (!seen[key]++) print;
}' > contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf


##make a bedfile in order to use for coverage analysis of the contig ends alignments
cat contig_ends/${prefix}.${tipsize2}kb_ends.fa.fai | awk '{print $1"\t1\t"$2}' | bedtools makewindows -w 10 -s 5 -b - > contig_ends/${prefix}.${tipsize2}kb_ends.10bpwindow.bed

##create a bed file from the alignments
cat contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf  | awk '{if($3 > $4) {print $1"\t"$4"\t"$3} else {print $1"\t"$3"\t"$4}}' > contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf.bed

##get all regions with at least XX coverage and larger than XXkb
##to get the minimum coverage we could estimate this by taking the number of contigs with +50bp telomeric sequences found and want at least 75% that coverage (unless 0 then make it 4)
goodtel=$( cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  | awk '{if($3-$2 > 50) print}' | awk -F "\t" '{print $1}' | wc -l )
covmin=$( echo "${goodtel}" | awk -v covthreshold="$covthreshold" '{if($1 != "0") {print $1*covthreshold} else {print "4"}}' )
##size min set to 2kb as default 
bedtools coverage -a contig_ends/${prefix}.${tipsize2}kb_ends.10bpwindow.bed -b contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf.bed > contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.cov.bed
echo "contig;start;end" | tr ';' '\t' > contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed
cat contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.cov.bed | awk -v covmin="$covmin" '{if($4 > covmin) print}' | bedtools merge -d 10 | awk -v sizemin="$sizemin"  '{if($3-$2 > sizemin) print}' >> contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed

##we shall now have identified any repeats that a common across these contig ends
##now we want to extract a reference sequence for the repeat based on all the coordinates identifed
tail -n+2 contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed | bedtools getfasta -fi contig_ends/${prefix}.${tipsize2}kb_ends.fa -bed - -fo contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa

##then generate clusters based on a 80/80/80 threshold
##and selecting a representative per cluster (the longest sequence meeting the threshold)
cd-hit-est -i contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa -o contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500
##now we have the representative and clusters, we want to select the cluster with the great number of alignments (removing some smaller regions)
##then we can BLAST for our representative across the regions again and get a better indicator of its appearance.
cluster=$( awk '/Cluster/{c++} />/{print "Cluster"c-1}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | sort | uniq -c | sort -nr  | head -n1 | awk '{print $2}' )
awk '/Cluster/{c++} />/{print "Cluster"c-1"\t"$0}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | grep "^${cluster}" | awk '{print $4}' | sed 's/\.\.\.//' | sed 's/>//' > subtelomeric_repeats/${prefix}.repeat_cluster.txt
rep=$( awk '/Cluster/{c++} />/{print "Cluster"c-1"\t"$0}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | grep "^${cluster}" | grep "\*$" | awk '{print $4}' | sed 's/\.\.\.//' | sed 's/>//' )
##extract just the rep
seqkit grep -p "${rep}" contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa > subtelomeric_repeats/${prefix}.repeat_rep.fa
seqkit grep -f subtelomeric_repeats/${prefix}.repeat_cluster.txt contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa > subtelomeric_repeats/${prefix}.repeat_cluster.fa

##now use the representative to BLAST against the ends again and create a new bed file for the coords
##just need some filters...
##an alignment of at least 80% identity and 500bb
blastn  -subject contig_ends/${prefix}.${tipsize2}kb_ends.fa -query subtelomeric_repeats/${prefix}.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.tsv
##create bed from nonredundant positions
echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed
cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10  >> subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed

##and the whole genome
if [[ $compressed == "yes" ]]
then
bgzip -dc ${assembly} | blastn -subject - -query subtelomeric_repeats/${prefix}.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv
else
blastn -subject ${assembly} -query subtelomeric_repeats/${prefix}.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv
fi
echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed
cat subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10 >> subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed
##get bed file for the whole genome so we can plot the positions of the repeats alongside it
echo "contig;start;end" | tr ';' '\t' > ${prefix}.genome.bed
cat ${assembly}.fai | awk '{print $1"\t1\t"$2}'  >> ${prefix}.genome.bed


cat Finderplots_tlhcrs.R | sed "s/SAMPLE/${prefix}/g" > plotting_Rscripts/${prefix}.R

Rscript plotting_Rscripts/${prefix}.R

done
