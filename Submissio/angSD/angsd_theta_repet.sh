#!/bin/bash

# OLD # qsub -cwd -V -N 'heteroz_eng' -j y -q h0809.q -l hostname=hercules09 -b y bash angsd_theta.sh Deng_reads.bam Dsil_genome.fasta

# UPDATED # qsub -cwd -V -N 'chr_U1_repet' -j y -q h11.q -l hostname=hercules11 -pe ompi64h11 24 -b y bash angsd_theta_repet.sh Deng_reads.bam Dsil_genome.fasta 'chr' 'min' 'max'


## SOFT path
path_angsd_H=/soft/angsd
path_angsd_A=/home/adria/Program_bioinfo/angsd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_angsd=$path_angsd_H
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Parameters
ncpu=50
maxdepth=$5
mindepth=$4
window_size=50000
step_size=50000

# angsd output folder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out_name=${3##*/}_mn$4_mx$5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

outfolder=output/${out_name}_RepeatMasker
mkdir -p $outfolder

# Input files
# Remember to index them with:
# samtools index *.bam
# samtools faidx *.fasta
bam_file=$1
ref_genome=$2

# scaffold of interest: A -rf formatted file from RepeatMasker.
# See angsd wiki on -rf formatting: must be "chr:1-10 \n chr:20-30 \n etc"

scf_id=$3

# "${path_variable##*/} # Talla l'string del path i usa solament l'últim directori"
# Si el path és 1/2/3, cridant com a dalt recull solament --> 3.

echo '#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#'
echo '#'
echo "# Estimation of Site Allele Frequency"
echo "# BAM file: ${bam_file##*/}"
echo "# Anc genome: ${ref_genome##*/}"
echo "# Ref genome: ${ref_genome##*/}"
echo "# SAF output file: ${outfolder}/${out_name}.glob.saf.idx"
echo "# Selected RepeatMasker File: ${scf_id}"
echo '#'
echo '#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#|#'

###########################################
# - Treure el SAF per a tot el genoma
# - Filtratge basat en la qualitat dels reads i en la profunditat del coverage
# - Ref/Anc: DsilV2.3

$path_angsd/angsd \
	-GL 2 \
	-out $outfolder/$out_name.glob \
	-nThreads $ncpu \
	-doSaf 1 \
	-doCounts 1 \
	-anc $ref_genome \
	-ref $ref_genome \
	-C 50 \
	-baq 1 \
	-minMapQ 20 \
	-minQ 20 \
	-uniqueOnly 1 \
	-remove_bads 1 \
	-only_proper_pairs 1 \
	-setMinDepth $mindepth \
	-setMaxDepth $maxdepth \
	-i $bam_file \
	-rf $scf_id

##########################################
# Per a un cromosoma concret d'interès
# - Treure el SAF
# - Treure els 10 GL amb -doGlf 4

ofolder=$outfolder/glftest
mkdir -p $ofolder

# No cal fer-ho; comment out:
if false; then

$path_angsd/angsd \
  -GL 2 \
  -out $ofolder/$out_name.glftest \
  -nThreads $ncpu \
  -doSaf 1 \
  -doCounts 1 \
  -doGlf 4 \
  -anc $ref_genome \
  -ref $ref_genome \
  -rf $scf_id \
  -C 50 \
  -baq 1 \
  -minMapQ 20 \
  -minQ 20 \
  -uniqueOnly 1 \
  -remove_bads 1 \
  -only_proper_pairs 1 \
  -setMinDepth $mindepth \
  -setMaxDepth $maxdepth \
  -i $bam_file

$path_angsd/misc/realSFS print $ofolder/$out_name.glftest.saf.idx -rf $scf_id > $ofolder/$out_name.glftest.saf.idx.txt
fi

########################################3

# Calculate the unfolded SFS
echo "## SFSunfolded estimation"
echo "## - output file: ${out_name}.glob.sfs"

$path_angsd/misc/realSFS \
	$outfolder/$out_name.glob.saf.idx \
	-P $ncpu \
	> $outfolder/$out_name.glob.sfs


# Calculate the folded SFS
echo "## SFSfolded estimation"
echo "## - output file: ${out_name}.glob.FOLDED.sfs"

$path_angsd/misc/realSFS \
	$outfolder/$out_name.glob.saf.idx \
	-fold 1 \
	-P $ncpu \
	> $outfolder/$out_name.glob.FOLDED.sfs


# angsd parameters
# -GL 2: # GATK GL
# -anc: # ancestral fasta
# -ref: # reference fasta
# -bam : # bam file
# -doGlf 4: veure tots els possibles genotype likelihoods.
# -fold 1: * NB the ancestral state needs to be supplied for the full SFS, but you can use the -fold 1 to estimate the folded SFS and then use the reference as ancestral.
# -doCounts 1 # The doCounts 1 options for allele counts is needed in order to determine the most common base.
#			   If multiple individuals are used the four bases are counted across individuals. Ns or filtered based are ignored
# -dumpCounts 1 # Printing Counts per site. '1' Print overall depth in the .pos file. This depth is the sum of reads covering a sites for all individuals. 
# 				The first column is the chromosome, the second is the position, the third is the total depth. 
# -doMajorMinor 1 # Infer major and minor from GL

# -setMinDepth NUM # (If total depth is larger then site is removed from analysis. -1 indicates no filtering). !! Supply "-doCounts 1"
# -setMaxDepth NUM # (If total depth is smaller then site is removed from analysis. -1 indicates no filtering). !! Supply "-doCounts 1"

# -C 50 # adjust mapQ for excessive mismatches (as SAMtools), !! Supply "-ref"
# -baq 1 # adjust qscores around indels (1=normal baq 2= extended(as SAMtools)), !! Supply "-ref"

# -uniqueOnly 1 # Remove reads that have multiple best hits. 0 no (default), 1 remove.
# -remove_bads 1 # Discard 'bad' reads, (flag & 512)
# -only_proper_pairs 1 # Only use reads where the mate could be mapped	


ofolder=$outfolder/Results_stats
mkdir -p $ofolder

### Calculate global theta
echo "## Step1. Calculate global theta ##"
echo "## - Using Unfolded SFS - ##"
echo "## -- Output file: ${out_name}.glob.sfs"

$path_angsd/misc/realSFS \
    saf2theta \
    $outfolder/$out_name.glob.saf.idx \
    -P $ncpu \
    -outname $outfolder/$out_name.glob.theta_per_site \
    -sfs $outfolder/$out_name.glob.sfs # unfolded SFS

 echo "## - Using Folded SFS - ##"
 echo "## -- Output file: ${out_name}.glob.FOLDED.sfs"
 $path_angsd/misc/realSFS \
     saf2theta \
     $outfolder/$out_name.glob.saf.idx \
     -outname $outfolder/$out_name.glob.FOLDED.theta_per_site \
     -fold 1 \
     -sfs $outfolder/$out_name.glob.FOLDED.sfs # folded SFS

####Calculate theta and summary statistics (slidding window analysis) (UNFOLDED)
echo "## Step2. Calculate the thetas for each site ##"
echo "## -- Output file: ${out_name}.glob.theta_per_site*"

echo "## 2.1 Estimate Tajimas D and other statistics from global theta."
$path_angsd/misc/thetaStat do_stat \
    $outfolder/$out_name.glob.theta_per_site.thetas.idx \
    -outnames $ofolder/$out_name.stat_global

# no admet els paràmetres -r i -P

echo "## 2.2 Estimate Tajimas D and other statistics from selected scaffold."
$path_angsd/misc/thetaStat do_stat \
    $outfolder/$out_name.glob.region.theta_per_site.thetas.idx \
    -outnames $ofolder/$out_name.stat_region

# Calculate thetas and statistics per window!!!!
echo "## 2.3 Estimate thetas and statistics per window"
echo "## -- Output file: ${out_name}.win.gz"

$path_angsd/misc/thetaStat do_stat \
    $outfolder/$out_name.glob.theta_per_site.thetas.idx \
    -win $window_size \
    -step $step_size \
    -outnames $ofolder/$out_name.win \
    -type 2

# # Print theta values for a given chromosome
$path_angsd/misc/thetaStat print $outfolder/$out_name.glob.theta_per_site.thetas.idx > $ofolder/chromosome.globthetas.txt

