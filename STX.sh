# Script is based on: https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh


#####----------- DOWNLOADING RAW READS FROM BASESPACE ---------###############
# Note: You must have a BaseSpace account 
# Download packages necessary to access our reads from BaseSpace
# In TACC Terminal:
cd
cd bin
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs
chmod u+x bs
bs auth # authenticate, should give you a url to copy to browser
bs list projects # Should list project numbers to download

# Should see Job ID listed
# For STX coral projext, Job Name is: JA23218. The ID# should be the second column.

# Make a new directory in Scratch to store your reads
cds 
mkdir STX
cd STX
echo "bs download project -i '391841321' -o /scratch/07090/kblack/STX --extension=fastq.gz" >bs
ls6_launcher_creator.py -j bs -n bs -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 12 -q normal
sbatch bs.slurm

# Move fastq files from subdirectories to top-level directory
find . -name '*.fastq.gz' -exec mv -it . {} +

# Backup raw reads to Stockyard (aka Work2)
mkdir $STOCKYARD/STX 
cp *fastq.gz $STOCKYARD/STX







#####----------- DOWNLOADING RAW READS SHARED DIRECTORY ON TACC---------###############
# If no BaseSpace account, and the reads are in a shared directory, you can just copy them over:
# The shared directory must have read/write permission set for all users in the same group
# All users should also uncomment umask 027 in their bashrc

cd $SCRATCH/STX
cp $STOCKYARD/STX2/*fastq.gz .
ll # Should see all the fastq files now in your scratch directory





#####----------- Deduplicate and rename reads ---------###############

# Unzip fastq files for demultiplexing
>gunz
for F in `ls *gz`; do echo "gunzip $F" >> gunz; done
ls6_launcher_creator.py -j gunz -n gunz -t 0:15:00 -w 48 -N 2 -a IBN21018 -e kblack@utexas.edu 
sbatch gunz.slurm

# Concatenate two lanes of reads for each sample
echo "ngs_concat.pl 'fastq' '(.+)_S'" > concat
ls6_launcher_creator.py -j concat -n concat -t 2:00:00 -w 6 -N 2 -a IBN21018 -e kblack@utexas.edu -q normal
sbatch concat.slurm

# Trim, deduplicate and demultiplex samples
2bRAD_trim_launch_dedup.pl fq sampleID=1 > trims
cat trims | sed 's/$/ minBCcount=100000/' > trims2
ls6_launcher_creator.py -j trims2 -n trims2 -t 2:00:00 -N 4 -w 8 -a IBN21018 -e kblack@utexas.edu -q normal
sbatch trims2.slurm
# 45 samples unsorted

# Need to upload a separate table (SampleID.txt) with actual file names. 
# Open a new terminal window locally and secure copy (scp) the text file from your local drive to our working TACC directory of deduplicated files
scp STX_SampleID.txt kblack@ls6.tacc.utexas.edu:/scratch/07090/kblack/STX

# Check that the file is now in TACC, then run this line to change file names
# column 1 of SampleID.txt will be the deduplicated file name, column 2 will be whatever you want to change it to
awk -F'\t' 'system("mv " $1 " " $2)' STX_SampleID.txt


# Split samples into separate directories by species
# scp STX_SpeciesID.txt to our working TACC directory (using a separate terminal window)
awk -F'\t' 'system("echo " $2)' STX_SpeciesID.txt
awk -F'\t' 'system("mkdir -p " $2)' STX_SpeciesID.txt
awk -F'\t' 'system("mv " $1 " " $2)' STX_SpeciesID.txt

# Backup sorted fastq directories to Stockyard (Change the name of your folder to your spcies)
cp -R FLKeys /work/07090/kblack/STX2



# Go into your species' folder
cd Bonnetheads







#####----------- PREPROCESSING YOUR READS ---------###############

# Install cutadapt: 
cdh
pip install --user cutadapt
cp .local/bin/cutadapt ~/bin



# for reference-based analysis: trimming poor quality bases off ends:
>trimse
for file in *.fastq; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.fastq/}.trim $file > ${file}_trimlog.txt" >> trimse;
done
ls6_launcher_creator.py -j trimse -n trimse -t 0:30:00 -w 48 -N 2 -a IBN21018 -e kblack@utexas.edu 
sbatch trimse.slurm
# do we have expected number of *.trim files created?
ls -l *.trim | wc -l
# If yes, skip to part: "IF YOUR SPECIES HAS A REFERENCE GENOME"


# for de novo analysis: removing reads with qualities at ends less than Q15
>trimse
for file in *.fastq; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.fastq/}.trim $file > ${file}_trimlog.txt" >> trimse;
done
ls6_launcher_creator.py -j trimse -n trimse -t 0:30:00 -w 48 -N 2 -a IBN21018 -e kblack@utexas.edu  
sbatch trimse.slurm
# do we have expected number of *.trim files created?
ls -l *.trim | wc -l
# If yes, skip to part: "IF YOUR SPECIES DOES NOT HAVE A REFERENCE GENOME, WE WILL NEED TO MAKE ONE DE NOVO"







#####----------- IF YOUR SPECIES HAS A REFERENCE GENOME ---------###############


#import reference genome for your species from NCBI
cd $STOCKYARD
cd db
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/679/065/GCA_024679065.1_ASM2467906v1/GCA_024679065.1_ASM2467906v1_genomic.fna.gz .
gunzip GCA_024679065.1_ASM2467906v1_genomic.fna.gz







#####----------- IF YOUR SPECIES DOES NOT HAVE A REFERENCE GENOME, WE WILL NEED TO MAKE ONE DE NOVO ---------###############

# 'uniquing' ('stacking') individual trimmed fastq reads:
ls *.trim | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii
ls6_launcher_creator.py -j unii -n unii -t 2:00:00 -w 48 -N 2 -a IBN21018 -e kblack@utexas.edu
sbatch unii.slurm
# Done! do you have .uni for all your samples?... 
ls -l *.uni | wc -l  

# collecting common tags (= major alleles)
# merging uniqued files (set minInd to >10, or >10% of total number of samples, whichever is greater)
mergeUniq.pl uni minInd=10 >all.uniq
# discarding tags that have more than 7 observations without reverse-complement
awk '!($3>7 && $4==0) && $2!="seq"' all.uniq >all.tab
# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab > all.fasta
# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

# making fake reference genome (of 10 chromosomes) out of major-allele tags
concatFasta.pl fasta=cdh_alltags.fas num=10







#####------------------------ BUILDING YOUR GENOME ---------###############

# We will need to import the symbiont genomes and concatenate them onto our coral genome
# This way, we can split out and remove any symbiont reads, and keep only coral reads
# We can copy the symbiont genomes from Carly's shared directory (/work/06909/cbscott/sym_references/):
cp /work/06909/cbscott/sym_references/sym* .
cat  cdh_alltags_cc.fasta symA_genomic_cc.fasta symB_genomic_cc.fasta symC_genomic_cc.fasta symD_genomic_cc.fasta > Coral_symABCD.fasta
grep ">chr" Coral_symABCD.fasta # If it shows chr20-23 then you have all 4 symbionts

# Index your genome
export GENOME_FASTA=cdh_alltags_cc.fasta
echo "bowtie2-build $GENOME_FASTA $GENOME_FASTA" >btb2
ls6_launcher_creator.py -j btb2 -n btb2 -t 1:00:00 -w 48 -N 4 -a IBN21018 -e kblack@utexas.edu  
sbatch btb2.slurm
# This line is short and can be run on log-in node:
samtools faidx $GENOME_FASTA







#####---------------------------- MAPPING TO YOUR GENOME -----------------###############

# Mapping your reads to the genome, converting to bams, indexing
export GENOME_FASTA=cdh_alltags_cc.fasta
>maps
for F in `ls *.trim`; do
echo "bowtie2 --local --no-unal -x $GENOME_FASTA -U ${F} -S ${F/.trim/}.sam && samtools sort -O bam -o ${F/.trim/}.bam ${F/.trim/}.sam && samtools index ${F/.trim/}.bam">>maps
done
ls6_launcher_creator.py -j maps -n maps -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 48 -N 4 -q normal
sbatch maps.slurm

# Get mapping stats
# Do this in idev:
>alignmentRates
for F in `ls *trim`; do 
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done
#scp alignmentrates

# Check how many reads mapped to each genome
conda activate samtools
>idx
for F in `ls *.bam`; do
echo "samtools idxstats ${F} >> idxstats">>idx
done
ls6_launcher_creator.py -j idx -n idx -a IBN21018 -e kblack@utexas.edu -t 1:00:00
sbatch idx.slurm
# scp idxstats






#####---------------------------- REMOVING SYMBIONT READS -----------------###############

# Split out coral reads to keep (remove reads mapping to symbiont genomes)
# First, make an index file for coral reads only
samtools faidx cdh_alltags_cc.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' cdh_alltags_cc.fasta.fai > Coral.bed

# Make a new directory for splitting out your reads
mkdir split
mv *bam* split # move your bams over
cp *bed split # copy your bed file over
cd split # working in subdirectory now

>split
for F in `ls *.bam`; do
echo "samtools view -L Coral.bed -o ${F/001.bam/}_Coral.bam $F" >>split
done
ls6_launcher_creator.py -j split -n split -a IBN21018 -e kblack@utexas.edu -t 00:10:00 -w 48 -q normal
sbatch split.slurm

# Reindex new bams names
>maps2
for F in `ls *_Coral.bam`; do
echo "samtools index $F">>maps2
done
ls6_launcher_creator.py -j maps2 -n maps2 -a IBN21018 -e kblack@utexas.edu -t 0:10:00  
sbatch maps2.slurm

# Move your bam files with only coral reads back to your working directory
cp *_Coral.bam* /scratch/07090/kblack/STX
cd /scratch/07090/kblack/STX

# Remove _Coral from all file names
mkdir rename
cp *_Coral.bam rename
cd rename
for filename in *.bam; do 
    [ -f "$filename" ] || continue
    mv "$filename" "${filename//_Coral/}"

done







#####-------------------------------- GENOTYPING YOUR SAMPLES --------------------###############

# quality assessment, removing bams with log(coverage)<3SD
# also imposing minimum number of individuals(MI) a locus must be seen in (genotyping rate cutoff - 50%)
# Use angsd installed locally (not conda version)
module load Rstats
# Change -minInd to 50% of your sample size
FILTERSQ="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1000 -minInd 43"
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo "ls *.bam > bams && angsd -b bams -GL 1 $FILTERSQ $TODOQ -P 1 -out dd && Rscript ~/bin/plotQC.R prefix=dd >qualRanks">a0
ls6_launcher_creator.py -j a0 -n a0 -a IBN21018 -e kblack@utexas.edu -t 1:00:00  
sbatch a0.slurm
# look at quality of reads in dd.pdf, check how many bams retained in bams.qc

# Initial IBS production, detecting and removing clones:
# Change STX in the next few lines to your species name
# Set minInd to 75-80% of your total number of bams
export MinIndPerc=0.75
FILTERS0="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd $MI -snp_pval 1e-5 -minMaf 0.05"
TODO0="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"
echo 'export NIND=`cat bams.qc | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.qc -GL 1 $FILTERS0 $TODO0 -P 12 -out STX && Rscript ~/bin/detect_clones.R bams.qc GBR.ibsMat 0.15">a1
ls6_launcher_creator.py -j a1 -n a1 -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 1
sbatch a1.slurm
# Look for clones in hctree.pdf, check how many bams retained in bams.nr

# Final IBS production:
# Change STX in the next few lines to your species name
# Set minInd to 75-80% of your total number of bams
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 1000 -snp_pval 1e-5 -minMaf 0.05"
TODO1="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"
echo 'cat bams.nr | sort > bams.NR && mv bams.NR bams.nr && export NIND2=`cat bams.nr | wc -l`; export MI2=`echo "($NIND2*$MinIndPerc+0.5)/1" | bc`' >calc2
echo "source calc2 && angsd -b bams.nr -GL 1 $FILTERS1 $TODO1 -P 12 -out STX && Rscript ~/bin/pcaStructure.R STX.ibsMat > pcaStruc.txt">a2
ls6_launcher_creator.py -j a2 -n a2 -a IBN21018 -e kblack@utexas.edu -t 2:00:00 -w 1 
sbatch a2.slurm
# Check how many sites retained:
NSITES=`zcat myresult.mafs.gz | wc -l`
echo $NSITES


# scp *Mat, *qopt and bams files to laptop, use angsd_ibs_pca.R to plot PCA and admixturePlotting_v4.R to plot ADMIXTURE
