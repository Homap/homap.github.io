---
layout: post
title: "Mapping to distant genomes"
subtitle: "How distant can we map"
background: '/img/posts/01.jpg'
---

I need to map genomes of very distantly related organism to another. How should I do that?

# Today is Tuesday February 8th 2022

# Script to prepare the genomes and run BUSCO analysis

# gunzip the reference genomes
cd /proj/snic2020-16-269/private/homap/green_algae/data/genomes
for i in GCA*
do 
gunzip $i/${i}_genomic.fna.gz
done

# Change the carriage return into Unix
dos2unix GCA_900128745.1_TOSAG39-1_assembly_report.txt
for i in ../genomes/GCA*
do 
echo $i
dos2unix $i/*report.txt
done

# Create a softlink of all fasta files in fasta_softlinks
cd /proj/snic2020-16-269/private/homap/green_algae/data/fasta_softlinks
for i in ../genomes/GCA*
do 
assembly_name=$(echo ${i}_genomic.fna | cut -f3 -d "/")
echo $assembly_name >> assembly_name_list.txt
ln -s ${i}/${assembly_name}
species=$(grep 'Organism name:' $i/*report.txt | sed -e 's/# Organism name: //g' -e 's/ (green algae)//g' -e 's/ (nom. nud.)//g' -e "s/'/s/g")
echo $species >> assembly_name_list.txt
first_part=$(echo $species | cut -f1 -d " " | cut -c1-3)
second_part=$(echo $species | cut -f2 -d " "| cut -c1-3)
third_part=$(echo $species | cut -f3 -d " ")
if [ "$second_part" = "sp." ]
then
 code_name="${first_part}${second_part^}${third_part^}.genomic.fna"
 echo $code_name >> assembly_name_list.txt
 mv $assembly_name $code_name
else
  code_name="${first_part}${second_part^}.genomic.fna"
  echo $code_name >> assembly_name_list.txt
  mv $assembly_name $code_name
fi
done

for i in ../genomes/lassb_abudhabi_genomes/*.fa
do
echo $i
assembly_name=$(echo $i | cut -f4 -d "/")
first_part=$(echo $i | cut -f4 -d "/" | cut -f1 -d "_" | cut -c1-3)
second_part=$(echo $i | cut -f4 -d "/" | cut -f2 -d "_" | cut -f1 -d "." | cut -c1-3)
code_name=${first_part}${second_part^}.genomic.fna
echo $assembly_name $code_name >> assembly_name_list.txt
ln -s ${i}
mv $assembly_name $code_name
done

# Run BUSCO
mkdir -p slurm_job_IDs
mkdir -p slurm_job_IDs/busco
for genome in ../data/fasta_softlinks/*.fna
do
echo $genome
dir_name=$(echo $genome | cut -f4 -d '/' | sed 's/.genomic.fna//g')
echo $dir_name
sbatch -J $dir_name busco.sh $genome $dir_name > slurm_job_IDs/busco/${dir_name}.out
done

# Some genomes do not belong to Chlorophyceae
# Euglena_viridis.fa EugVir.genomic.fna : Euglena viridis belongs to Euglenoidea and therefore, only 0.8% BUSCOs of chlorophyta were found in it.
# The BUSCO run continued for a over 17 hours. I stop it and remove this from further analysis since it does not belong to Chlorophyceae.
# Lingulodinium_polyedra.fa LinPol.genomic.fna: Lingulodinium polyedra belong to Dinophyceae, only 0.1% BUSCOs of chlorophyta were found in it. 
# I made a mistake initially to include this by miss-reading the class name as Chlorophyceae. Job run took over 17 hours, I stop it and remove it from further analysis.
# Also removing the genomes from data/genomes/lassb_abudhabi_genomes and also from  data/fasta_softlinks/assembly_name_list.txt.
# Lastly, ChlSp.ICE-L has been running for over 18 hours but it still has not produced any summary file. I let this run but continue with the rest of the analysis now since only 1:47 hour is left
# and it might fail so I cancel the job and remove its run folder. I will re-run it with longer hours but for now, let's continue with the rest of the analysis.

# When running BUSCOs, it outputs the result directories starting with "run_" in the bin folder. First copy all directories into the result under busco/busco_output and then remove them from the
# bin directory.
mkdir -p ../result/busco/busco_output
cp -r run_* ../result/busco/busco_output # This takes a very long time so let's do it later but NEED TO BE DONE!

# Get the BUSCO plot for the genomes
mkdir -p ../result/busco/busco_summaries
for f in run_*/short_summary*
do
echo $f
cp $f ../result/busco/busco_summaries
done

export BUSCO_HOME=BUSCO-with-BLASTp
python3 ${BUSCO_HOME}/BUSCO_plot.py -wd ../result/busco/busco_summaries 

# The figure was produced in result/busco/busco_summaries, however, it needed modifications. I removed the busco_figure.png from the directory, copied the busco_figure.R into my local computer,
# opened it in Rstudio and made the modifications for font and figure size there. The resulting figure is busco_short_summaries_plot.pdf located under Dropbox/green_algae/result/busco.

# Extract complete BUSCO genes
# We need to find all the "complete" BUSCO genes predicted from all of the genomes. Since 1519 genes have unique ids, we can quickly get a list of all gene names that are complete and use
# them from downstream analysis.
for f in run_*/full_table_*
do
grep -v "^#" ${f} | awk '$2=="Complete" {print $1}' >> complete_busco_ids.txt
done

# To filter out the complete BUSCOs that are only present less than 3 genomes, we can use uniq and awk commands
sort complete_busco_ids.txt |uniq -c > complete_busco_ids_with_counts.txt
awk '$NF > 2 {print $2}' complete_busco_ids_with_counts.txt > final_busco_ids.txt

# Now let's copy the genes (both nucleotides and amino acids) to different directory. Note that the gene names are just identified as BUSCO id with either faa or fna extension. If you try
# to copy them directly to a single directory, you will likely overwrite all the files and end up with only the last set of files. You can give them a unique name and then merge them all together,
# writing them to a single busco id file. Before you do this, you will also need to edit the sequence name (fasta header) so that it includes organism identifier in it.

mkdir -p ../result/busco/busco_aa
mkdir -p ../result/busco/busco_nt

# This took several hours, it is better to separate into two loops, one for .faa and one for .fna so it take half of the time.
for dir in  run_*/single_copy_busco_sequences
do
  abbrv=$(echo $dir | cut -f1 -d "/" | cut -f2 -d "_")
  echo $abbrv
  for dirfile in ${dir}/*.faa
  do
    #echo $dirfile
    file=$(echo $dirfile | cut -f3 -d "/")
    cp $dirfile ../result/busco/busco_aa/${abbrv}_${file}
    sed -e 's/^>/>'${abbrv}'|/g' -e 's/:.*//g' ../result/busco/busco_aa/${abbrv}_${file} | tr '[:lower:]' '[:upper:]' > ../result/busco/busco_aa/${abbrv}_${file}.1
    mv -f ../result/busco/busco_aa/${abbrv}_${file}.1 ../result/busco/busco_aa/${abbrv}_${file}
    done
  for dirfile in ${dir}/*.fna
  do
  #echo $dirfile
  file=$(echo $dirfile | cut -f3 -d "/")
  cp $dirfile ../result/busco/busco_nt/${abbrv}_${file}
  sed -e 's/^>/>'${abbrv}'|/g' -e 's/:.*//g' ../result/busco/busco_nt/${abbrv}_${file} | tr '[:lower:]' '[:upper:]' > ../result/busco/busco_nt/${abbrv}_${file}.1
  mv -f ../result/busco/busco_nt/${abbrv}_${file}.1 ../result/busco/busco_nt/${abbrv}_${file}
  done
done

# Next step is to generate single file for each of the BUSCO id. We will do this with a simple bash loop. 
# Remember the file final_busco_ids.txt we created previously? we need that for this purpose:
mkdir -p ../result/busco/busco_allSP_AA_fasta
mkdir -p ../result/busco/busco_allSP_nt_fasta
while read line
do
echo ${line}
cat ../result/busco/busco_aa/*_${line}.faa >> ../result/busco/busco_allSP_AA_fasta/${line}_aa.fasta
cat ../result/busco/busco_nt/*_${line}.fna >> ../result/busco/busco_allSP_nt_fasta/${line}_nt.fasta
done < final_busco_ids.txt

# Running one alignment job for each BUSCO ID using MAFFT
mkdir -p ../result/mafft_alignment
mkdir -p slurm_job_IDs/mafft_alignment

for aafasta in ../result/busco/busco_allSP_AA_fasta/*_aa.fasta
do
buscoID=$(echo $aafasta | cut -f5 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID mafft.sh ${aafasta} ../result/mafft_alignment/${buscoID}_aa.aln > slurm_job_IDs/mafft_alignment/${buscoID}.out
done

# Gene tree reconstruction
mkdir -p ../result/raxml_aa
mkdir -p slurm_job_IDs/raxml_aa

# Running raxml script from its own directory since I cannot figure out how to specify output directory path. # On Thursday 10th of February, I ran Raxml on about 30 genes because
# I gave 16 cores per each job and if they each take 1 hour, that would be 16000 core hours so I was afraid I would use up all resources. This is something to discuss.
cp raxml.sh ../result/raxml_aa

cd ../result/raxml_aa

for aa_aln in ../mafft_alignment/subsample/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done


# Species tree reconstruction
# Downloading and installing ASTRAL
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
./make.sh
# Test if it works:
java -jar astral.5.7.8.jar -i Astral/test_data/song_primates.424.gene.tre
# It works! LALA :D
mkdir -p ../result/astral

cat ../result/raxml_aa/RAxML_bestTree.* >> ../result/astral/busco_genes_ML.tree
# Before we proceed with the species tree, we need to fix the title (or the taxa name) in the gene trees. As even single
# tax in these gene trees are identified as genspp|buscoid, they all have unique ids and it will be considered as separate
# taxa. We don't want this. For this, we need to remove the buscoid from the name.

sed -i 's/|[^:]*//g' busco_genes_ML.tree 

# Finally, let's generate a species tree using ASTRAL coalescent-based species tree estimation program. ASTRAL is a 
# fast method for estimating species trees from multiple genes which is statistically consistent, and can run on multiple gene trees efficiently.
export ASTRAL_HOME=ASTRAL

java -jar ${ASTRAL_HOME}/astral.5.7.8.jar --input ../result/astral/busco_genes_ML.tree --output ../result/astral/busco_genes_ASTRAL_spp.tree

# A newick format tree is produced. I upload it into http://www.evolgenius.info/evolview to get a tree figure.

# We end here on Thursday 10th of February - Continue after discussion.

# We are now March 4th and continue with the analyses with full dataset

cp raxml.sh ../result/raxml_aa

cd ../result/raxml_aa

# Create batches
for i in {1..20}
do
mkdir -p batch.${i}
done
# moving files to batch.1 to batch.20 in counts of 60 manually by changing batch.1 to batch.20 in the loop below
i=0
for f in *aln
do
  ((i=i+1))
  echo $i
  if (( $i < 60 ))
  then
    echo $f
    mv $f batch.20
  else
    break
  fi
done

# Submit jobs 1 batch at a time from batch.1 to batch.2 to prevent overuse of resources
for aa_aln in ../mafft_alignment/batch.1/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done

for aa_aln in ../mafft_alignment/batch.2/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done

for aa_aln in ../mafft_alignment/batch.3/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done

for aa_aln in ../mafft_alignment/batch.4/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done

for aa_aln in ../mafft_alignment/batch.5/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done

for aa_aln in ../mafft_alignment/batch.6/*_aa.aln
do
buscoID=$(echo $aa_aln | cut -f4 -d "/" | cut -f1 -d "_")
echo $buscoID
sbatch -J $buscoID raxml.sh ${aa_aln} ${buscoID} > ../../bin/slurm_job_IDs/raxml_aa/${buscoID}.out
done

cat ../result/raxml_aa/RAxML_bestTree.* >> ../result/astral/busco_genes_ML.tree
# Before we proceed with the species tree, we need to fix the title (or the taxa name) in the gene trees. As even single
# tax in these gene trees are identified as genspp|buscoid, they all have unique ids and it will be considered as separate
# taxa. We don't want this. For this, we need to remove the buscoid from the name.

sed -i 's/|[^:]*//g' ../result/astral/busco_genes_ML.tree 

# Finally, let's generate a species tree using ASTRAL coalescent-based species tree estimation program. ASTRAL is a 
# fast method for estimating species trees from multiple genes which is statistically consistent, and can run on multiple gene trees efficiently.
export ASTRAL_HOME=ASTRAL

java -jar ${ASTRAL_HOME}/astral.5.7.8.jar --input ../result/astral/busco_genes_ML.tree --output ../result/astral/busco_genes_ASTRAL_spp.tree

# A newick format tree is produced. I upload it into http://www.evolgenius.info/evolview to get a tree figure.

# Analysis of gene family evolution and diversification


# One extra tip: try adding --p-raxml-version AVX2 to your command. This will greatly speed up your run time.

# Let's take 4 species
# Sphaeropleales: sister group to Chlamymonadales
# Multicellular: Tetradesmus acuminatus - GCA_902809745.2_Scenedesmus-acuminatus-SAG-38.81
# Unicellular: Scenedesmus vacuolatus - GCA_004764505.1_UFZ-Svac-1.0

# Multicellular: Scenedesmus quadricauda: GCA_002317545.1_ASM231754v1
# Unicellular: Monoraphidium neglectum: GCA_000611645.1_mono_v1

# Let's take 4 species with gene annotation (gff)
# C. reinhardtii : GCA_000002595.3_Chlamydomonas_reinhardtii_v5.5_genomic.fna GCF_000002595.2
cp ../genomes/GCA_000002595.3_Chlamydomonas_reinhardtii_v5.5/GCA_000002595.3_Chlamydomonas_reinhardtii_v5.5_protein.faa.gz .
# Tetrabaena socialis : GCA_002891735.1
cp ../genomes/GCA_002891735.1_TetSoc1/GCA_002891735.1_TetSoc1_protein.faa.gz .
# Gonium Pectorale : GCA_001584585.1
cp ../genomes/GCA_001584585.1_ASM158458v1/GCA_001584585.1_ASM158458v1_protein.faa.gz .
# Volvox reticuliferus : GCA_019650235.1
cp ../genomes/GCA_019650235.1_Vretifemale_1.0/GCA_019650235.1_Vretifemale_1.0_protein.faa.gz .

# decompress all protein files
gunzip *gz
# change names 
# each .fasta file must have a name in the form 'xxxx.fasta' where xxxx is a three or four letter unique taxon code.  For example: hsa.fasta or eco.fasta
mv GCA_000002595.3_Chlamydomonas_reinhardtii_v5.5_protein.faa crei.fasta
mv GCA_002891735.1_TetSoc1_protein.faa tsoc.fasta
mv GCA_001584585.1_ASM158458v1_protein.faa gpec.fasta
mv GCA_019650235.1_Vretifemale_1.0_protein.faa vret.fasta

# Create OrthoMCL sqlite database
if [ ! $ORTHO_DB  ] ; then
  ORTHO_DB="orthoMCL.DB.`date +%Y.%m.%d`"; 
fi

echo "dbVendor=sqlite 
dbConnectString=DBI:SQLite:database=${ORTHO_DB}.sqlite
dbLogin=none
dbPassword=none
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5" > orthomcl.config

# drop db to get a clean slate
if [ -e "$ORTHO_DB.sqlite" ] ; then
  echo "Deleting $ORTHO_DB.sqlite"
  rm $ORTHO_DB.sqlite
fi

orthomclInstallSchema orthomcl.config

# Create clean fasta for orthoMCL
module load bioinfo-tools
module load OrthoMCL/2.0.9

DIR=`pwd`
UNIQID=1
# let orthoMCL clean the FASTAs: 
if [ ! -d cleanseq ] ; then
  mkdir cleanseq
  cd cleanseq
  for file in $DIR/proteomes/*.fasta
   do
   echo $file
   base=`basename $file .fasta`
   echo $base
   ## which field has a uniq ID? => 2
   # $base : a three or four letter unique abbreviation for the taxon
   # $file : fasta file
   # $UNIQID : the field that contains the protein ID
   orthomclAdjustFasta $base $file $UNIQID # Create an OrthoMCL compliant .fasta file, by adjusting definition lines.
  done
  cd $DIR
  # Create goodProteins.fasta containing all good proteins and rejectProteins.fasta containing all rejects.
  # input_dir :leanseq - a directory containing a set of .fasta files
  # 10: minimum allowed length of proteins
  # 10 : maximum percent stop codons.
  orthomclFilterFasta cleanseq 10 10 
fi

## format goodProteins BLAST DB
if [ ! -e $DIR/goodProteins.fasta.pin ] ; then
  makeblastdb -in $DIR/goodProteins.fasta -title "OrthoMCLPeps" -parse_seqids -dbtype prot
fi

SPLIT_DEFAULT=50

if [ ! $SPLIT  ] ; then 
  SPLIT=$SPLIT_DEFAULT
fi

## split up goodProteins into smaller FASTA
MAX=`expr $SPLIT - 1`
if [ ! -d $DIR/good_proteins_split ] ; then
  mkdir $DIR/good_proteins_split
fi
if [ ! -e $DIR/good_proteins_split/goodProteins_part_${MAX}.fasta ] ; then
  cd $DIR/good_proteins_split
  cp -s $DIR/goodProteins.fasta .
  split_fasta.pl $SPLIT goodProteins.fasta
  rm -f $DIR/good_proteins_split/goodProteins.fasta
fi 
cd $DIR 

# Set some default values:
BLAST_THREADS=4
BLAST_TIME="10:00:00"
ACCOUNT="snic2022-22-149"
EMAIL="homa.papoli_yazdi@biol.lu.se"

## get ready for blast
## consider adding check for qsub:  if [ `which qsub` != '' ] ; then
if [ ! -d $DIR/blastp_out ] ; then
  mkdir $DIR/blastp_out
#!/bin/bash -l
#SBATCH -A snic2022-22-149
#SBATCH -p core -n 4
#SBATCH -t 10:00:00
#SBATCH -J run_split_blast_array
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load blast/2.12.0+

PART=$1

export DIR=/proj/snic2020-16-269/private/homap/green_algae/data
blastp -query good_proteins_split/goodProteins_part_${PART}.fasta -db goodProteins.fasta -num_threads 4 -outfmt 6 -out $DIR/blastp_out/goodProteins_part_${PART}_vs_goodProteins.BLASTP.tab -evalue 1e-3

# Submit BLAST jobs
for PART in {0..49}; do
    echo submitting job for part $PART
    sbatch $DIR/run_split_blast_array.sh $PART 
done

# Wednesday the 30th, start here
# This will reformat the BLAST data into something to load into the database
if [ ! -f goodProteins.BLASTP.bpo ]; then
 cat blastp_out/*tab > goodProteins.BLASTP.tab
 orthomclBlastParser goodProteins.BLASTP.tab cleanseq > goodProteins.BLASTP.bpo
fi

# Load the data into the DB
orthomclLoadBlast orthomcl.config goodProteins.BLASTP.bpo


finish.orthomcl.sh

# Running orthoMCL
#!/bin/bash -l

#SBATCH -A snic2022-22-149
#SBATCH -p core -n 4
#SBATCH -t 4-00:00:00
#SBATCH -J finish.orthomcl
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL


module load bioinfo-tools
module load OrthoMCL #/2.0.9
# module load mcl/14-137   
# module load sqlite/3.16.2
# module load blast/2.5.0+
# module load perl/5.18.4
# module load BioPerl/1.6.924_Perl5.18.4

cd $DIR 
 
# now run the ortholog/paralog initial finding
rm -rf pairs pairs.log 
orthomclPairs orthomcl.config pairs.log cleanup=no

# Dump out the ortholog groups and the mclInput file
orthomclDumpPairsFiles orthomcl.config

# Run mcl for clustering
mcl mclInput  --abc -I 1.5 -o mclOutput.I15.out

mcl mclInput  --abc -I 2 -o mclOutput.I2.out

mcl mclInput  --abc -I 3 -o mclOutput.I3.out

# convert the MCL clusters into OrthoMCL groups
orthomclMclToGroups OG 1 < mclOutput.I15.out > mclGroups.I15.table
orthomclMclToGroups OG 1 < mclOutput.I2.out > mclGroups.I2.table
orthomclMclToGroups OG 1 < mclOutput.I3.out > mclGroups.I3.table
# And finally obtain the singleton files (so everything in your goodProteins.fasta file that was not included in the clusters in "mclGroups.I15.table"):
orthomclSingletons goodProteins.fasta mclGroups.I15.table >> singletons.txt
# And you are done! do a final count of the final clusters obtained yielded for your dataset:

#cut -d":" -f1 mclGroups.I15.table | wc -l;
# The number you obtained should match what was written in the log file (in my case called "mcl_ortho.output") for the "mcl" command you ran hopefully in the cluster...
#grep "clusters found" mcl_ortho.output


# OrthoFinder/2.5.2
# 


# For species selected, create a fasta file with protein sequences
# If protein file exists from the assembly files, then download it
# If only cds exists, translate it to protein
# If only gff exists, extract CDS and translate to protein
# If protein coding gbff exist, convert to gff
# If no annotation exists, annotate!

# List all possible genomes and transcriptomes
# Have all possible genomes and transcriptomes, with good quality and as complete as possible protein set

# Do a new phylogenetic analysis, a proper one with checking everything so we can decide on the species needed and so we can order them.

# Of those with transcriptome data, make a transcriptome assembly and use for phylogeny

# Two phylogenies: 1 with multi and uni together and 1 containing only the unicellular organisms - Think about using IQ3 if it goes faster.
# 312 data points exist of Ulvophycea short read RNA samples
# 703 data points exist of Trebouxiophyceae short read RNA samples
# 669 data points exist of Sphaeropleales SRA

# Are some genes evolve faster in some lineages compared to others?

# Now, let's do a phylostratiography analysis

# Using the Dollop program from PHYLIP v3.69
# Turn the gene family output into presence-absence matrix for Dollop input and run Dollop
# Download Phylip
wgt https://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz
tar -zxvf phylip-3.697.tar.gz
cd phylip-3.697/src/
mv Makefile.unx Makefile
make install

# Executables are now in phylip-3.697/exe
# Run Dollop
python get_sp_frequency.py mclGroups.I15.table > dollop_input
sed -i 's/\(\w\{4\}\)\(.*\)/\1      \t\2/g' dollop_input

phylip-3.697/exe/dollop

# Important step here is to transform the output to the one with number of genes gained and lost

# Make phylogenetic tree of the species with species name - One multi only, only uni only and one all together

# Run BUSCO for Chara Braunii
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/427/395/GCA_003427395.1_Cbr_1.0 ../data/genomes

mkdir -p klebsormidium_Scaffolds_v1.0
cd klebsormidium_Scaffolds_v1.0
wget http://www.plantmorphogenesis.bio.titech.ac.jp/~algae_genome_project/klebsormidium/kf_download/120824_klebsormidium_Scaffolds_v1.0.fna
wget http://www.plantmorphogenesis.bio.titech.ac.jp/~algae_genome_project/klebsormidium/kf_download/171026_klebsormidium_v1.1.gff
wget http://www.plantmorphogenesis.bio.titech.ac.jp/~algae_genome_project/klebsormidium/kf_download/160614_klebsormidium_v1.1_AA.fasta
wget http://www.plantmorphogenesis.bio.titech.ac.jp/~algae_genome_project/klebsormidium/kf_download/160614_klebsormidium_v1.1_CDS.fasta
wget http://www.plantmorphogenesis.bio.titech.ac.jp/~algae_genome_project/klebsormidium/kf_download/160614_klebsormidium_v1.1_RNA.fasta

ln -s ../genomes/klebsormidium_Scaffolds_v1.0/120824_klebsormidium_Scaffolds_v1.0.fna
ln -s ../genomes/GCA_003427395.1_Cbr_1.0/GCA_003427395.1_Cbr_1.0_genomic.fna 

mv 120824_klebsormidium_Scaffolds_v1.0.fna KleFla.genomic.fna
mv GCA_003427395.1_Cbr_1.0_genomic.fna ChaBra.genomic.fna

for genome in ../data/fasta_softlinks/KleFla.genomic.fna ../data/fasta_softlinks/ChaBra.genomic.fna
do
echo $genome
dir_name=$(echo $genome | cut -f4 -d '/' | sed 's/.genomic.fna//g')
echo $dir_name
sbatch -J $dir_name busco.sh $genome $dir_name > slurm_job_IDs/busco/${dir_name}.out
done

zcat GCA_003427395.1_Cbr_1.0/GCA_003427395.1_Cbr_1.0_genomic.gff.gz | grep -v "#" | awk '$3=="gene"' | wc
zcat GCA_021605115.1_Astre_guber_v1.0/GCA_021605115.1_Astre_guber_v1.0_genomic.gff.gz | grep -v "#" | awk '$3=="gene"' | wc
zcat GCA_019650175.1_Vafri_1.0/GCA_019650175.1_Vafri_1.0_genomic.gff.gz | grep -v "#" | awk '$3=="gene"' | wc
zcat GCA_000143455.1_v1.0/GCA_000143455.1_v1.0_genomic.gff.gz | grep -v "#" | awk '$3=="gene"' | wc
zcat GCA_019650235.1_Vretifemale_1.0/GCA_019650235.1_Vretifemale_1.0_genomic.gff.gz | grep -v "#" | awk '$3=="gene"' | wc
zcat GCA_905146915.1_Ostreobium_1D_genome/GCA_905146915.1_Ostreobium_1D_genome_genomic.gff.gz | grep -v "#" | awk '$3=="gene"' | wc

Just make a list of all transcriptomes avaialble 
Number of unique species?
Where it is placed in the phylogeny?