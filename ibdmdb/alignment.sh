
# Trim
mkdir metaT_trimmed/
parallel -j 1 --xapply trim_galore --paired -fastqc -j 8 -o ./metaT_trimmed/ ::: ./metaT/*_1.fastq.gz ::: ./metaT/*_2.fastq.gz
mkdir metaG_trimmed/
parallel -j 1 --xapply trim_galore --paired -fastqc -j 8 -o ./metaG_trimmed/ ::: ./metaG/*_1.fastq.gz ::: ./metaG/*_2.fastq.gz
    
cd metaG_trimmed/
ls -1 *_1_val_1.fq.gz | rev | cut -c15- | rev | sort | uniq > ../metaG_samples.txt
cd ../

cd metaT_trimmed/
ls -1 *_1_val_1.fq.gz | rev | cut -c15- | rev | sort | uniq > ../metaT_samples.txt
cd ../


#################
# WoL metaG
#################

## Create output directories
mkdir WoL_metaG
mkdir WoL_metaG/samfiles
mkdir WoL_metaG/bowfiles
# Align fastqfiles to custom database
while read -r sample; do bowtie2 -1 ./metaG_trimmed/${sample}_1_val_1.fq.gz -2 ./metaG_trimmed/${sample}_2_val_2.fq.gz -x ./WoL/WoLr1 -p 8 --no-unal --no-head -S ./WoL_metaG/samfiles/$sample\_Custom_index.sam 2> ./WoL_metaG/bowfiles/$sample\_Custom_index.bow; done <./metaG_samples.txt
# Calculate genomes coverages
python ./calculate_coverages.py	-i ./WoL_metaG/samfiles -o ./WoL_metaG/coverages.txt -d ./metadata.tsv
    
#################
# metaG
#################

## Create output directories
mkdir metaG_aligned
mkdir metaG_aligned/samfiles
mkdir metaG_aligned/bowfiles
# Align fastqfiles to custom database
while read -r sample; do bowtie2 -1 ./metaG_trimmed/${sample}_1_val_1.fq.gz -2 ./metaG_trimmed/${sample}_2_val_2.fq.gz -x ./index_by_genome/whole_genome_index -p 8 --no-unal --no-head -S ./metaG_aligned/samfiles/$sample\_Custom_index.sam 2> ./metaG_aligned/bowfiles/$sample\_Custom_index.bow; done <./metaG_samples.txt
# Calculate genomes coverages
python ./calculate_coverages.py	-i ./metaG_aligned/samfiles -o ./metaG_aligned/coverages.txt -d ./custom_metadata.tsv

woltka classify -i ./metaG_aligned/samfiles -o ./metaG_aligned/samfiles/Custom_index_counts.tsv

# # Creation of a subset of the WoL database that contains only genomes with >2% genome coverage:
# for id in $(cat ./metaG_aligned/gotu_above_2percent.txt); do grep ${id} ./fasta/merged/773strains.fna -A1 >> ./subset/subset.concat.fna; done
# # Build new subset database index
# bowtie2-build -f ./subset/subset.concat.fna ./subset/subset -p 8

#################
# metaT
#################
## Create output directories
mkdir metaT_aligned
mkdir metaT_aligned/samfiles
mkdir metaT_aligned/bowfiles
# Align fastqfiles to custom database
while read -r sample; do bowtie2 -1 ./metaT_trimmed/${sample}_1_val_1.fq.gz -2 ./metaT_trimmed/${sample}_2_val_2.fq.gz -x ./index_by_gene/annotated_genome_index -p 8 --no-unal --no-head -S ./metaT_aligned/samfiles/$sample\_Custom_index.sam 2> ./metaT_aligned/bowfiles/$sample\_Custom_index.bow; done <./metaT_samples.txt
# Counts
woltka classify -i ./metaT_aligned/samfiles -o ./metaT_aligned/samfiles/Custom_index_counts.tsv