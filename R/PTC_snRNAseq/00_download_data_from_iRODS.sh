############### THYROID TUMOUR - paed PTC - Y24 snRNAseq - cellranger302 ###############
#irods iget -r /seq/illumina/runs/46/46320/cellranger/cellranger302_count_46320_SB_Thy_R13236839_GRCh38-1_2_0/ /lustre/scratch117/casm/team274/mt22/Data/thyroid_10X/cellranger302_GRCh38-1_2_0/cellranger302_count_46320_SB_Thy_R13236839_GRCh38-1_2_0/

#irods iget -r /seq/illumina/runs/46/46320/cellranger/cellranger302_count_46320_SB_Thy_R13236840_GRCh38-1_2_0/ /lustre/scratch117/casm/team274/mt22/Data/thyroid_10X/cellranger302_GRCh38-1_2_0/cellranger302_count_46320_SB_Thy_R13236840_GRCh38-1_2_0/

#irods iget -r /seq/illumina/runs/46/46320/cellranger/cellranger302_count_46320_SB_Thy_R13236841_GRCh38-1_2_0/ /lustre/scratch117/casm/team274/mt22/Data/thyroid_10X/cellranger302_GRCh38-1_2_0/cellranger302_count_46320_SB_Thy_R13236841_GRCh38-1_2_0/

#irods iget -r /seq/illumina/runs/46/46320/cellranger/cellranger302_count_46320_SB_Thy_R13236842_GRCh38-1_2_0/ /lustre/scratch117/casm/team274/mt22/Data/thyroid_10X/cellranger302_GRCh38-1_2_0/cellranger302_count_46320_SB_Thy_R13236842_GRCh38-1_2_0/

#irods iget -r /seq/illumina/runs/46/46320/cellranger/cellranger302_count_46320_SB_Thy_R13236843_GRCh38-1_2_0/ /lustre/scratch117/casm/team274/mt22/Data/thyroid_10X/cellranger302_GRCh38-1_2_0/cellranger302_count_46320_SB_Thy_R13236843_GRCh38-1_2_0/

#irods iget -r /seq/illumina/runs/46/46320/cellranger/cellranger302_count_46320_SB_Thy_R13236844_GRCh38-1_2_0/ /lustre/scratch117/casm/team274/mt22/Data/thyroid_10X/cellranger302_GRCh38-1_2_0/cellranger302_count_46320_SB_Thy_R13236844_GRCh38-1_2_0/





############### THYROID TUMOUR - paed PTC - snRNAseq - cellranger710 ###############
samples=(

# Donor Y24 
'/seq/illumina/runs/46/46320/cellranger/cellranger710_count_46320_SB_Thy_R13236844_GRCh38-2020-A'
'/seq/illumina/runs/46/46320/cellranger/cellranger710_count_46320_SB_Thy_R13236843_GRCh38-2020-A'
'/seq/illumina/runs/46/46320/cellranger/cellranger710_count_46320_SB_Thy_R13236842_GRCh38-2020-A'
'/seq/illumina/runs/46/46320/cellranger/cellranger710_count_46320_SB_Thy_R13236841_GRCh38-2020-A'
'/seq/illumina/runs/46/46320/cellranger/cellranger710_count_46320_SB_Thy_R13236840_GRCh38-2020-A'
'/seq/illumina/runs/46/46320/cellranger/cellranger710_count_46320_SB_Thy_R13236839_GRCh38-2020-A'

# Donor Y46 
'/seq/illumina/runs/48/48713/cellranger/cellranger710_count_48713_CG_SB_NB14695467_GRCh38-2020-A'
'/seq/illumina/runs/48/48713/cellranger/cellranger710_count_48713_CG_SB_NB14695466_GRCh38-2020-A'
'/seq/illumina/runs/48/48853/cellranger/cellranger710_count_48853_CG_SB_NB14664105_GRCh38-2020-A'
'/seq/illumina/runs/48/48853/cellranger/cellranger710_count_48853_CG_SB_NB14664106_GRCh38-2020-A'
'/seq/illumina/runs/48/48853/cellranger/cellranger710_count_48853_CG_SB_NB14664107_GRCh38-2020-A'
'/seq/illumina/runs/48/48853/cellranger/cellranger710_count_48853_CG_SB_NB14664108_GRCh38-2020-A'
'/seq/illumina/runs/48/48853/cellranger/cellranger710_count_48853_CG_SB_NB14664109_GRCh38-2020-A'
'/seq/illumina/runs/48/48853/cellranger/cellranger710_count_48853_CG_SB_NB14664110_GRCh38-2020-A'
)

outDir=~/0_Projects_git_repo/FetalThyroidAtlas/Data/thyroid_10X/
mkdir -p $outDir

for irod_path in ${samples[@]}; do
   sampleName=$(basename $irod_path)
   mkdir -p $outDir/$sampleName
    
    irods iget -r $irod_path/filtered_feature_bc_matrix $outDir/$sampleName/filtered_feature_bc_matrix
    irods iget -r $irod_path/raw_feature_bc_matrix $outDir/$sampleName/raw_feature_bc_matrix
    
#     irods iget -r $irod_path/possorted_genome_bam.bam $outDir/$sampleName/possorted_genome_bam.bam
#     irods iget -r $irod_path/possorted_genome_bam.bam.bai $outDir/$sampleName/possorted_genome_bam.bam.bai


done