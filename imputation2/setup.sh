#!/bin/bash

CURRENT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Get imputation reference data and related programs
mkdir "$CURRENT/impute_tools"
cd "$CURRENT/impute_tools"

# Get plink 1.9
wget -P plink http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200219.zip
cd "plink"
unzip plink_linux_x86_64_20200219.zip
rm plink_linux_x86_64_20200219.zip
cd "$CURRENT/impute_tools"

# Get impute2
wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz
gunzip impute_v2.3.2_x86_64_static.tgz
tar -xvf impute_v2.3.2_x86_64_static.tar
rm impute_v2.3.2_x86_64_static.tar

# Get shapeit
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
tar zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
rm shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
mv bin shapeit

# Get gtools
wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz
tar zxvf gtool_v0.7.5_x86_64.tgz
rm gtool_v0.7.5_x86_64.tgz

# Get the reference from 1000genomes
wget https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
gunzip ALL_1000G_phase1integrated_v3_impute.tgz
tar xf ALL_1000G_phase1integrated_v3_impute.tar
rm ALL_1000G_phase1integrated_v3_impute.tar
wget https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_annotated_legends.tgz
gunzip ALL_1000G_phase1integrated_v3_annotated_legends.tgz
tar xf ALL_1000G_phase1integrated_v3_annotated_legends.tar
rm ALL_1000G_phase1integrated_v3_annotated_legends.tar
mv ALL_1000G_phase1integrated_v3_annotated_legends/* ALL_1000G_phase1integrated_v3_impute/
rmdir ALL_1000G_phase1integrated_v3_annotated_legends

# Link to the X-chr
cd ALL_1000G_phase1integrated_v3_impute
ln -s genetic_map_chrX_nonPAR_combined_b37.txt genetic_map_chrX_combined_b37.txt
ln -s ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz ALL_1000G_phase1integrated_v3_chrX_impute.hap.gz
ln -s ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz ALL_1000G_phase1integrated_v3_chrX_impute.legend.gz

# Modify path in config.yml
replace=$(echo $CURRENT | sed 's_/_\\/_g')
sed -i -e "s/\/path\/to/$replace/" "$CURRENT/config.yml"

echo -e "\033[1mDone! Check path in $CURRENT/config.yml and installed tools in $CURRENT/impute_tools before running the imputation!\033[0m"
