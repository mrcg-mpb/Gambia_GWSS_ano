#!/bin/bash
# This script is used to set up the environment for the IDB project.
# It installs the required packages and sets up the database.
# Usage: ./IDB_script.sh

#check if beagle is installled

if ! command -v eagle &> /dev/null
then
    echo "eagle could not be found, please install it first."
    exit
fi

# install refined
wget https://faculty.washington.edu/browning/refined-ibd/refined-ibd.17Jan20.102.jar
mv refined-ibd.17Jan20.102.jar refined.jar
# add it to the path
echo "refined installed successfully."


INPUT="3R_sweep_region_renamed.vcf.gz"

# Create a map file
echo "Creating genetic map file..."
bcftools query -f '%CHROM\t%POS\t%ID\n' $INPUT > anopheles_genetic_map.txt
./generate_map.py anopheles_genetic_map.txt anopheles_map.txt

# compress the map file
echo "Compressing genetic map file..."
bgzip anopheles_map.txt
bgzip ano_eagle_format.txt  

#-----OPTIONAL STEP-----
# Reformat chrom names vcf to eagle format
echo "Reformatting VCF file..."
echo "3R    3" > chrom_mapping.txt
bcftools annotate --rename-chrs chrom_mapping.txt $INPUT -Oz -o 3R_28,5_renamed_2.vcf.gz
#--END OPTIONAL STEP-----
# Run eagle
echo "Running eagle..."
#eagle   --vcf=3R_28,5_renamed_2.vcf.gz    --geneticMapFile=ano_eagle_format.txt.gz   --outPrefix phased_data

eagle   --vcf=$INPUT   --geneticMapFile=ano_eagle_format.txt.gz   --outPrefix phased_data

# Run refined
echo "Running refined..."
java -jar refined.jar gt=phased_data.vcf.gz  window=40 length=0.05 out=phased_data

echo "Refined completed successfully."

