#!/usr/bin/env python

import os
import argparse
import pysam
import pandas as pd

# Parse command line options of this script
parser = argparse.ArgumentParser()
parser.add_argument("--vcf", "-v",required=True, help="Genotype VCF for sample (gzipped with index).")
parser.add_argument("--input","-i", required=True, help="Directory containing BAM file with mapped RNA-seq reads")
args = parser.parse_args()

# Create output dyrectory 
os.makedirs("./Output_OxPhos/", exist_ok=True)
        
# Function that selects and annotates the reads from the bam file tha matches with input variants      
def read_counter(bam,vcf,output_name):
    
    # Create the lists that will define the annotations of the reads in the output dataframe
    Chrom = list()            
    Pos = list()
    SNP = list()
    REF = list()
    Nucleotide = list()
    Rname = list()
    Output = list()
    Gene_ID = list()
    Gene_Symbol = list()
    
    # Iterate over each variant in the VCF file
    for variant in vcf:
    
        # Change the definition of chromosomes to match the format of the files used     
        #chrom = 'chr' + variant.chrom 
        chrom = variant.chrom
        pos=variant.pos
        
        # Get the reads that map to this variant adapted to different formats to define chromosomes
        try:
            reads = bam.fetch(chrom, pos-1, pos)
        except:
            if chrom=="MT":
                chrom="chrM"
            elif chrom!="MT":
                chrom="chr"+chrom
            reads = bam.fetch(chrom, pos-1, pos)
        
        # Iterate over each read in the bam file
        for read in reads:
                        
            # Select only matching reads to annotate them
            if pos -1 in read.get_reference_positions(full_length=True):
                Nucleotide.append(read.query_sequence[read.get_reference_positions(full_length=True).index(pos-1)])
                Pos.append(variant.pos)
                SNP.append(variant.alts[0])
                REF.append(variant.ref[0])
                Rname.append(read.query_name)
                Chrom.append(variant.chrom)
                # Define the genotype harboring the read
                if read.query_sequence[read.get_reference_positions(full_length=True).index(pos-1)]==variant.alts[0]:
                    Output.append(1)
                elif read.query_sequence[read.get_reference_positions(full_length=True).index(pos-1)]==variant.ref[0]:
                    Output.append(0)
                else:
                    Output.append(2)
    
    # Close bam file       
    bam.close()
        
    # Create a dataframe with annotated matching reads and write them in a txt file
    dict={'Chrom':Chrom,'Pos':Pos, 'Ref':REF, 'Seq':Nucleotide, 'Output':Output, 'Rname':Rname}
    data=pd.DataFrame(dict)
    data.to_csv(output_name, sep='\t', index=False)

# Iterate through the input directory to select the bam files.
for filename in os.listdir(args.input):
    if filename.endswith(".bam"):
    
    	# Create the output name for the txt file to be written to the output directory
        output_name='./Output_OxPhos/Cell_'+filename.split(".")[0]+'.output.txt'
        print(filename)
        
        # Load the vcf/vcf.gz file
        if args.vcf.endswith(".gz"):
            with gzip.open(args.vcf, 'rb') as ungz_vcf:
                vcf=pysam.VariantFile(ungz_vcf)
        elif args.vcf.endswith(".vcf"):
            vcf=pysam.VariantFile(args.vcf)
                 # Load the bam file from input directory
        if args.input.endswith("/"):
            bam = pysam.AlignmentFile(args.input+filename, 'rb')
        elif args.input.endswith("/") == False:
            bam = pysam.AlignmentFile(args.input+"/"+filename, 'rb')
        
        # We retrieve all reads that define the genotype
        read_counter(bam,vcf,output_name)

# Close vcf file 
vcf.close()
 
