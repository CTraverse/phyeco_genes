# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:58:56 2017

@author: chuck
"""



import os, gzip
from sys import argv


script, midas_db_path, sorted_phyeco_map, phyeco_genes, metagenome_path, = argv
#Build dictionary of representative genomes as keys that call phyeco_id and gene_id

phyeco = list(open(sorted_phyeco_map, 'r'))

del phyeco[0] #remove header

phyeco_map = {} #Master dictionary for phyeco genes

for phy in phyeco: #loop through each phyeco entry
    phy_split = phy.strip().split('\t')
    species_name, phyeco_id, gene_id = phy_split[3], phy_split[4], phy_split[0]
    if species_name not in phyeco_map: #If the species name is not already entered into the dictionary, enter it
        phyeco_map[species_name] = {} #New dict entry (key) is species name and is populated with a blank dictionary

    if phyeco_id not in phyeco_map[species_name]:#If the phyeco ID is not already in the dictionary for the species your on, add the phyeco id
        phyeco_map[species_name][phyeco_id] = [] #Phyeco id is now the key for the species dictionary and is populated with a blank list

    phyeco_map[species_name][phyeco_id].append(gene_id) #Append every gene in the phyeco.map.sorted file to its corresponding species and phyeco ID
                                                        #Stored as a list because some species have multiple genes associated with each phyeco ID. So, store them all


#List of phyeco genes we want to look at
phyeco_genes_list = []
with open(phyeco_genes, 'r') as phy:
    for gene in phy:
        phyeco_genes_list.append(gene.strip())


#Build dictionary of genome of interest with gene_id as keys that call scaffold, gene start, gene stop, and strand
##I wrote this as a function to clean up the master loop below. 
def get_features(path_to_features): #Uses the path to a specific species features file as the input

    dict_of_features = {} #Create blank dicitonary
    with gzip.open(path_to_features, 'rb') as genome_features: #Open the .gzipped features file specified in the path
        #Everything under the "with gzip.open..." with an indent will occur while the gzipped file is being read
        for feature in genome_features: #loop through each gene in the features file
            if 'gene_id' in feature: #Skip the first line
                pass
            else:
                feature_split = feature.strip().split('\t')
                dict_of_features[feature_split[0]] = [feature_split[1], feature_split[2], feature_split[3], feature_split[4]] #use the gene name as the dictionary key, and populate with a list of the gene features of interest

    return dict_of_features # Return the populated dictionary to the master loop below


###Path to snps output
snps_path = metagenome_path + "/snps/output/"

snps_genomes_list = []
for filename in [gz for gz in os.listdir(snps_path) if gz.endswith("gz")]:
    snps_genomes_list.append(filename)

#snps_genomes_list = next(os.walk(snps_path))[2] #This function creates a list of all the names of files in the specified path. Using the [2] Ignores directory names in that path
rep_genomes_dir = midas_db_path + "/rep_genomes/" #Path to representatve genomes directory


for snps_genome in snps_genomes_list: ##Iterate through genomes with enough coverage to call snps. Only genomes with sufficient coverage will appear in this directory
    genome_directory = rep_genomes_dir + snps_genome.strip('.snps.gz') #Build path to representative genome.
    features_dict = get_features(genome_directory + '/genome.features.gz') #Grab gene features of current genome. See function above
    new_genome_dir = metagenome_path + "/marker_genes_snps/" + snps_genome.strip('.snps.gz') + "/"
    if not os.path.exists(os.path.dirname(new_genome_dir)):
        os.makedirs(new_genome_dir)

    with gzip.open(genome_directory + '/genome.fna.gz') as check_genome_id: #Need to figure out which gene from the phyeco.map.sorted file corresponds to the rep genome
        description = check_genome_id.readline() #Just grab the first line and then close the file by putting nothing else indented under the "with gzip..."

    for phyeco_gene in phyeco_genes_list: #Loop through each marker gene

        current_gene_out = open(new_genome_dir + snps_genome.strip('.snps.gz') + "." + phyeco_gene, "w" )

        potential_genes_list = phyeco_map[snps_genome.strip('.snps.gz')][phyeco_gene] #Now, you have a specific species with a specific phyeco ID. Grab the genes associated with dicitonary entry
        for gene in potential_genes_list: #Loop through the list of genes associated with that species and phyeco ID
            genome_identifier = gene.split('.peg')[0] #Split by '.peg' because everything before .peg is the genome identifier
            if genome_identifier in description: #Is the specific genome identifier associated with your phyeco gene found in the first line of your representative genome? If so, you found the correct gene
                ref_gene = gene
                break #Break this loop because you found your gene

        gene_data_list = [] #This will contain the snp data for your speces/phyeco id gene
        phyeco_gene_features = features_dict[ref_gene] #Grab the gene features for your gene from the features dictionary that contains all genes and locations for your representative genome
        scaffold, gene_start, gene_end, strand = phyeco_gene_features[0], int(phyeco_gene_features[1]), int(phyeco_gene_features[2]), phyeco_gene_features[3] #Parse out that features information

        with gzip.open(snps_path + snps_genome) as snps_file: #Open your snps file
            for base in snps_file: #loop through each base in the genome
                if scaffold in base: #If you're on the correct scaffold, do the stuff below. Otherwise, go to the next base
                    base_split = base.strip().split('\t') #This is used to split the information based on the tabs and the .strip() gets rid of the \n regular expression at the end of the line
                    current_loc = int(base_split[1]) #Current reference location within scaffold
                    if current_loc >= gene_start and current_loc <= gene_end: #If your current location is greater than or equal to the start of the gene and less than or equal to the end, then you're inside the gene of interest
                        current_gene_out.write(base)
                        gene_data_list.append(base.split) #Apppend that information to the list you made above
                        print base.strip() #Print current line

