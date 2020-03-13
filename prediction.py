import pandas as pd
import os
import argparse
import pickle
import sys
import math
import numpy as np

from sklearn.metrics import auc
from numpy import trapz


parser = argparse.ArgumentParser(description='.')
parser.add_argument('--exp-index', type=int, dest='exp_index', default=5,
                    help='experiment index number')
parser.add_argument('--pheno-cutoff', type=int, dest='pheno_cutoff', default=5,
                    help='minimum number of gene in a gene set')
args = parser.parse_args()

print(args)

inputfile_directory = "./data/"

library_files = {
	"GO_Biological_Process_2018": "../../Data/Gene_set_libraries/GO_Biological_Process_2018.txt",
	"GO_Molecular_Function_2018": "../../Data/Gene_set_libraries/GO_Molecular_Function_2018.txt",
	"ENCODE_Histone_Modifications_2015": "../../Data/Gene_set_libraries/ENCODE_Histone_Modifications_2015.txt",
	"ENCODE_TF_ChIP-seq_2015": "../../Data/Gene_set_libraries/ENCODE_TF_ChIP-seq_2015.txt",
	"GWAS_Catalog_2019": "../../Data/Gene_set_libraries/GWAS_Catalog_2019.txt",
	"KEGG_2019_Human": "../../Data/Gene_set_libraries/KEGG_2019_Human.txt",
	"KEGG_2019_Mouse": "../../Data/Gene_set_libraries/KEGG_2019_Mouse.txt",
	"MGI_Mammalian_Phenotype_2013": "../../Data/Gene_set_libraries/MGI_Mammalian_Phenotype_2013.txt",
	"MGI_Mammalian_Phenotype_2017": "../../Data/Gene_set_libraries/MGI_Mammalian_Phenotype_2017.txt",
	"MGI_Mammalian_Phenotype_Level_3": "../../Data/Gene_set_libraries/MGI_Mammalian_Phenotype_Level_3.txt",
	"MGI_Mammalian_Phenotype_Level_4": "../../Data/Gene_set_libraries/MGI_Mammalian_Phenotype_Level_4.txt",
	"MGI_Mammalian_Phenotype_Level_4_2019": "../../Data/Gene_set_libraries/MGI_Mammalian_Phenotype_Level_4_2019.txt",
	

}

def load(ratio):
	with open(inputfile_directory+"exp{}_{}.pkl".format(args.exp_index, ratio), "rb") as f:
		df = pickle.load(f)
	return df

def fill_na_diagonal(df):
	tmp = df.values
	np.fill_diagonal(tmp, None)
	tmp_df = pd.DataFrame(tmp)
	tmp_df.index = df.index
	tmp_df.columns = df.columns

	return tmp_df


def prediction(df, geneset_gene_dict, gene_geneset_dict):
	
	print("Predicting gene functions...")
	prediction_dict = dict()
	
	for geneset in geneset_gene_dict.keys():

		prediction_geneset = df.loc[:, geneset_gene_dict[geneset]].mean(axis=1)
		prediction_geneset_dict = prediction_geneset.to_dict()
		prediction_dict[geneset] = prediction_geneset_dict
		
	prediction_df = pd.DataFrame.from_dict(prediction_dict)
	
	return prediction_df


def load_library(library_name, gene_list):
	print("Loading libraries...")
	geneset_gene_dict = dict()
	gene_geneset_dict = dict()
	with open(library_files[library_name], "r") as f:
		lines = f.read().strip().split("\n")
		for line in lines:
			splited = line.strip().split("\t")
			genes = [x for x in splited[2:] if x != "" and x in gene_list ]
			if len(genes) >= args.pheno_cutoff:
				geneset_gene_dict[splited[0]] = genes
				for gene in genes:
					if gene in gene_geneset_dict:
						gene_geneset_dict[gene].append(splited[0])
					else:
						gene_geneset_dict[gene] = [splited[0]]
	return geneset_gene_dict, gene_geneset_dict

def sort_df (df):
	df = df.reindex(sorted(df.columns), axis=1)
	return df

def scale_vector(v, start, end):

	scaled_vector = (v - min(v)) / max(v - min(v)) * (end - start) + start

	return scaled_vector


def evaluation(prediction_df, gene_geneset_dict):
	print("Evaluating...")
	result = dict()
	for gene in gene_geneset_dict.keys():
		if len(gene_geneset_dict[gene]) >= args.pheno_cutoff:
			tmp_prediction_df = prediction_df.loc[gene, :]
			tmp_prediction_values = tmp_prediction_df.values
			
			v1 = [x in gene_geneset_dict[gene] for x in tmp_prediction_df.index]
			sorted_idx = sorted(range(len(tmp_prediction_values)), key=lambda k: tmp_prediction_values[k], reverse=True)
			ov1 = [v1[x] for x in sorted_idx]

			cumulative = np.cumsum(ov1)
			scaled_y = scale_vector(cumulative, 0, 1)
			scaled_x = scale_vector(np.arange(0, len(cumulative), 1), 0, 1)
			auc = trapz(scaled_y, scaled_x)
			result[gene] = auc
		# break
	# print(result)
	return result
		# break


def main():
	with open("./result.txt", "w") as f:
		for library in library_files.keys():
			for ratio in range(0, 11, 1):
				coexp_df = load(ratio)
				coexp_df = fill_na_diagonal(coexp_df)

				print(coexp_df.shape)
				gene_list = list(coexp_df.index)
				
				geneset_gene_dict, gene_geneset_dict = load_library(library, gene_list)
				print(len(geneset_gene_dict.keys()), "sets")
				
				prediction_df = prediction(coexp_df, geneset_gene_dict, gene_geneset_dict)
				prediction_df = sort_df(prediction_df)
				auc_dict = evaluation(prediction_df, gene_geneset_dict)
				f.write(",".join(map(str, [library, ratio, np.mean(list(auc_dict.values()))])))
				f.write("\n")
				f.flush()

				with open("./data/auc_exp{}_{}_{}.pkl".format(args.exp_index, ratio, library), "wb") as fpk:
					pickle.dump(auc_dict, fpk)

		# break

if __name__ == "__main__":
	main()