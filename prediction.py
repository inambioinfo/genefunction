import pandas as pd
import os
import argparse
import pickle
import sys
import numpy as np

parser = argparse.ArgumentParser(description='.')
parser.add_argument('--exp-index', dest='exp_index', default=5,
                    help='experiment index number')
parser.add_argument('--pheno-cutoff', dest='pheno_cutoff', default=5,
                    help='minimum number of gene in a gene set')
args = parser.parse_args()

print(args)

inputfile_directory = "./data/"

library_files = {
	"GO_BP": "../../Data/Gene_set_libraries/GO_Biological_Process_2018.txt",


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
			for gene in splited[2:]:
				if gene in gene_geneset_dict:
					gene_geneset_dict[gene].append(splited[0])
				else:
					gene_geneset_dict[gene] = [splited[0]]
	return geneset_gene_dict, gene_geneset_dict


def evaluation():


	for gene in gene_geneset_dict.keys():
		if gene_geneset_dict[gene] >= args.phen_count:
			
			
#scaling the prediction because is returns the best results like this
spred = scale(prediction)

auc = c()
for (current_gene in sgenes) {
    geneprob = spred[current_gene,]
    v1 = colnames(spred) %in% reversegmt[[current_gene]] 

    oo = rev(order(geneprob))
    ov1 = v1[oo]

    cumulative = cumsum(ov1)
    scaled_y = scale_vector(cumulative, 0, 1)
    scaled_x= scale_vector(1:length(cumulative), 0, 1)
    
    auc[current_gene] = trapz(scaled_x, scaled_y)
    #plot((1:length(cumulative))/length(cumulative),cumulative/max(cumulative), type="l", main=paste(current_gene, round(auc[current_gene], digits=2)))
    #abline(0,1, col="blue")
}

def main():
	for ratio in range(0, 11, 1):
		coexp_df = load(ratio)
		coexp_df = fill_na_diagonal(coexp_df)
		gene_list = list(coexp_df.index)
		
		geneset_gene_dict, gene_geneset_dict = load_library("GO_BP", gene_list)
		print(len(geneset_gene_dict.keys()), "sets")
		
		prediction_df = prediction(coexp_df, geneset_gene_dict, gene_geneset_dict)

		evaluation()
		break

if __name__ == "__main__":
	main()