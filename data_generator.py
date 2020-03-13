import pandas as pd
import argparse
import sys
import time
import pickle

import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri

parser = argparse.ArgumentParser(description='.')
parser.add_argument('--exp-index', type=int, dest='exp_index', default=5,
                    help='experiment index number')
parser.add_argument('--corr-method', dest='corr_method', default="pearson",
                    help='correlation method ex pearson/spearman')
parser.add_argument('--from-raw', type=bool, dest='fram_raw', default=False,
                    help='boolean generate coexp matrix from raw data')
parser.add_argument('--test', type=bool, dest='test', default=False,
                    help='boolean test mode')
args = parser.parse_args()

print(args)


output_directory = "./data/"
coexp_file = {
	
	# "proteomicsdb": "../ProteomicsDB/coexp_proteomics_df_filtered_proteomicsdb.csv",
	# "ccle": "../ProteomicsDB/coexp_proteomics_df_filtered_CCLE.csv",
	# "proteomics": "../ProteomicsDB/coexp_proteomics_df_filtered_ccle+proteomics_new.csv",	
	# "transcriptome": "../../Data/ARCHS4/human_correlation.rda",
	# "proteomics_normalized": "../ProteomicsDB/coexp_proteomics_df_filtered_ccle+proteomics_normalized.csv"
	"proteomicsdb": "./data/coexp_proteomics_df_filtered_proteomicsdb.csv",
	"ccle": "./data/coexp_proteomics_df_filtered_CCLE.csv",
	"proteomics": "./data/coexp_proteomics_df_filtered_ccle+proteomics_new.csv",	
	"transcriptome": "./data/human_correlation.rda",
	"proteomics_normalized": "./data/coexp_proteomics_df_filtered_ccle+proteomics_normalized.csv"

}

raw_data_file = {
	"proteomicsdb": "./data/proteomics_matrix_filtered_proteomicsdb.csv",
	"ccle": "./data/ProteomicsDB/proteomics_matrix_filtered_CCLE.csv",
}

def read_rda(filename):

	pandas2ri.activate()
	data = robjects.r['load'](filename)
	cc = robjects.r['cc']
	df = pd.DataFrame(cc)
	r('col<-colnames(cc)')
	r('idx<-rownames(cc)')
	df.columns = robjects.r['col']
	df.index = robjects.r['idx']
	
	return df

def filter_common_genes(df1, df2):
	gene_list1 = set(list(df1.index))
	gene_list2 = set(list(df2.index))

	common_gene_list = gene_list1.intersection(gene_list2)

	return df1.loc[common_gene_list, common_gene_list], df2.loc[common_gene_list, common_gene_list]




def generate_coexp_matrix(ratio, from_raw=False):

	print("Generating co-expression matrix...")
	if from_raw == True:
		print("Generating co-expression matrix from raw datasets...")
		proteomicsdb = pd.read_csv(raw_data_file["proteomicsdb"], index_col=0)
		ccle = pd.read_csv(raw_data_file["ccle"], index_col = 0)
		
		if args.exp_index == 8:
			# z-norm before concatenation
			print("Normalization...")
			proteomicsdb_t = proteomicsdb.transpose().copy()
			z_norm_proteomicsdb = ((proteomicsdb_t-proteomicsdb_t.mean())/proteomicsdb_t.std()).transpose()
			
			ccle_t = ccle.transpose().copy()
			z_norm_ccle = ((ccle_t-ccle_t.mean())/ccle_t.std()).transpose()
			
			proteomics = z_norm_proteomicsdb.join(z_norm_ccle, how='inner')

		else:
			print("Concat...")
			proteomics = proteomicsdb.join(ccle, how='inner')

		print("Conducting coexp func...")
		coexp_proteomics = coexp_matrix(proteomics, method=args.corr_method)

		print("Saving...")
		coexp_proteomics.to_csv(coexp_file["proteomics"])
		# coexp_proteomicsdb.to_csv(coexp_file["proteomicsdb"])

		# coexp_ccle= coexp_matrix(ccle, method=args.corr_method)
		# coexp_ccle.to_csv(coexp_file["ccle"])
# proteomics_normalized


	else:
		print("Generating co-expression matrix from pre-processed datasets...")
		# coexp_proteomicsdb = pd.read_csv(coexp_file["proteomicsdb"], index_col=0)
		# coexp_ccle = pd.read_csv(coexp_file["ccle"], index_col=0)
		strt_time = time.time()
		coexp_proteomics = pd.read_csv(coexp_file["proteomics"], index_col=0)
		print(time.time()-strt_time)
	strt_time = time.time()
	coexp_transcriptomics = read_rda(coexp_file["transcriptome"])
	print(time.time()-strt_time)
	# filter common genes
	
	if args.exp_index in [5, 8]:
		coexp_transcriptomics, coexp_proteomics = filter_common_genes(coexp_transcriptomics, coexp_proteomics)
		coexp_df = ratio * coexp_proteomics + (1-ratio) * coexp_transcriptomics
	elif args.exp_index == 2:
		coexp_df = coexp_transcriptomics

	# elif args.exp_index == 8:




	print(coexp_df.head())
	print(coexp_df.shape)
	return coexp_df


def coexp_matrix(df, method="pearson"):
	
	# return df.T.corr(method=method)
	corr = np.corrcoef(df.values)
	result_df = pd.DataFrame(corr)
	result_df.index = df.index
	result_df.columns = df.index
	return result_df



def save(obj, ratio):
	print("Saving co-expression matrix...")
	output_fname = output_directory+"exp{}_{}.pkl".format(args.exp_index, ratio)
	# with open(output_fname, "wb") as f:
	# 	pickle.dump(obj, f)
	obj.to_pickle(output_fname)

	return output_fname


def main():

	for ratio in range(0, 11, 1):
		coexp_df = generate_coexp_matrix(ratio*0.1, args.fram_raw)
		if args.test == True:
			output_fname = save(coexp_df.iloc[:100, :100], ratio)
		else:
			output_fname = save(coexp_df, ratio)

		print("Saved! {}".format(output_fname))

		if args.exp_index != 5:
			break


if __name__ == "__main__":

	main()