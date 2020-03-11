import pandas as pd
import argparse
import sys
import time
import pickle



import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri

parser = argparse.ArgumentParser(description='.')
parser.add_argument('--exp-index', dest='exp_index', default=5,
                    help='experiment index number')
parser.add_argument('--corr-method', dest='corr_method', default="pearson",
                    help='correlation method ex pearson/spearman')
parser.add_argument('--from-raw', dest='fram_raw', default=False,
                    help='boolean generate coexp matrix from raw data')

args = parser.parse_args()

print(args)


output_directory = "./data/"
coexp_file = {
	
	"proteomicsdb": "../ProteomicsDB/coexp_proteomics_df_filtered_proteomicsdb.csv",
	"ccle": "../ProteomicsDB/coexp_proteomics_df_filtered_CCLE.csv",
	"proteomics": "../ProteomicsDB/coexp_proteomics_df_filtered_ccle+proteomics.csv",	
	"transcriptome": "../../Data/ARCHS4/human_correlation.rda",

}

raw_data_file = {
	"proteomicsdb": "../ProteomicsDB/proteomics_matrix_filtered_2020_02_21.tsv",
	"ccle": "../ProteomicsDB/proteomics_matrix_filtered_CCLE.tsv",
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
		proteomicsdb = pd.read_csv(raw_data_file["proteomicsdb"], sep="\t", index_col=0)
		coexp_proteomicsdb = coexp_matrix(proteomicsdb, method=args.corr_method)
		coexp_proteomicsdb.to_csv(coexp_file["proteomicsdb"])

		ccle = pd.read_csv(raw_data_file["ccle"], sep="\t", index_col = 0)
		coexp_ccle= coexp_matrix(ccle, method=args.corr_method)
		coexp_ccle.to_csv(coexp_file["ccle"])
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
	coexp_transcriptomics, coexp_proteomics = filter_common_genes(coexp_transcriptomics, coexp_proteomics)
	if args.exp_index == 5:
		coexp_df = ratio * coexp_proteomics + (1-ratio) * coexp_transcriptomics

	print(coexp_df.head())
	print(coexp_df.shape)
	return coexp_df


def coexp_matrix(df, method="pearson"):
	
	return df.T.corr(method=method)



def save(obj, ratio):
	print("Saving co-expression matrix...")
	output_fname = output_directory+"exp{}_{}.pkl".format(args.exp_index, ratio)
	with open(output_fname, "wb") as f:
		pickle.dump(obj, f)

	return output_fname


def main():
	for ratio in range(0, 11, 1):
		coexp_df = generate_coexp_matrix(ratio*0.1, args.fram_raw)
		output_fname = save(coexp_df, ratio)

		print("Saved! {}".format(output_fname))


if __name__ == "__main__":

	main()