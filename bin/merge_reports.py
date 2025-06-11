#!/usr/bin/env python

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

import pandas as pd
import glob
import re
import os
import sys
import errno
import argparse


def parse_args(args=None):
	Description = "Merges svbenchmark or truvari bench reports from multiple samples"
	Epilog = "Example usage: python merge_reports.py file1 file2 file3 -o merged_table.csv -b truvari/svbenchmark/wittyer/happy/sompy -v snv/indel -a germline/somatic "

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument("inputs", nargs="+", help="List of files to merge")
	parser.add_argument("--output", "-o", required=True, help="Output file")
	parser.add_argument("--bench", "-b", required=True, help="svbenchmark/truvari/happy/sompy")
	parser.add_argument("--vartype", "-v", required=True, help="Variant type: snv,indel,sv,small")
	parser.add_argument("--analysis", "-a", required=True, help="Analysis type: germline,somatic")

	return parser.parse_args(args)


## SVanalyzer results
def get_svbenchmark_results(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Define regular expressions to extract the values
	DTP_pattern = re.compile(r'Number of detected true variants \(.*\): (\d+)')
	FN_pattern = re.compile(r'Number of undetected true variants \(.*\): (\d+)')
	PTP_pattern = re.compile(r'Number of predictions that are true \(.*\): (\d+)')
	FP_pattern = re.compile(r'Number of false positives \(.*\): (\d+)')
	recall_pattern = re.compile(r'Recall\s*\(.*?\):\s*(\d+(?:\.\d+)?)%')
	precision_pattern = re.compile(r'Precision\s*\(.*?\):\s*(\d+(?:\.\d+)?)%')
	f1_pattern = re.compile(r'F1 \(.*\): ([\d\.]+(?:e[+-]?\d+)?)')

	# Iterate over each table file
	for file in file_paths:
		# Read the table into a DataFrame
		filename = os.path.basename(file)
		with open(file, 'r') as f:
			text = f.read()

		# Search for matches in the text
		DTP_match = DTP_pattern.search(text)
		FN_match = FN_pattern.search(text)
		PTP_match = PTP_pattern.search(text)
		FP_match = FP_pattern.search(text)
		recall_match = recall_pattern.search(text)
		precision_match = precision_pattern.search(text)
		f1_match = f1_pattern.search(text)

		# Initialize a dictionary to store the data
		data = {
			'Tool': [filename.split(".")[0]],
			'File': filename,
			'TP_base': [int(DTP_match.group(1)) if DTP_match else 'NA'],
			'FP': [int(FP_match.group(1)) if FP_match else 'NA'],
			'TP_comp': [int(DTP_match.group(1)) if DTP_match else 'NA'],
			'FN': [int(FN_match.group(1)) if FN_match else 'NA'],
			'Recall': [float(recall_match.group(1))/100 if recall_match else 'NA'],
			'Precision': [float(precision_match.group(1))/100 if precision_match else 'NA'],
			'F1': [float(f1_match.group(1)) if f1_match else 'NA']}

		df = pd.DataFrame(data)

		merged_df = pd.concat([merged_df, df], ignore_index=True)

	return merged_df

## Truvari results

def get_truvari_results(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
	# Read the json into a DataFrame
		filename = os.path.basename(file)
		with open(file, 'r') as f:
			data = pd.read_json(f)

			relevant_data = {
				"Tool": filename.split(".")[0],
				"File": filename,
				"TP_base": int(data["TP-base"].iloc[0]),
				"TP_comp": int(data["TP-comp"].iloc[0]),
				"FP": int(data["FP"].iloc[0]),
				"FN": int(data["FN"].iloc[0]),
				"Precision": float(data["precision"].iloc[0]),
				"Recall": float(data["recall"].iloc[0]),
				"F1": float(data["f1"].iloc[0])}

		df = pd.DataFrame([relevant_data])
		merged_df = pd.concat([merged_df, df], ignore_index=True)

	return merged_df

def get_wittyer_results(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	for file in file_paths:
	# Read the json into a DataFrame
		filename = os.path.basename(file)
		with open(file, 'r') as f:
			data = pd.read_json(f)

			relevant_data = []
			for sample in data['PerSampleStats']:
				for stats in sample['OverallStats']:
					relevant_data.append({
						"Tool": filename.split(".")[0],
						"File": filename,
						"StatsType": stats["StatsType"],
						"TP_base": int(stats["TruthTpCount"]) if pd.notna(stats["TruthTpCount"]) else 0,
						"TP_comp": int(stats["QueryTpCount"]) if pd.notna(stats["QueryTpCount"]) else 0,
						"FP": int(stats["QueryFpCount"]) if pd.notna(stats["QueryFpCount"]) else 0,
						"FN": int(stats["TruthFnCount"]) if pd.notna(stats["TruthFnCount"]) else 0,
						"Precision": float(stats["Precision"]) if pd.notna(stats["Precision"]) else float('nan'),
						"Recall": float(stats["Recall"]) if pd.notna(stats["Recall"]) else float('nan'),
						"F1": float(stats["Fscore"]) if pd.notna(stats["Fscore"]) else float('nan')
					})

		df = pd.DataFrame(relevant_data)
		merged_df = pd.concat([merged_df, df], ignore_index=True)

	return merged_df

def get_rtgtools_results(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
		filename = os.path.basename(file)

		with open(file, 'r') as f:
			lines = f.readlines()

		# Extract header
		header = lines[0].strip().split()

		# Extract data
		data = []
		for line in lines[2:]:
			data.append(line.strip().split())

		# Create DataFrame
		df = pd.DataFrame(data, columns=header)
		df['Tool'] = filename.split(".")[0]
		df['File'] = filename
		df_redesigned = df[['Tool', 'File', 'Threshold','True-pos-baseline','True-pos-call','False-pos','False-neg','Precision','Sensitivity','F-measure']]
		df_redesigned.columns = ['Tool','File', 'Threshold','TP_base','TP_call','FP','FN','Precision','Recall','F1']
		# Convert relevant columns to integers, handling potential NaN values
		int_columns = ['TP_base', 'FN', 'TP_call', 'FP']
		float_columns = ['Recall','Precision','F1']
		df_redesigned[int_columns] = df_redesigned[int_columns].fillna(0).astype(int)
		df_redesigned[float_columns] = df_redesigned[float_columns].fillna(0).astype(float)

		merged_df = pd.concat([merged_df, df_redesigned], ignore_index=True)

	return merged_df

def get_happy_results(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
		filename = os.path.basename(file)

		df = pd.read_csv(file)

		df['Tool'] = filename.split(".")[0]
		df['File'] = filename

		df_redesigned = df[['Tool', 'File', 'Type','Filter','TRUTH.TOTAL','TRUTH.TP','TRUTH.FN','QUERY.TOTAL','QUERY.FP','QUERY.UNK','FP.gt','FP.al','METRIC.Recall','METRIC.Precision','METRIC.Frac_NA','METRIC.F1_Score','TRUTH.TOTAL.TiTv_ratio','QUERY.TOTAL.TiTv_ratio','TRUTH.TOTAL.het_hom_ratio','QUERY.TOTAL.het_hom_ratio']]
		df_redesigned.columns = ['Tool', 'File', 'Type','Filter','TP_base','TP','FN','TP_call','FP','UNK','FP_gt','FP_al','Recall','Precision','Frac_NA','F1','TRUTH_TiTv_ratio','QUERY_TiTv_ratio','TRUTH_het_hom_ratio','QUERY_het_hom_ratio']

		# Convert relevant columns to integers, handling potential NaN values
		int_columns = ['TP_base', 'TP', 'FN', 'TP_call', 'FP', 'UNK', 'FP_gt', 'FP_al']
		float_columns = ['Recall','Precision','Frac_NA','F1','TRUTH_TiTv_ratio','QUERY_TiTv_ratio','TRUTH_het_hom_ratio','QUERY_het_hom_ratio']
		df_redesigned[int_columns] = df_redesigned[int_columns].fillna(0).astype(int)
		df_redesigned[float_columns] = df_redesigned[float_columns].fillna(0).astype(float)

		# Concatenate with the merged DataFrame
		merged_df = pd.concat([merged_df, df_redesigned], ignore_index=True)

	return merged_df

def get_intersect_results(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
		filename = os.path.basename(file)

		df = pd.read_csv(file)

		df['Tool'] = filename.split(".")[0]
		df['File'] = filename

		# Convert relevant columns to integers, handling potential NaN values
		int_columns = ['TP_base','TP_comp', 'FN', 'FP']
		float_columns = ['Recall','Precision','F1']
		df[int_columns] = df[int_columns].fillna(0).astype(int)
		df[float_columns] = df[float_columns].fillna(0).astype(float)

		# Concatenate with the merged DataFrame
		merged_df = pd.concat([merged_df, df], ignore_index=True)

	return merged_df

def get_sompy_results(file_paths, vartype):
# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
		filename = os.path.basename(file)

		df = pd.read_csv(file)

		df['Tool'] = filename.split(".")[0]
		df['File'] = filename

		df_redesigned = df[['Tool','File', 'type','total.truth','tp','fn','total.query','fp','unk','recall','precision','recall_lower','recall_upper','recall2','precision_lower','precision_upper','na','ambiguous','fp.region.size','fp.rate']]
		df_redesigned.columns = ['Tool', 'File', 'Type','TP_base','TP','FN','TP_call','FP','UNK','Recall','Precision','recall_lower','recall_upper','recall2','precision_lower','precision_upper','na','ambiguous','fp.region.size','F1']
		# Convert relevant columns to integers, handling potential NaN values
		int_columns = ['TP_base', 'TP', 'FN', 'TP_call', 'FP', 'UNK']
		float_columns = ['Recall','Precision','recall_lower','recall_upper','recall2','precision_lower','precision_upper','na','ambiguous','fp.region.size','F1']
		df_redesigned[int_columns] = df_redesigned[int_columns].fillna(0).astype(int)
		df_redesigned[float_columns] = df_redesigned[float_columns].fillna(0).astype(float)

		merged_df = pd.concat([merged_df, df_redesigned], ignore_index=True)

	if vartype == "snv":
		merged_df1 = merged_df[merged_df["Type"] == 'SNVs']
	elif vartype == "indel":
		merged_df1 = merged_df[merged_df["Type"] == "indels"]
	else:
		merged_df1 = merged_df[merged_df["Type"] == "records"]

	if vartype == "snv":
		merged_df2 = merged_df[merged_df["Type"].str.contains(r'SNVs.')]
	elif vartype == "indel":
		merged_df2 = merged_df[merged_df["Type"].str.contains(r"indels.")]
	else:
		merged_df2 = merged_df[merged_df["Type"].str.contains(r"records.")]

	return merged_df1,merged_df2

def main(args=None):
	args = parse_args(args)

	if args.bench == "truvari":
		summ_table = get_truvari_results(args.inputs)

	elif args.bench == "svbenchmark":
		summ_table = get_svbenchmark_results(args.inputs)

	elif args.bench == "wittyer":
		summ_table = get_wittyer_results(args.inputs)

	elif args.bench == "rtgtools":
		summ_table = get_rtgtools_results(args.inputs)

	elif args.bench == "happy":
		summ_table = get_happy_results(args.inputs)

	elif args.bench == "intersect":
		summ_table = get_intersect_results(args.inputs)

	elif args.bench == "sompy":
		summ_table,summ_table2 = get_sompy_results(args.inputs,args.vartype)
		summ_table2.reset_index(drop=True, inplace=True)
		summ_table2.to_csv(args.output + ".regions.csv", index=False)

	else:
		raise ValueError('only results from intersect, truvari, svbenchmark, wittyer, rtgtools, happy or sompy tools can be merged')

	summ_table.reset_index(drop=True, inplace=True)
	summ_table.to_csv(args.output + ".summary.csv", index=False)

if __name__ == "__main__":
	sys.exit(main())
