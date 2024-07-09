#!/usr/bin/env python

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
def get_svbenchmark_resuls(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Define regular expressions to extract the values
	DTP_pattern = re.compile(r'Number of detected true variants \(.*\): (\d+)')
	FN_pattern = re.compile(r'Number of undetected true variants \(.*\): (\d+)')
	PTP_pattern = re.compile(r'Number of predictions that are true \(.*\): (\d+)')
	FP_pattern = re.compile(r'Number of false positives \(.*\): (\d+)')
	recall_pattern = re.compile(r'Recall \(.*\): (\d+\.\d+)%')
	precision_pattern = re.compile(r'Precision \(.*\): (\d+\.\d+)%')
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
			'TP_base': [DTP_match.group(1) if DTP_match else 'NA'],
			'FP': [FP_match.group(1) if FP_match else 'NA'],
			'TP_comp': [DTP_match.group(1) if DTP_match else 'NA'],
			'FN': [FN_match.group(1) if FN_match else 'NA'],
			'Recall': [float(recall_match.group(1))/100 if recall_match else 'NA'],
			'Precision': [float(precision_match.group(1))/100 if precision_match else 'NA'],
			'F1': [float(f1_match.group(1)) if f1_match else 'NA']}

		df = pd.DataFrame(data)

		merged_df = pd.concat([merged_df, df])

	return merged_df

## Truvari results

def get_truvari_resuls(file_paths):
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
				"TP_base": data["TP-base"].iloc[0],
				"TP_comp": data["TP-comp"].iloc[0],
				"FP": data["FP"].iloc[0],
				"FN": data["FN"].iloc[0],
				"Precision": data["precision"].iloc[0],
				"Recall": data["recall"].iloc[0],
				"F1": data["f1"].iloc[0]}

		df = pd.DataFrame([relevant_data])
		merged_df = pd.concat([merged_df, df])

	return merged_df

def get_wittyer_resuls(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
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
						"StatsType": stats["StatsType"],
						"TP_base": stats["TruthTpCount"],
						"TP_comp": stats["QueryTpCount"],
						"FP": stats["QueryFpCount"],
						"FN": stats["TruthFnCount"],
						"Precision": stats["Precision"],
						"Recall": stats["Recall"],
						"F1": stats["Fscore"]}
					)

		df = pd.DataFrame(relevant_data)
		merged_df = pd.concat([merged_df, df])

	return merged_df

def get_rtgtools_resuls(file_paths):
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
		df_redesigned = df[['Tool', 'Threshold','True-pos-baseline','True-pos-call','False-pos','False-neg','Precision','Sensitivity','F-measure']]
		df_redesigned.columns = ['Tool', 'Threshold','TP_base','TP_call','FP','FN','Precision','Recall','F1']

		merged_df = pd.concat([merged_df, df_redesigned])
	return merged_df

def get_happy_resuls(file_paths):
	# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
		filename = os.path.basename(file)

		df = pd.read_csv(file)

		df['Tool'] = filename.split(".")[0]

		df_redesigned = df[['Tool', 'Type','Filter','TRUTH.TOTAL','TRUTH.TP','TRUTH.FN','QUERY.TOTAL','QUERY.FP','QUERY.UNK','FP.gt','FP.al','METRIC.Recall','METRIC.Precision','METRIC.Frac_NA','METRIC.F1_Score','TRUTH.TOTAL.TiTv_ratio','QUERY.TOTAL.TiTv_ratio','TRUTH.TOTAL.het_hom_ratio','QUERY.TOTAL.het_hom_ratio']]
		df_redesigned.columns = ['Tool', 'Type','Filter','TP_base','TP','FN','TP_call','FP','UNK','FP_gt','FP_al','Recall','Precision','Frac_NA','F1','TRUTH_TiTv_ratio','QUERY_TiTv_ratio','TRUTH_het_hom_ratio','QUERY_het_hom_ratio']

		merged_df = pd.concat([merged_df, df_redesigned])

	return merged_df

def get_sompy_resuls(file_paths, vartype):
# Initialize an empty DataFrame to store the merged data
	merged_df = pd.DataFrame()

	# Iterate over each table file
	for file in file_paths:
		filename = os.path.basename(file)

		df = pd.read_csv(file)

		df['Tool'] = filename.split(".")[0]
		df_redesigned = df[['Tool','type','total.truth','tp','fn','total.query','fp','unk','recall','precision','recall_lower','recall_upper','recall2','precision_lower','precision_upper','na','ambiguous','fp.region.size','fp.rate']]
		df_redesigned.columns = ['Tool','Type','TP_base','TP','FN','TP_call','FP','UNK','Recall','Precision','recall_lower','recall_upper','recall2','precision_lower','precision_upper','na','ambiguous','fp.region.size','fp.rate']

		merged_df = pd.concat([merged_df, df_redesigned])

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

	#check if the files are from svanalyzer or truvari

	if args.analysis == "germline":
		if args.bench == "truvari":
			summ_table = get_truvari_resuls(args.inputs)

		elif args.bench == "svbenchmark":
			summ_table = get_svbenchmark_resuls(args.inputs)

		elif args.bench == "wittyer":
			summ_table = get_wittyer_resuls(args.inputs)

		elif args.bench == "rtgtools":
			summ_table = get_rtgtools_resuls(args.inputs)

		elif args.bench == "happy":
			summ_table = get_happy_resuls(args.inputs)
		else:
			raise ValueError('Only truvari/svbenchmark/wittyer/rtgtools/happy results can be merged for germline analysis!!')

		summ_table.reset_index(drop=True, inplace=True)
		summ_table.to_csv(args.output + ".summary.csv", index=False)

	elif args.analysis == "somatic":
		if args.bench == "sompy":
			summ_table,summ_table2 = get_sompy_resuls(args.inputs,args.vartype)
		else:
			raise ValueError('Only sompy results can be merged for somatic analysis!!')

		## reset index
		summ_table.reset_index(drop=True, inplace=True)
		summ_table2.reset_index(drop=True, inplace=True)

		# Save the merged DataFrame to a new CSV file
		summ_table.to_csv(args.output + ".summary.csv", index=False)
		summ_table2.to_csv(args.output + ".regions.csv", index=False)
	else:
		raise ValueError('Analysis must be germline or somatic')

if __name__ == "__main__":
	sys.exit(main())
