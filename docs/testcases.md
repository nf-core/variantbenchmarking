# nf-core/variantbenchmarking: Test cases

This pipeline is able to benchmark various type of analysis. Below, explanations of some use cases given in [tests](../conf/tests/) will be explained.

## Case 1 : Germline Structural Variant benchmarking

### Config file

- [test config](../conf/tests/germline_sv.config)

### Analysis

- We are using 3 different publicly available structural variant calls for HG002 to benchmark against Genome in a Bottle HG002 SV analysis. Only chromosome 21 will be used for the analysis.
- Preprocessing includes normalization enabling left alignment of the variants, splitting multi allelic variants, and deduplication.
- Size filtering by >30 bases applied to test files to ensure equal SV size in benchmarking.
- Filtering out BND and TRA type of variants
- Truvari, SVbenchmark and Wittyer are used as benchmarking method. Tools spesific parameters are given in the corresponding samplesheet.

### Results

_Truvari_

| Tool  | File                           | TP_base | TP_comp | FP  | FN  | Precision | Recall  | F1      |
| ----- | ------------------------------ | ------- | ------- | --- | --- | --------- | ------- | ------- |
| test1 | test1.HG002.manta.summary.json | 254     | 254     | 75  | 584 | 0.77204   | 0.30310 | 0.43530 |
| test2 | test2.HG002.lumpy.summary.json | 38      | 38      | 20  | 800 | 0.65517   | 0.04535 | 0.08482 |
| test3 | test3.HG002.delly.summary.json | 50      | 50      | 27  | 788 | 0.64935   | 0.05967 | 0.10929 |

_SVbenchmark_

| Tool  | File                     | TP_base | FP  | TP_comp | FN  | Recall | Precision | F1                |
| ----- | ------------------------ | ------- | --- | ------- | --- | ------ | --------- | ----------------- |
| test1 | test1.HG002.manta.report | 195     | 142 | 195     | 643 | 0.2327 | 0.5684    | 0.330207687187869 |
| test2 | test2.HG002.lumpy.report | 30      | 29  | 30      | 808 | 0.0358 | 0.5       | 0.066815144766147 |
| test3 | test3.HG002.delly.report | 18      | 60  | 18      | 820 | 0.0215 | 0.2208    | 0.039150460593654 |

_Wittyer_

| Tool  | File                   | StatsType | TP_base | TP_comp | FP       | FN    | Precision             | Recall               | F1                  |
| ----- | ---------------------- | --------- | ------- | ------- | -------- | ----- | --------------------- | -------------------- | ------------------- |
| test1 | test1.HG002.manta.json | Event     | 222     | 192     | 137      | 616   | 0.583586626139817     | 0.26491646778042904  | 0.364410474749288   |
| test1 | test1.HG002.manta.json | Base      | 21082   | 21082   | 43165    | 31669 | 0.32813983532305      | 0.39965119144660705  | 0.36038222875604703 |
| test2 | test2.HG002.lumpy.json | Event     | 31      | 30      | 28       | 807   | 0.5172413793103441    | 0.036992840095465    | 0.069047442274853   |
| test2 | test2.HG002.lumpy.json | Base      | 18845   | 18845   | 41765    | 33906 | 0.310922290051146     | 0.35724441242820004  | 0.332477659865385   |
| test3 | test3.HG002.delly.json | Event     | 22      | 21      | 52       | 816   | 0.28767123287671204   | 0.026252983293556003 | 0.048114976046656   |
| test3 | test3.HG002.delly.json | Base      | 14271   | 14271   | 27535165 | 38480 | 0.0005180142344830001 | 0.270535155731644    | 0.001034048497678   |

The number of TPs found in SVbenchmark is significanly lower than Truvari and Wittyer yet this is not primarly as methodological differences but also because of differences of the parameters defining SV comparions. Therefore, it is highly important to set meaningful parameters before starting to perform benchmarks. Furthermore, precisions in Truvari are higher than two other methods as we used more relaxed parameters for truvari.

## Case 2 : Germline BND (Breakends representation of structural variants) benchmarking

### Config file

- [test config](../conf/tests/germline_bnd.config)

### Analysis

- 4 different publicly available structural variant calls for HG002 are being benchmarked against Genome in a Bottle HG002 SV analysis. Only chromosome 21 will be used for the analysis.
- Size filtering by >30 bases applied to test files to ensure equal SV size in benchmarking.
- svdecompose is used to convert SVTYPEs to "BND"
- RTGtools bndeval is used as benchmarking method.

### Results

| Tool  | File                          | Threshold | TP_base | TP_call | FP   | FN  | Precision | Recall | F1     |
| ----- | ----------------------------- | --------- | ------- | ------- | ---- | --- | --------- | ------ | ------ |
| test1 | test1.HG002.manta.summary.txt | None      | 163     | 165     | 1302 | 295 | 0.1125    | 0.3559 | 0.1709 |
| test2 | test2.HG002.lumpy.summary.txt | None      | 48      | 48      | 279  | 410 | 0.1468    | 0.1048 | 0.1223 |
| test3 | test3.HG002.delly.summary.txt | None      | 52      | 54      | 338  | 406 | 0.1378    | 0.1135 | 0.1245 |
| test4 | test4.HG002.svaba.summary.txt | None      | 16      | 16      | 49   | 442 | 0.2462    | 0.0349 | 0.0612 |

Note that RTGtools bndeval is only configured to benchmark SVTYPE=BND. Eventough we decomposed SV tpes to BND, SVs with other types can remain and they might be matched with the tool. That is why both precison and recall values are signififcanly lower than SVBenchmark and Truvari. bndeval should only be used for BND SV analysis with corresponding truths.

## Case 3 : Lifting over one test file for Structural Variant benchmarking

### Config file

- [test config](../conf/tests/liftover_test.config)

### Analysis

- Besides to the 4 test cases, we are adding another test case here whom is generated through GRCh38 reference instead of GRCh37. This test case will be lifted over to GRCh37 by the pipeline.
- Preprocessing includes normalization enabling left alignment of the variants, splitting multi allelic variants, and deduplication.
- Size filtering by >30 bases applied to test files to ensure equal SV size in benchmarking.
- Truvari is used as benchmarking method. Tools spesific parameters are given in the corresponding samplesheet.

### Results

_Truvari_

| Tool  | File                            | TP_base | TP_comp | FP  | FN  | Precision | Recall | F1     |
| ----- | ------------------------------- | ------- | ------- | --- | --- | --------- | ------ | ------ |
| test1 | test1.HG002.manta.summary.json  | 83      | 83      | 113 | 94  | 0.4235    | 0.4689 | 0.4450 |
| test2 | test2.HG002.lumpy.summary.json  | 24      | 24      | 65  | 153 | 0.2697    | 0.1356 | 0.1805 |
| test3 | test3.HG002.delly.summary.json  | 11      | 11      | 52  | 166 | 0.1746    | 0.0621 | 0.0917 |
| test4 | test4.HG002.svaba.summary.json  | 8       | 8       | 22  | 169 | 0.2667    | 0.0452 | 0.0773 |
| test5 | test5.HG002.dragen.summary.json | 8       | 8       | 19  | 169 | 0.2963    | 0.0452 | 0.0784 |

Now we can also take test5 into account. test5 originally includes small variants, since we filtered them by lenght they are not involved to the analysis.

## Case 4 : Germline Small Variant benchmarking

### Config file

- [test config](../conf/tests/germline_small.config)

### Analysis

- Now we are using two HG002 variant calls from a germline analysis for small variant benchmarking. Chromosome 21 is extracted both from the test and truth cases.
- We normalize, deduplicate and prepy to preprocess the variants
- Both RTGtools vcfeval and hap.py methods are used for benchmarking.

### Results

_Hap.py_

| Tool  | File                             | Type  | Filter | TP_base | TP    | FN  | TP_call | FP  | UNK | FP_gt | FP_al | Recall   | Precision | Frac_NA | F1       | TRUTH_TiTv_ratio | QUERY_TiTv_ratio   | TRUTH_het_hom_ratio | QUERY_het_hom_ratio |
| ----- | -------------------------------- | ----- | ------ | ------- | ----- | --- | ------- | --- | --- | ----- | ----- | -------- | --------- | ------- | -------- | ---------------- | ------------------ | ------------------- | ------------------- |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | ALL    | 7675    | 7581  | 94  | 7889    | 88  | 0   | 10    | 42    | 0.987752 | 0.988845  | 0.0     | 0.988299 | 0.0              | 0.0                | 1.087144089732528   | 1.2469600463231036  |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | PASS   | 7675    | 7567  | 108 | 7841    | 53  | 0   | 2     | 31    | 0.985928 | 0.993241  | 0.0     | 0.989571 | 0.0              | 0.0                | 1.087144089732528   | 1.2372389791183294  |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | ALL    | 45057   | 44701 | 356 | 45380   | 622 | 0   | 30    | 33    | 0.992099 | 0.986294  | 0.0     | 0.989188 | 2.14825682945574 | 2.1273428886438808 | 1.434554690877648   | 1.4639947865754317  |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | PASS   | 45057   | 44602 | 455 | 44767   | 110 | 0   | 10    | 5     | 0.989902 | 0.997543  | 0.0     | 0.993708 | 2.14825682945574 | 2.140281966753174  | 1.434554690877648   | 1.434936350777935   |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | ALL    | 7675    | 7212  | 463 | 7741    | 330 | 0   | 252   | 61    | 0.939674 | 0.95737   | 0.0     | 0.94844  | 0.0              | 0.0                | 1.087144089732528   | 1.3940917661847894  |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | PASS   | 7675    | 7212  | 463 | 7741    | 330 | 0   | 252   | 61    | 0.939674 | 0.95737   | 0.0     | 0.94844  | 0.0              | 0.0                | 1.087144089732528   | 1.3940917661847894  |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | ALL    | 45057   | 44703 | 354 | 45708   | 949 | 0   | 32    | 75    | 0.992143 | 0.979238  | 0.0     | 0.985648 | 2.14825682945574 | 2.1275314723590584 | 1.434554690877648   | 1.4739605889995668  |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | PASS   | 45057   | 44703 | 354 | 45708   | 949 | 0   | 32    | 75    | 0.992143 | 0.979238  | 0.0     | 0.985648 | 2.14825682945574 | 2.1275314723590584 | 1.434554690877648   | 1.4739605889995668  |

_RTGtools_

| Tool  | File                             | Threshold | TP_base | TP_call | FP   | FN  | Precision | Recall | F1     |
| ----- | -------------------------------- | --------- | ------- | ------- | ---- | --- | --------- | ------ | ------ |
| test6 | test6.HG002.strelka.summary.txt  | 30.000    | 52159   | 52183   | 405  | 569 | 0.9923    | 0.9892 | 0.9908 |
| test6 | test6.HG002.strelka.summary.txt  | None      | 52277   | 52301   | 708  | 451 | 0.9866    | 0.9914 | 0.989  |
| test7 | test7.HG002.bcftools.summary.txt | None      | 51916   | 51943   | 1270 | 812 | 0.9761    | 0.9846 | 0.9803 |

_Aardvark_

| Tool | File | Caller | Comparison | Region | Filter | Type | TP_base | TP_comp | FP | FN | Precision | Recall | F1 | FN_gt | FP_gt |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| test6 | test6.HG002.strelka.summary.tsv | strelka | GT | ALL | ALL | ALL | 53153 | 52739 | 662 | 449 | 0.9876 | 0.9916 | 0.9896 | 43 | 15 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | GT | ALL | ALL | Snv | 45084 | 44796 | 592 | 356 | 0.9870 | 0.9921 | 0.9895 | 30 | 2 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | GT | ALL | ALL | Insertion | 4011 | 3934 | 25 | 50 | 0.9937 | 0.9875 | 0.9906 | 5 | 6 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | GT | ALL | ALL | Deletion | 4058 | 4003 | 45 | 43 | 0.9889 | 0.9894 | 0.9891 | 8 | 7 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | GT | ALL | ALL | Indel | 0 | 6 | 0 | 0 | 1.0000 | 0.0000 | 0.0000 | 0 | 0 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | GT | ALL | ALL | JointIndel | 8069 | 7943 | 70 | 93 | 0.9913 | 0.9885 | 0.9899 | 13 | 13 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | BASEPAIR | ALL | ALL | ALL | 188294 | 186596 | 1618 | 1698 | 0.9914 | 0.9910 | 0.9912 | 0 | 0 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | BASEPAIR | ALL | ALL | Snv | 127176 | 126367 | 1243 | 900 | 0.9903 | 0.9929 | 0.9916 | 0 | 0 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | BASEPAIR | ALL | ALL | Insertion | 28894 | 28136 | 260 | 705 | 0.9908 | 0.9756 | 0.9832 | 0 | 0 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | BASEPAIR | ALL | ALL | Deletion | 32476 | 32042 | 406 | 336 | 0.9875 | 0.9897 | 0.9886 | 0 | 0 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | BASEPAIR | ALL | ALL | Indel | 0 | 62 | 0 | 0 | 1.0000 | 0.0000 | 0.0000 | 0 | 0 |
| test6 | test6.HG002.strelka.summary.tsv | strelka | BASEPAIR | ALL | ALL | JointIndel | 61370 | 60240 | 666 | 1041 | 0.9891 | 0.9830 | 0.9860 | 0 | 0 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT | ALL | ALL | ALL | 53153 | 52625 | 952 | 812 | 0.9822 | 0.9847 | 0.9835 | 288 | 17 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT | ALL | ALL | Snv | 45084 | 44809 | 908 | 352 | 0.9801 | 0.9922 | 0.9861 | 31 | 6 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT | ALL | ALL | Insertion | 4011 | 3886 | 17 | 221 | 0.9956 | 0.9449 | 0.9696 | 129 | 4 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT | ALL | ALL | Deletion | 4058 | 3930 | 27 | 239 | 0.9932 | 0.9411 | 0.9664 | 128 | 7 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT | ALL | ALL | JointIndel | 8069 | 7816 | 44 | 460 | 0.9944 | 0.9430 | 0.9680 | 257 | 11 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR | ALL | ALL | ALL | 188314 | 183388 | 2100 | 4926 | 0.9887 | 0.9738 | 0.9812 | 0 | 0 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR | ALL | ALL | Snv | 127176 | 126558 | 1804 | 848 | 0.9859 | 0.9933 | 0.9896 | 0 | 0 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR | ALL | ALL | Insertion | 28894 | 26509 | 155 | 2342 | 0.9942 | 0.9189 | 0.9551 | 0 | 0 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR | ALL | ALL | Deletion | 32476 | 30477 | 243 | 1942 | 0.9921 | 0.9402 | 0.9654 | 0 | 0 |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR | ALL | ALL | JointIndel | 61370 | 56986 | 398 | 4284 | 0.9931 | 0.9302 | 0.9606 | 0 | 0 |

While hap.py reports the metrics for SNPs and INDELs separately, RTGtools merges all variant types into a single aggregate score. Aardvark takes the most granular approach, breaking down performance not only by highly specific variant types (SNV, Insertion, Deletion, JointIndel) but also by the matching criteria (Genotype GT vs. BASEPAIR).

As you can see, even for small variant benchmarking, the exact results vary depending on the tool. For example, looking at the unfiltered test6 (Strelka) overall error counts:

- Aardvark (GT, ALL): 449 FNs / 662 FPs

- Hap.py (SNP + INDEL, ALL): 450 FNs / 710 FPs

- RTGtools (No threshold): 451 FNs / 708 FPs

These slight but distinct differences—especially in the False Positive (FP) counts—emerge from how each tool's internal engine handles variant normalization, complex indel alignment, and allele string matching. This highlights exactly why standardizing the benchmarking tool within a pipeline is just as critical as standardizing the variant caller itself.

## Case 5 : Lifting over truth vcf for Germline Small Variant benchmarking

### Config file

- [test config](../conf/tests/liftover_truth.config)

### Analysis

- We are using the same test samples used in _Case 4_ but in this case, we will try to becnhamark them against a truth file which is generated with GRCh37 using lifting over strategy.
- Both RTGtools vcfeval and hap.py methods are used for benchmarking.

### Results

_Hap.py_

| Tool  | File                             | Type  | Filter | TP_base | TP    | FN  | TP_call | FP  | UNK | FP_gt | FP_al | Recall   | Precision | Frac_NA | F1       | TRUTH_TiTv_ratio | QUERY_TiTv_ratio | TRUTH_het_hom_ratio | QUERY_het_hom_ratio |
| ----- | -------------------------------- | ----- | ------ | ------- | ----- | --- | ------- | --- | --- | ----- | ----- | -------- | --------- | ------- | -------- | ---------------- | ---------------- | ------------------- | ------------------- |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | ALL    | 7737    | 7560  | 177 | 7864    | 78  | 0   | 10    | 45    | 0.977123 | 0.990081  | 0.0     | 0.983559 | 0.0              | 0.0              | 1.053963            | 1.249636            |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | PASS   | 7737    | 7545  | 192 | 7817    | 43  | 0   | 1     | 33    | 0.975184 | 0.994499  | 0.0     | 0.984747 | 0.0              | 0.0              | 1.053963            | 1.239219            |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | ALL    | 44530   | 44467 | 63  | 44898   | 375 | 0   | 26    | 36    | 0.998585 | 0.991648  | 0.0     | 0.995104 | 2.148685         | 2.140509         | 1.434783            | 1.446479            |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | PASS   | 44530   | 44389 | 141 | 44564   | 122 | 0   | 9     | 7     | 0.996834 | 0.997262  | 0.0     | 0.997048 | 2.148685         | 2.148195         | 1.434783            | 1.431323            |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | ALL    | 7737    | 7196  | 541 | 7737    | 335 | 0   | 252   | 63    | 0.930076 | 0.956702  | 0.0     | 0.943201 | 0.0              | 0.0              | 1.053963            | 1.398046            |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | PASS   | 7737    | 7196  | 541 | 7737    | 335 | 0   | 252   | 63    | 0.930076 | 0.956702  | 0.0     | 0.943201 | 0.0              | 0.0              | 1.053963            | 1.398046            |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | ALL    | 44530   | 44451 | 79  | 45236   | 731 | 0   | 29    | 72    | 0.998226 | 0.983840  | 0.0     | 0.990981 | 2.148685         | 2.139269         | 1.434783            | 1.460743            |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | PASS   | 44530   | 44451 | 79  | 45236   | 731 | 0   | 29    | 72    | 0.998226 | 0.983840  | 0.0     | 0.990981 | 2.148685         | 2.139269         | 1.434783            | 1.460743            |

_RTGtools_

| Tool  | File                             | Threshold | TP_base | TP_call | FP   | FN  | Precision | Recall | F1     |
| ----- | -------------------------------- | --------- | ------- | ------- | ---- | --- | --------- | ------ | ------ |
| test6 | test6.HG002.strelka.summary.txt  | 24.000    | 51961   | 51986   | 301  | 303 | 0.9942    | 0.9942 | 0.9942 |
| test6 | test6.HG002.strelka.summary.txt  | None      | 52025   | 52050   | 449  | 239 | 0.9914    | 0.9954 | 0.9934 |
| test7 | test7.HG002.bcftools.summary.txt | None      | 51648   | 51673   | 1059 | 616 | 0.9799    | 0.9882 | 0.984  |

As you can see, lifting over files from GRCh37 to GRCh38 worked quite well, there is not too much difference.

## Case 6 : GH4GH Germline Small Variant benchmarking

- [test config](../conf/tests/test_ga4gh.config)

### Analysis

- Again, using the same test samples used in _Case 4_, we will apply GA4GH best practice benchmark method, running hap.py tool as benchmark engine.

### Results

_Hap.py_

| Tool  | File                             | Type  | Filter | TP_base | TP    | FN  | TP_call | FP  | UNK   | FP_gt | FP_al | Recall   | Precision | Frac_NA  | F1       | TRUTH_TiTv_ratio  | QUERY_TiTv_ratio   | TRUTH_het_hom_ratio | QUERY_het_hom_ratio |
| ----- | -------------------------------- | ----- | ------ | ------- | ----- | --- | ------- | --- | ----- | ----- | ----- | -------- | --------- | -------- | -------- | ----------------- | ------------------ | ------------------- | ------------------- |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | ALL    | 8074    | 7970  | 104 | 14216   | 79  | 6577  | 25    | 12    | 0.987119 | 0.989658  | 0.462648 | 0.988387 | 0.0               | 0.0                | 1.3221167673281564  | 1.589951690821256   |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | PASS   | 8074    | 7949  | 125 | 13225   | 47  | 5633  | 15    | 7     | 0.984518 | 0.993809  | 0.425936 | 0.989142 | 0.0               | 0.0                | 1.3221167673281564  | 1.4346123727486295  |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | ALL    | 45084   | 44727 | 357 | 68565   | 619 | 23185 | 31    | 46    | 0.992081 | 0.98636   | 0.338146 | 0.989212 | 2.146346569893224 | 1.857571339304312  | 1.4364461738002594  | 1.8474484707446808  |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | PASS   | 45084   | 44628 | 456 | 52951   | 108 | 8184  | 11    | 3     | 0.989886 | 0.997588  | 0.154558 | 0.993722 | 2.146346569893224 | 2.042154835745463  | 1.4364461738002594  | 1.351699622306154   |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | ALL    | 8074    | 7602  | 472 | 12415   | 295 | 4906  | 263   | 22    | 0.941541 | 0.960714  | 0.395167 | 0.951031 | 0.0               | 0.0                | 1.3221167673281564  | 1.6045832397214108  |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | PASS   | 8074    | 7602  | 472 | 12415   | 295 | 4906  | 263   | 22    | 0.941541 | 0.960714  | 0.395167 | 0.951031 | 0.0               | 0.0                | 1.3221167673281564  | 1.6045832397214108  |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | ALL    | 45084   | 44729 | 355 | 73238   | 943 | 27530 | 43    | 153   | 0.992126 | 0.979369  | 0.375898 | 0.985706 | 2.146346569893224 | 1.8436966309579257 | 1.4364461738002594  | 2.03402967752632    |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | PASS   | 45084   | 44729 | 355 | 73238   | 943 | 27530 | 43    | 153   | 0.992126 | 0.979369  | 0.375898 | 0.985706 | 2.146346569893224 | 1.8436966309579257 | 1.4364461738002594  | 2.03402967752632    |

Don't forget that the only difference between the cases in 4,5 and 6 are the methods applied the same sample sets and still the results are different.

## Case 7 : Germline Small Variant benchmarking with stratifications - TODO

- [test config](../conf/tests/test_happy.config)

### Analysis

- Using the same test samples used in _Case 4_, we will apply stratifications to hap.py.

## Case 8 : Somatic SNV Variant benchmarking

- [test config](../conf/tests/somatic_snv.config)

### Analysis

- Benchmarking somatic variants is possible using this pipeline. In order to demonstare SNV type of benchmarking, 3 example somatic variant calls will be benchmarked against SEQC2 truth SNV calls.
- Preprocessing includes normalization enabling left alignment of the variants, deduplication and prepy (only for som.py).
- Using som.py (a version of hap.py tuned specially for somatic benchmarking) and rtgtools vcfeval with "--squash-ploidy" parameter on.

### Results

_Som.py_

| Tool   | File                            | Type | TP_base | TP   | FN    | TP_call | FP   | UNK | Recall             | Precision          | recall_lower          | recall_upper       | recall2            | precision_lower    | precision_upper    | na  | ambiguous | fp.region.size | F1                    |
| ------ | ------------------------------- | ---- | ------- | ---- | ----- | ------- | ---- | --- | ------------------ | ------------------ | --------------------- | ------------------ | ------------------ | ------------------ | ------------------ | --- | --------- | -------------- | --------------------- |
| test8  | test8.SEQC2.freebayes.stats.csv | SNVs | 39447   | 7    | 39440 | 1808    | 1801 | 0   | 0.0001774532917585 | 0.0038716814159292 | 7.937745863543826e-05 | 0.0003483898507791 | 0.0001774532917585 | 0.0017334020987412 | 0.0075866570604652 | 0.0 | 0.0       | 2875001522.0   | 0.0003393528057204051 |
| test9  | test9.SEQC2.mutect2.stats.csv   | SNVs | 39447   | 6    | 39441 | 6       | 0    | 0   | 0.0001521028215073 | 1.0                | 6.348950120714873e-05 | 0.0003135023155461 | 0.0001521028215073 | 0.5407418735600996 | 1.0                | 0.0 | 0.0       | 2875001522.0   | 0.0003041593795147879 |
| test10 | test10.SEQC2.strelka.stats.csv  | SNVs | 39447   | 1419 | 38028 | 2604    | 1185 | 0   | 0.0359723172864856 | 0.5449308755760369 | 0.0341685494392541    | 0.0378441736192306 | 0.0359723172864856 | 0.525763248380803  | 0.5639986848401537 | 0.0 | 0.0       | 2875001522.0   | 0.06748947706356555   |

_RTGtools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP   | FN    | Precision | Recall | F1     |
| ------ | --------------------------------- | --------- | ------- | ------- | ---- | ----- | --------- | ------ | ------ |
| test8  | test8.SEQC2.freebayes.summary.txt | 61.000    | 7       | 7       | 1436 | 39440 | 0.0049    | 0.0002 | 0.0003 |
| test8  | test8.SEQC2.freebayes.summary.txt | None      | 7       | 7       | 1801 | 39440 | 0.0039    | 0.0002 | 0.0003 |
| test9  | test9.SEQC2.mutect2.summary.txt   | None      | 8       | 8       | 11   | 39439 | 0.4211    | 0.0002 | 0.0004 |
| test10 | test10.SEQC2.strelka.summary.txt  | None      | 0       | 0       | 0    | 39447 |           | 0.0    |        |

The number of TPs found quite similar with som.py and RTGtools, yet test10 was not able to run properly with rtgtools as GT field is not reported in strelka.

## Case 9 : Somatic INDEL Variant benchmarking

- [test config](../conf/tests/somatic_indel.config)

### Analysis

- In order to demonstare INDEL type of benchmarking for somatic variants, 3 example somatic variant calls will be benchmarked against SEQC2 truth INDEL calls.
- This time we are not applying preprocessing to the test variants.
- Using som.py (a version of hap.py tuned specially for somatic benchmarking) and rtgtools vcfeval with "--squash-ploidy" parameter on.

### Results

_som.py_

| Tool   | File                            | Type   | TP_base | TP  | FN   | TP_call | FP  | UNK | Recall             | Precision          | recall_lower       | recall_upper       | recall2            | precision_lower    | precision_upper    | na  | ambiguous | fp.region.size | F1                    |
| ------ | ------------------------------- | ------ | ------- | --- | ---- | ------- | --- | --- | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | --- | --------- | -------------- | --------------------- |
| test8  | test8.SEQC2.freebayes.stats.csv | indels | 1625    | 1   | 1624 | 250     | 249 | 0   | 0.0006153846153846 | 0.004              | 0.0                | 0.0028727403989354 | 0.0006153846153846 | 0.0                | 0.0185415188410741 | 0.0 | 0.0       | 2875001522.0   | 0.0010666666666666435 |
| test9  | test9.SEQC2.mutect2.stats.csv   | indels | 1625    | 1   | 1624 | 1       | 0   | 0   | 0.0006153846153846 | 1.0                | 0.0                | 0.0028727403989354 | 0.0006153846153846 | 0.025              | 1.0                | 0.0 | 0.0       | 2875001522.0   | 0.0012300123001229705 |
| test10 | test10.SEQC2.strelka.stats.csv  | indels | 1625    | 55  | 1570 | 159     | 104 | 0   | 0.0338461538461538 | 0.3459119496855345 | 0.0258659328642734 | 0.0434839942668625 | 0.0338461538461538 | 0.2752802625605345 | 0.4220960778452732 | 0.0 | 0.0       | 2875001522.0   | 0.06165919282511203   |

_RTGtools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP  | FN   | Precision | Recall | F1     |
| ------ | --------------------------------- | --------- | ------- | ------- | --- | ---- | --------- | ------ | ------ |
| test8  | test8.SEQC2.freebayes.summary.txt | 99.000    | 1       | 1       | 148 | 1623 | 0.0067    | 0.0006 | 0.0011 |
| test8  | test8.SEQC2.freebayes.summary.txt | None      | 1       | 1       | 249 | 1623 | 0.004     | 0.0006 | 0.0011 |
| test9  | test9.SEQC2.mutect2.summary.txt   | None      | 1       | 1       | 9   | 1623 | 0.1       | 0.0006 | 0.0012 |
| test10 | test10.SEQC2.strelka.summary.txt  | None      | 0       | 0       | 0   | 1624 |           | 0.0    |        |

As the truth set contains very few variants, the number of matcheds quite small and rtgtools couldnt able to find any.

## Case 10 : Somatic SV Variant benchmarking

- [test config](../conf/tests/somatic_sv.config)

### Analysis

- In this case we are using nf-core/sarek variant calls from manta and tiddit. AWS stored files are being pulled automatically with the pipeline.
- We are filtering out extra contigs, and taking only the variants longer than 30bp.
- We use truvari as benchmarking method and all the method related parameters are in samplesheet.
- Please note that SEQC2 SV truth file used in this case is not validated, it is generated only for testing purposes using the dataset reported in (Talsania et all.)[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02816-6]

### Results

| Tool   | File                             | TP_base | TP_comp | FP    | FN  | Precision             | Recall               | F1                    |
| ------ | -------------------------------- | ------- | ------- | ----- | --- | --------------------- | -------------------- | --------------------- |
| test11 | test11.SEQC2.manta.summary.json  | 2       | 2       | 962   | 587 | 0.002074688796680498  | 0.003395585738539898 | 0.0025756600128782996 |
| test12 | test12.SEQC2.tiddit.summary.json | 38      | 38      | 12319 | 551 | 0.0030751800598850854 | 0.06451612903225806  | 0.005870539162675731  |

## Case 11 : Somatic CNV Variant benchmarking

- [test config](../conf/tests/somatic_cnv.config)

### Analysis

- We are again using nf-core/sarek variant calls from cnvkit, ascat, controlfreec and manta.
- We are filtering out extra contigs and splitting multi-allelic sites.
- We use truvari, wittyer and _intersect_ methods. _intersect_ (using bedtools) is not a benchmarking method but rather to intersect BED regions given for test and truth files. truvari and wittyer will be only applied if input test is in VCF format reming that CNV results are tend to report in BED similar formats.
- The truth file used here is SEQC2 SV truth file.

### Results

_RTG Tools_

| Tool   | File                            | Caller | Threshold | TP_base | TP_comp | FP  | FN  | Precision | Recall | F1     |
| ------ | ------------------------------- | ------ | --------- | ------- | ------- | --- | --- | --------- | ------ | ------ |
| test11 | test11.SEQC2.manta.summary.txt  | manta  | None      | 28      | 28      | 156 | 334 | 0.1522    | 0.0773 | 0.1026 |
| test14 | test14.SEQC2.cnvkit.summary.txt | cnvkit | None      | 151     | 151     | 916 | 211 | 0.1415    | 0.4171 | 0.2113 |

_Truvari_

| Tool   | File                             | Caller | TP_base | TP_comp | FP  | FN   | Precision            | Recall                | F1                   |
| ------ | -------------------------------- | ------ | ------- | ------- | --- | ---- | -------------------- | --------------------- | -------------------- |
| test11 | test11.SEQC2.manta.summary.json  | manta  | 5       | 5       | 695 | 1114 | 0.007142857142857143 | 0.004468275245755138  | 0.005497526113249038 |
| test14 | test14.SEQC2.cnvkit.summary.json | cnvkit | 8       | 8       | 96  | 1111 | 0.07692307692307693  | 0.0071492403932082215 | 0.013082583810302535 |

_Witty.er_

| Tool   | File                     | Caller | StatsType | TP_base | TP_comp | FP        | FN        | Precision | Recall | F1  |
| ------ | ------------------------ | ------ | --------- | ------- | ------- | --------- | --------- | --------- | ------ | --- |
| test14 | test14.SEQC2.cnvkit.json | cnvkit | Event     | 0       | 0       | 104       | 1032      | 0.0       | 0.0    |     |
| test14 | test14.SEQC2.cnvkit.json | cnvkit | Base      | 0       | 0       | 525553404 | 115803512 | 0.0       | 0.0    |     |
| test11 | test11.SEQC2.manta.json  | manta  | Event     | 0       | 0       | 304       | 1032      | 0.0       | 0.0    |     |
| test11 | test11.SEQC2.manta.json  | manta  | Base      | 0       | 0       | 64244018  | 149089348 | 0.0       | 0.0    |     |

_Intersect_

| Tool             | Caller             | Analysis           | File                                    | TP_base | TP_comp | FN   | FP   | Precision          | Recall             | F1                 |
| ---------------- | ------------------ | ------------------ | --------------------------------------- | ------- | ------- | ---- | ---- | ------------------ | ------------------ | ------------------ |
| test13           | controlfreec_stats | Bedtools_Intersect | test13.SEQC2.controlfreec_stats.csv     | 95      | 95      | 1407 | 344  | 0.2164009111617312 | 0.0632490013315579 | 0.0978876867594023 |
| test14           | cnvkit_stats       | Bedtools_Intersect | test14.SEQC2.cnvkit_stats.csv           | 123     | 123     | 1379 | 259  | 0.3219895287958115 | 0.0818908122503328 | 0.1305732484076433 |
| test15           | ascat_stats        | Bedtools_Intersect | test15.SEQC2.ascat_stats.csv            | 64      | 64      | 1438 | 119  | 0.3497267759562841 | 0.0426098535286285 | 0.0759643916913946 |
| test11_converted | manta              | Bedtools_Intersect | test11.SEQC2.manta.converted_stats.csv  | 1       | 1       | 1501 | 1063 | 0.0009398496240601 | 0.0006657789613848 | 0.0007794232268121 |
| test14_converted | cnvkit             | Bedtools_Intersect | test14.SEQC2.cnvkit.converted_stats.csv | 90      | 90      | 1412 | 207  | 0.303030303030303  | 0.0599201065246338 | 0.1000555864369094 |

Only CNVkit was reporting the variants in VCF format that is why truvari, wittyer and rgtools was running only that test sample. Wittyer was not proprely running as variants are not sequence resolved.
We also have two type of results for CNVkit for intersection analysis. \*\_converted one show the result from the VCF input being converted to BED for ntersection analysis while the other result is using input BED file directly.

## Case 12 : Concordance analysis

In order to do concordance comparison, we dont need to provide truth VCF! as we are comparing test VCFs with each other. This type of analysis can be usefull for pairwise comparisons. Concordance analysis is only avilable for small variants for now.

### Config file

- [test config](../conf/tests/concordance.config)

### Analysis

- Now we are using two HG002 variant calls from germline variant calling analysis to perform concordance analysis.
- Chromosome 21 is extracted both from the test and truth cases.
- We normalize, split multiallelics and deduplicate variants to preprocess the variants

### Results

_GATK4 Concordance_

| Tool        | File                    | Type  | TP    | FP    | FN   | Precision | Recall | F1                 |
| ----------- | ----------------------- | ----- | ----- | ----- | ---- | --------- | ------ | ------------------ |
| test7-test6 | test7-test6.summary.tsv | SNP   | 52698 | 20579 | 279  | 0.995     | 0.719  | 0.8347782963827304 |
| test7-test6 | test7-test6.summary.tsv | INDEL | 12129 | 1108  | 1885 | 0.865     | 0.916  | 0.8897697922515441 |

Concordance statistics are based on tool comparisons, the pairs are generated randomly for each run therefore the order of the comparison may change. Here, test6 is considered as BASE whole test7 is being considered as COMP set. Therefore, be aware that FP, FN and Precisison and Recall may switch.

> [!NOTE]
> As a final note, it is highly encouraged to run all those test and investigate comparisons, plots and tables created. This pipeline serves for multiple tools enabling efficient benchmark analysis.

## Case 13 : Ensembling test files to generate a proxy truth for small variants

Truth/Golden set of samples to use as a baseline for benchmarking analysis is a challenge especially for somatic studies. One of the common approaches reseacrhes uses is to generate a proxy set of truth using input (test) VCF files, first ensembling the variants and then filtering them down with a majority vote system.

### Config file

- [test config](../conf/tests/somatic_snv_ensemble.config)

### Analysis

- Benchmarking somatic variants is possible using this pipeline. In order to demonstare SNV type of benchmarking, 3 example somatic variant calls will be used.
- Truth generation is done accordingly: ensemble_truth is used as 2, so the variant is kept only if it exist in 2 or more callers (test vcfs)
- Preprocessing includes normalization enabling left alignment of the variants.
- Using som.py.

### Results

_Som.py_

| Tool   | File                            | Type | TP_base | TP  | FN  | TP_call | FP   | UNK | Recall             | Precision          | recall_lower       | recall_upper | recall2            | precision_lower    | precision_upper    | na  | ambiguous | fp.region.size | F1                   |
| ------ | ------------------------------- | ---- | ------- | --- | --- | ------- | ---- | --- | ------------------ | ------------------ | ------------------ | ------------ | ------------------ | ------------------ | ------------------ | --- | --------- | -------------- | -------------------- |
| test8  | test8.truth.freebayes.stats.csv | SNVs | 7       | 7   | 0   | 1898    | 1891 | 0   | 1.0                | 0.0036880927291886 | 0.5903836027749967 | 1.0          | 1.0                | 0.0016511339351389 | 0.0072275976268013 | 0.0 | 0.0       | 46709983.0     | 0.007349081364829357 |
| test9  | test9.truth.mutect2.stats.csv   | SNVs | 7       | 6   | 1   | 8       | 2    | 0   | 0.8571428571428571 | 0.75               | 0.4992025383809224 | 1.0          | 0.8571428571428571 | 0.4083764894266625 | 0.9440331102753412 | 0.0 | 0.0       | 46709983.0     | 0.7999999999999999   |
| test10 | test10.truth.strelka.stats.csv  | SNVs | 7       | 7   | 0   | 5251    | 5244 | 0   | 1.0                | 0.001333079413445  | 0.5903836027749967 | 1.0          | 1.0                | 0.0005964718908893 | 0.0026156384477055 | 0.0 | 0.0       | 46709983.0     | 0.002662609357169911 |

## Case 14 : Ensembling test files to generate a proxy truth for structual variants

This is another example to generate a proxy truth file using test files for structural variants.

### Config file

- [test config](../conf/tests/somatic_sv_ensemble.config)

### Analysis

- In this case we are using nf-core/sarek variant calls from manta and tiddit. AWS stored files are being pulled automatically with the pipeline.
- We are filtering out extra contigs, and taking only the variants longer than 30bp.
- We use truvari as benchmarking method and all the method related parameters are in samplesheet.
- Truth generation is done accordingly: ensemble_truth is used as 2, so the variant is kept only if it exist in 2 or more callers (test vcfs)

### Results

| Tool   | File                             | TP_base | TP_comp | FP   | FN   | Precision          | Recall              | F1                  |
| ------ | -------------------------------- | ------- | ------- | ---- | ---- | ------------------ | ------------------- | ------------------- |
| test11 | test11.truth.manta.summary.json  | 333     | 333     | 631  | 6447 | 0.3454356846473029 | 0.04911504424778761 | 0.08600206611570248 |
| test12 | test12.truth.tiddit.summary.json | 3344    | 3344    | 9013 | 3436 | 0.2706158452698875 | 0.49321533923303834 | 0.349480064795945   |
