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
- We normalize and deduplicate to preprocess the variants
- RTGtools vcfeval, hap.py and aardvark methods are used for benchmarking.

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

| Tool  | File                             | Caller   | Comparison | Region | Filter | Type       | TP_base | TP_comp | FP   | FN   | Precision | Recall | F1     | FN_gt | FP_gt |
| ----- | -------------------------------- | -------- | ---------- | ------ | ------ | ---------- | ------- | ------- | ---- | ---- | --------- | ------ | ------ | ----- | ----- |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | GT         | ALL    | ALL    | ALL        | 53153   | 52739   | 662  | 449  | 0.9876    | 0.9916 | 0.9896 | 43    | 15    |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | GT         | ALL    | ALL    | Snv        | 45084   | 44796   | 592  | 356  | 0.9870    | 0.9921 | 0.9895 | 30    | 2     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | GT         | ALL    | ALL    | Insertion  | 4011    | 3934    | 25   | 50   | 0.9937    | 0.9875 | 0.9906 | 5     | 6     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | GT         | ALL    | ALL    | Deletion   | 4058    | 4003    | 45   | 43   | 0.9889    | 0.9894 | 0.9891 | 8     | 7     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | GT         | ALL    | ALL    | Indel      | 0       | 6       | 0    | 0    | 1.0000    | 0.0000 | 0.0000 | 0     | 0     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | GT         | ALL    | ALL    | JointIndel | 8069    | 7943    | 70   | 93   | 0.9913    | 0.9885 | 0.9899 | 13    | 13    |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | BASEPAIR   | ALL    | ALL    | ALL        | 188294  | 186596  | 1618 | 1698 | 0.9914    | 0.9910 | 0.9912 | 0     | 0     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | BASEPAIR   | ALL    | ALL    | Snv        | 127176  | 126367  | 1243 | 900  | 0.9903    | 0.9929 | 0.9916 | 0     | 0     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | BASEPAIR   | ALL    | ALL    | Insertion  | 28894   | 28136   | 260  | 705  | 0.9908    | 0.9756 | 0.9832 | 0     | 0     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | BASEPAIR   | ALL    | ALL    | Deletion   | 32476   | 32042   | 406  | 336  | 0.9875    | 0.9897 | 0.9886 | 0     | 0     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | BASEPAIR   | ALL    | ALL    | Indel      | 0       | 62      | 0    | 0    | 1.0000    | 0.0000 | 0.0000 | 0     | 0     |
| test6 | test6.HG002.strelka.summary.tsv  | strelka  | BASEPAIR   | ALL    | ALL    | JointIndel | 61370   | 60240   | 666  | 1041 | 0.9891    | 0.9830 | 0.9860 | 0     | 0     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT         | ALL    | ALL    | ALL        | 53153   | 52625   | 952  | 812  | 0.9822    | 0.9847 | 0.9835 | 288   | 17    |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT         | ALL    | ALL    | Snv        | 45084   | 44809   | 908  | 352  | 0.9801    | 0.9922 | 0.9861 | 31    | 6     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT         | ALL    | ALL    | Insertion  | 4011    | 3886    | 17   | 221  | 0.9956    | 0.9449 | 0.9696 | 129   | 4     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT         | ALL    | ALL    | Deletion   | 4058    | 3930    | 27   | 239  | 0.9932    | 0.9411 | 0.9664 | 128   | 7     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | GT         | ALL    | ALL    | JointIndel | 8069    | 7816    | 44   | 460  | 0.9944    | 0.9430 | 0.9680 | 257   | 11    |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR   | ALL    | ALL    | ALL        | 188314  | 183388  | 2100 | 4926 | 0.9887    | 0.9738 | 0.9812 | 0     | 0     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR   | ALL    | ALL    | Snv        | 127176  | 126558  | 1804 | 848  | 0.9859    | 0.9933 | 0.9896 | 0     | 0     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR   | ALL    | ALL    | Insertion  | 28894   | 26509   | 155  | 2342 | 0.9942    | 0.9189 | 0.9551 | 0     | 0     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR   | ALL    | ALL    | Deletion   | 32476   | 30477   | 243  | 1942 | 0.9921    | 0.9402 | 0.9654 | 0     | 0     |
| test7 | test7.HG002.bcftools.summary.tsv | bcftools | BASEPAIR   | ALL    | ALL    | JointIndel | 61370   | 56986   | 398  | 4284 | 0.9931    | 0.9302 | 0.9606 | 0     | 0     |

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

- We are using the same test samples used in _Case 4_ but in this case, we will try to benchmark them against a truth file which is generated with GRCh37 using lifting over strategy.
- Preprocessing includes normalization enabling left alignment of the variants and deduplication.
- Both RTGtools vcfeval and hap.py methods are used for benchmarking.

### Results

_Hap.py_

| Tool  | File                             | Type  | Filter | TP_base | TP    | FN  | TP_call | FP  | UNK  | FP_gt | FP_al | Recall  | Precision | Frac_NA | F1      | TRUTH_TiTv_ratio | QUERY_TiTv_ratio | TRUTH_het_hom_ratio | QUERY_het_hom_ratio |
| ----- | -------------------------------- | ----- | ------ | ------- | ----- | --- | ------- | --- | ---- | ----- | ----- | ------- | --------- | ------- | ------- | ---------------- | ---------------- | ------------------- | ------------------- |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | ALL    | 7431    | 7264  | 167 | 11079   | 57  | 3538 | 10    | 35    | 0.97753 | 0.99244   | 0.31934 | 0.98493 | 0.0              | 0.0              | 1.04604             | 1.39238             |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | PASS   | 7431    | 7252  | 179 | 10765   | 31  | 3259 | 1     | 23    | 0.97591 | 0.99587   | 0.30274 | 0.98579 | 0.0              | 0.0              | 1.04604             | 1.34765             |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | ALL    | 42795   | 42757 | 38  | 47264   | 208 | 4248 | 9     | 23    | 0.99911 | 0.99517   | 0.08988 | 0.99713 | 2.15564          | 2.06047          | 1.43484             | 1.47336             |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | PASS   | 42795   | 42700 | 95  | 45949   | 102 | 3099 | 5     | 6     | 0.99778 | 0.99762   | 0.06744 | 0.9977  | 2.15564          | 2.10436          | 1.43484             | 1.42342             |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | ALL    | 7431    | 6915  | 516 | 10160   | 313 | 2735 | 243   | 56    | 0.93056 | 0.95785   | 0.26919 | 0.94401 | 0.0              | 0.0              | 1.04604             | 1.61431             |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | PASS   | 7431    | 6915  | 516 | 10160   | 313 | 2735 | 243   | 56    | 0.93056 | 0.95785   | 0.26919 | 0.94401 | 0.0              | 0.0              | 1.04604             | 1.61431             |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | ALL    | 42795   | 42740 | 55  | 47517   | 432 | 4295 | 15    | 40    | 0.99872 | 0.99001   | 0.09039 | 0.99434 | 2.15564          | 2.0699           | 1.43484             | 1.45788             |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | PASS   | 42795   | 42740 | 55  | 47517   | 432 | 4295 | 15    | 40    | 0.99872 | 0.99001   | 0.09039 | 0.99434 | 2.15564          | 2.0699           | 1.43484             | 1.45788             |

_RTGtools_

| Tool  | File                             | Threshold | TP_base | TP_call | FP  | FN  | Precision | Recall | F1     |
| ----- | -------------------------------- | --------- | ------- | ------- | --- | --- | --------- | ------ | ------ |
| test6 | test6.HG002.strelka.summary.txt  | 23.000    | 49977   | 49999   | 167 | 248 | 0.9967    | 0.9951 | 0.9959 |
| test6 | test6.HG002.strelka.summary.txt  | None      | 50020   | 50042   | 264 | 205 | 0.9948    | 0.9959 | 0.9953 |
| test7 | test7.HG002.bcftools.summary.txt | None      | 49658   | 49681   | 738 | 567 | 0.9854    | 0.9887 | 0.987  |

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

## Case 7 : Germline Small Variant benchmarking with stratifications

- [test config](../conf/tests/test_happy.config)

### Analysis

- Using the same test samples used in _Case 4_, we will apply stratifications to hap.py.

### Results

_Hap.py_

| Tool  | File                             | Type  | Filter | TP_base | TP   | FN  | TP_call | FP  | UNK  | FP_gt | FP_al | Recall  | Precision | Frac_NA | F1      | TRUTH_TiTv_ratio | QUERY_TiTv_ratio | TRUTH_het_hom_ratio | QUERY_het_hom_ratio |
| ----- | -------------------------------- | ----- | ------ | ------- | ---- | --- | ------- | --- | ---- | ----- | ----- | ------- | --------- | ------- | ------- | ---------------- | ---------------- | ------------------- | ------------------- |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | ALL    | 1536    | 1484 | 52  | 3339    | 30  | 1958 | 20    | 7     | 0.96615 | 0.97828   | 0.5864  | 0.97217 | 0.0              | 0.0              | 1.63918             | 1.55833             |
| test6 | test6.HG002.strelka.summary.csv  | INDEL | PASS   | 1536    | 1472 | 64  | 3072    | 16  | 1712 | 13    | 3     | 0.95833 | 0.98824   | 0.55729 | 0.97306 | 0.0              | 0.0              | 1.63918             | 1.43388             |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | ALL    | 155     | 154  | 1   | 919     | 6   | 753  | 1     | 3     | 0.99355 | 0.96386   | 0.81937 | 0.97848 | 0.64894          | 0.62257          | 1.34848             | 2.31408             |
| test6 | test6.HG002.strelka.summary.csv  | SNP   | PASS   | 155     | 142  | 13  | 417     | 2   | 270  | 1     | 1     | 0.91613 | 0.9864    | 0.64748 | 0.95    | 0.64894          | 0.58333          | 1.34848             | 1.18947             |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | ALL    | 1536    | 1243 | 293 | 2633    | 178 | 1322 | 171   | 7     | 0.80925 | 0.86423   | 0.50209 | 0.83583 | 0.0              | 0.0              | 1.63918             | 2.58383             |
| test7 | test7.HG002.bcftools.summary.csv | INDEL | PASS   | 1536    | 1243 | 293 | 2633    | 178 | 1322 | 171   | 7     | 0.80925 | 0.86423   | 0.50209 | 0.83583 | 0.0              | 0.0              | 1.63918             | 2.58383             |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | ALL    | 155     | 147  | 8   | 641     | 10  | 480  | 5     | 4     | 0.94839 | 0.93789   | 0.74883 | 0.94311 | 0.64894          | 0.69129          | 1.34848             | 1.27758             |
| test7 | test7.HG002.bcftools.summary.csv | SNP   | PASS   | 155     | 147  | 8   | 641     | 10  | 480  | 5     | 4     | 0.94839 | 0.93789   | 0.74883 | 0.94311 | 0.64894          | 0.69129          | 1.34848             | 1.27758             |

## Case 8 : Somatic SNV Variant benchmarking

- [test config](../conf/tests/somatic_snv.config)

### Analysis

- Benchmarking somatic variants is possible using this pipeline. In order to demonstare SNV type of benchmarking, 3 example somatic variant calls will be benchmarked against SEQC2 truth SNV calls.
- Preprocessing includes normalization enabling left alignment of the variants.
- Using som.py (a version of hap.py tuned specially for somatic benchmarking) and rtgtools vcfeval with "--squash-ploidy" parameter on.

### Results

_Som.py_

| Tool   | File                            | Type | TP_base | TP   | FN    | TP_call | FP   | UNK | Recall  | Precision | recall_lower | recall_upper | recall2 | precision_lower | precision_upper | na  | ambiguous | fp.region.size | F1      |
| ------ | ------------------------------- | ---- | ------- | ---- | ----- | ------- | ---- | --- | ------- | --------- | ------------ | ------------ | ------- | --------------- | --------------- | --- | --------- | -------------- | ------- |
| test8  | test8.SEQC2.freebayes.stats.csv | SNVs | 39560   | 7    | 39553 | 1898    | 1891 | 0   | 0.00018 | 0.00369   | 0.00008      | 0.00035      | 0.00018 | 0.00165         | 0.00723         | 0.0 | 0.0       | 2875001522     | 0.00034 |
| test9  | test9.SEQC2.mutect2.stats.csv   | SNVs | 39560   | 6    | 39554 | 8       | 2    | 0   | 0.00015 | 0.75      | 0.00006      | 0.00031      | 0.00015 | 0.40838         | 0.94403         | 0.0 | 0.0       | 2875001522     | 0.0003  |
| test10 | test10.SEQC2.strelka.stats.csv  | SNVs | 39560   | 1421 | 38139 | 5251    | 3830 | 0   | 0.03592 | 0.27062   | 0.03412      | 0.03779      | 0.03592 | 0.25873         | 0.28276         | 0.0 | 0.0       | 2875001522     | 0.06342 |

_RTGtools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP     | FN    | Precision | Recall | F1     |
| ------ | --------------------------------- | --------- | ------- | ------- | ------ | ----- | --------- | ------ | ------ |
| test8  | test8.SEQC2.freebayes.summary.txt | 61.000    | 7       | 7       | 1482   | 39553 | 0.0047    | 0.0002 | 0.0003 |
| test8  | test8.SEQC2.freebayes.summary.txt | None      | 7       | 7       | 1891   | 39553 | 0.0037    | 0.0002 | 0.0003 |
| test9  | test9.SEQC2.mutect2.summary.txt   | None      | 8       | 8       | 15     | 39552 | 0.3478    | 0.0002 | 0.0004 |
| test10 | test10.SEQC2.strelka.summary.txt  | None      | 1552    | 1552    | 111433 | 38008 | 0.0137    | 0.0392 | 0.0203 |

The number of TPs found quite similar with som.py and RTGtools.

## Case 9 : Somatic INDEL Variant benchmarking

- [test config](../conf/tests/somatic_indel.config)

### Analysis

- In order to demonstare INDEL type of benchmarking for somatic variants, 3 example somatic variant calls will be benchmarked against SEQC2 truth INDEL calls.
- This time we are not applying preprocessing to the test variants.
- Using som.py (a version of hap.py tuned specially for somatic benchmarking) and rtgtools vcfeval with "--squash-ploidy" parameter on.

### Results

_som.py_

| Tool   | File                            | Type   | TP_base | TP  | FN   | TP_call | FP  | UNK | Recall  | Precision | recall_lower | recall_upper | recall2 | precision_lower | precision_upper | na  | ambiguous | fp.region.size | F1      |
| ------ | ------------------------------- | ------ | ------- | --- | ---- | ------- | --- | --- | ------- | --------- | ------------ | ------------ | ------- | --------------- | --------------- | --- | --------- | -------------- | ------- |
| test8  | test8.SEQC2.freebayes.stats.csv | indels | 1922    | 1   | 1921 | 256     | 255 | 0   | 0.00052 | 0.00391   | 0.0          | 0.00243      | 0.00052 | 0.0             | 0.01811         | 0.0 | 0.0       | 2485467324     | 0.00092 |
| test9  | test9.SEQC2.mutect2.stats.csv   | indels | 1922    | 1   | 1921 | 1       | 0   | 0   | 0.00052 | 1.0       | 0.0          | 0.00243      | 0.00052 | 0.025           | 1.0             | 0.0 | 0.0       | 2485467324     | 0.00104 |
| test10 | test10.SEQC2.strelka.stats.csv  | indels | 1922    | 70  | 1852 | 445     | 375 | 0   | 0.03642 | 0.1573    | 0.02873      | 0.0455       | 0.03642 | 0.12574         | 0.1933          | 0.0 | 0.0       | 2485467324     | 0.05915 |

_RTGtools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP   | FN   | Precision | Recall | F1     |
| ------ | --------------------------------- | --------- | ------- | ------- | ---- | ---- | --------- | ------ | ------ |
| test8  | test8.SEQC2.freebayes.summary.txt | 99.000    | 1       | 1       | 148  | 1623 | 0.0067    | 0.0006 | 0.0011 |
| test8  | test8.SEQC2.freebayes.summary.txt | None      | 1       | 1       | 249  | 1623 | 0.004     | 0.0006 | 0.0011 |
| test9  | test9.SEQC2.mutect2.summary.txt   | None      | 1       | 1       | 9    | 1623 | 0.1       | 0.0006 | 0.0012 |
| test10 | test10.SEQC2.strelka.summary.txt  | None      | 61      | 61      | 4593 | 1563 | 0.0131    | 0.0376 | 0.0194 |

As the truth set contains very few variants, the number of matcheds quite small.

## Case 10 : Somatic SV Variant benchmarking

- [test config](../conf/tests/somatic_sv.config)

### Analysis

- In this case we are using nf-core/sarek variant calls from manta and tiddit. AWS stored files are being pulled automatically with the pipeline.
- We are filtering out extra contigs, and taking only the variants longer than 30bp.
- We use truvari as benchmarking method and all the method related parameters are in samplesheet.
- Please note that SEQC2 SV truth file used in this case is not validated, it is generated only for testing purposes using the dataset reported in (Talsania et all.)[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02816-6]

### Results

_Truvari

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

| Tool   | File                            | Threshold | TP_base | TP_call | FP  | FN  | Precision | Recall | F1     |
| ------ | ------------------------------- | --------- | ------- | ------- | --- | --- | --------- | ------ | ------ |
| test11 | test11.SEQC2.manta.summary.txt  | None      | 28      | 28      | 156 | 334 | 0.1522    | 0.0773 | 0.1026 |
| test14 | test14.SEQC2.cnvkit.summary.txt | None      | 151     | 151     | 916 | 211 | 0.1415    | 0.4171 | 0.2113 |

_Truvari_

| Tool   | File                             | Caller | TP_base | TP_call | FP  | FN   | Precision | Recall  | F1      |
| ------ | -------------------------------- | ------ | ------- | ------- | --- | ---- | --------- | ------- | ------- |
| test11 | test11.SEQC2.manta.summary.json  | manta  | 5       | 5       | 695 | 1114 | 0.00714   | 0.00447 | 0.0055  |
| test14 | test14.SEQC2.cnvkit.summary.json | cnvkit | 8       | 8       | 96  | 1111 | 0.07692   | 0.00715 | 0.01308 |

_Intersect_

| Tool             | File                                    | Caller       | Analysis           | TP_base | TP_comp | FP   | FN   | Precision | Recall  | F1      |
| ---------------- | --------------------------------------- | ------------ | ------------------ | ------- | ------- | ---- | ---- | --------- | ------- | ------- |
| test11_converted | test11.SEQC2.manta.converted_stats.csv  | manta        | Bedtools_Intersect | 1       | 1       | 1063 | 1501 | 0.00094   | 0.00067 | 0.00078 |
| test13           | test13.SEQC2.controlfreec_stats.csv     | controlfreec | Bedtools_Intersect | 95      | 95      | 344  | 1407 | 0.2164    | 0.06325 | 0.09789 |
| test14           | test14.SEQC2.cnvkit_stats.csv           | cnvkit       | Bedtools_Intersect | 123     | 123     | 259  | 1379 | 0.32199   | 0.08189 | 0.13057 |
| test14_converted | test14.SEQC2.cnvkit.converted_stats.csv | cnvkit       | Bedtools_Intersect | 90      | 90      | 207  | 1412 | 0.30303   | 0.05992 | 0.10006 |
| test15           | test15.SEQC2.ascat_stats.csv            | ascat        | Bedtools_Intersect | 64      | 64      | 119  | 1438 | 0.34973   | 0.04261 | 0.07596 |

Only CNVkit was reporting the variants in VCF format that is why truvari and rgtools was running only that test sample.
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

| Tool             | File                         | Type  | TP    | FP  | FN  | Precision | Recall | F1      |
| ---------------- | ---------------------------- | ----- | ----- | --- | --- | --------- | ------ | ------- |
| test6-test7      | test6-test7.summary.tsv      | INDEL | 7784  | 184 | 80  | 0.99      | 0.977  | 0.98346 |
| test6-test7      | test6-test7.summary.tsv      | SNP   | 44755 | 20  | 962 | 0.979     | 1.0    | 0.98939 |
| test6-test7-copy | test6-test7-copy.summary.tsv | INDEL | 7784  | 184 | 80  | 0.99      | 0.977  | 0.98346 |
| test6-test7-copy | test6-test7-copy.summary.tsv | SNP   | 44755 | 20  | 962 | 0.979     | 1.0    | 0.98939 |
| test7-test7-copy | test7-test7-copy.summary.tsv | INDEL | 7864  | 0   | 0   | 1.0       | 1.0    | 1.0     |
| test7-test7-copy | test7-test7-copy.summary.tsv | SNP   | 45717 | 0   | 0   | 1.0       | 1.0    | 1.0     |

Concordance statistics are based on tool comparisons, the pairs are generated randomly for each run therefore the order of the comparison may change. Here, test6 is considered as BASE whole test7 is being considered as COMP set. Therefore, be aware that FP, FN and Precisison and Recall may switch. test-7 vs a copy of test-7 are basically the same dataset compared each other.

## Case 13 : Ensembling test files to generate a proxy truth for small variants

Truth/Golden set of samples to use as a baseline for benchmarking analysis is a challenge especially for somatic studies. One of the common approaches reseacrhes uses is to generate a proxy set of truth using input (test) VCF files, first ensembling the variants and then filtering them down with a majority vote system.

### Config file

- [test config](../conf/tests/somatic_snv_ensemble.config)

### Analysis

- Benchmarking somatic variants is possible using this pipeline. In order to demonstare SNV type of benchmarking, 3 example somatic variant calls will be used.
- Truth generation is done accordingly: ensemble_truth is used as 2, so the variant is kept only if it exist in 2 or more callers (test vcfs)
- Preprocessing includes normalization enabling left alignment of the variants.
- Using som.py and rtgtools vcfeval.

### Results

_Som.py_

| Tool   | File                            | Type | TP_base | TP  | FN  | TP_call | FP   | UNK | Recall  | Precision | recall_lower | recall_upper | recall2 | precision_lower | precision_upper | na  | ambiguous | fp.region.size | F1      |
| ------ | ------------------------------- | ---- | ------- | --- | --- | ------- | ---- | --- | ------- | --------- | ------------ | ------------ | ------- | --------------- | --------------- | --- | --------- | -------------- | ------- |
| test8  | test8.truth.freebayes.stats.csv | SNVs | 7       | 7   | 0   | 1898    | 1891 | 0   | 1.0     | 0.00369   | 0.59038      | 1.0          | 1.0     | 0.00165         | 0.00723         | 0.0 | 0.0       | 46709983       | 0.00735 |
| test9  | test9.truth.mutect2.stats.csv   | SNVs | 7       | 6   | 1   | 8       | 2    | 0   | 0.85714 | 0.75      | 0.4992       | 1.0          | 0.85714 | 0.40838         | 0.94403         | 0.0 | 0.0       | 46709983       | 0.8     |
| test10 | test10.truth.strelka.stats.csv  | SNVs | 7       | 7   | 0   | 5251    | 5244 | 0   | 1.0     | 0.00133   | 0.59038      | 1.0          | 1.0     | 0.0006          | 0.00262         | 0.0 | 0.0       | 46709983       | 0.00266 |

_RTG Tools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP     | FN  | Precision | Recall | F1     |
| ------ | --------------------------------- | --------- | ------- | ------- | ------ | --- | --------- | ------ | ------ |
| test8  | test8.truth.freebayes.summary.txt | 0.000     | 152     | 152     | 1746   | 4   | 0.0801    | 0.9744 | 0.148  |
| test8  | test8.truth.freebayes.summary.txt | None      | 152     | 152     | 1746   | 4   | 0.0801    | 0.9744 | 0.148  |
| test9  | test9.truth.mutect2.summary.txt   | None      | 18      | 18      | 5      | 138 | 0.7826    | 0.1154 | 0.2011 |
| test10 | test10.truth.strelka.summary.txt  | None      | 154     | 154     | 112831 | 2   | 0.0014    | 0.9872 | 0.0027 |

The contrast in results stems from how each tool handles input filtering and regional masking. som.py is built specifically for somatic data, applying strict regional masking and aggressively filtering the baseline down to just 7 variants. RTG Tools lacks this somatic aware pre-filtering, evaluating a much broader and noisier baseline of 156 variants. Consequently, RTG Tools reports massively inflated false positives and false negatives because it evaluates a wider genomic footprint and struggles to reconcile complex somatic alleles without a specialized matching engine.

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

| Tool   | File                             | Caller | TP_base | TP_call | FP   | FN   | Precision | Recall  | F1      |
| ------ | -------------------------------- | ------ | ------- | ------- | ---- | ---- | --------- | ------- | ------- |
| test11 | test11.truth.manta.summary.json  | manta  | 333     | 333     | 631  | 6447 | 0.34544   | 0.04912 | 0.086   |
| test12 | test12.truth.tiddit.summary.json | tiddit | 3344    | 3344    | 9013 | 3436 | 0.27062   | 0.49322 | 0.34948 |

> [!NOTE]
> As a final note, it is highly encouraged to run all those test and investigate comparisons, plots and tables created. This pipeline serves for multiple tools enabling efficient benchmark analysis.
