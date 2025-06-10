# nf-core/variantbenchmarking: Test cases

This pipeline is able to benchmark various type of analysis. Below, explanations of some use cases given in [tests](../conf/tests/) will be explained.

## Case 1 : Germline Structural Variant benchmarking

### Config file

- [test config](../conf/tests/germline_sv.config)

### Analysis

- 3 different publicly available structural variant calls for HG002 are being benchmarked against Genome in a Bottle HG002 SV analysis. Only chromosome 21 will be used for the analysis.
- Preprocessing includes normalization enabling left alignment of the variants, splitting multi allelic variants, and deduplication. Using svync ensures correct formattings of the files.
- Size filtering by >30 bases applied to test files to ensure equal SV size in benchmarking.
- Filtering out BND and TRA type of variants
- Truvari, SVbenchmark and Wittyer are used as benchmarking method. Tools spesific parameters are given in the corresponding samplesheet.

### Results

_Truvari_

| Tool  | File                           | TP_base | TP_comp | FP  | FN   | Precision          | Recall            | F1                  |
| ----- | ------------------------------ | ------- | ------- | --- | ---- | ------------------ | ----------------- | ------------------- |
| test1 | test1.HG002.manta.summary.json | 211     | 211     | 9   | 1253 | 0.9590909090909091 | 0.144125683060109 | 0.250593824228028   |
| test2 | test2.HG002.lumpy.summary.json | 51      | 51      | 7   | 1413 | 0.8793103448275861 | 0.03483606557377  | 0.06701708278580801 |
| test3 | test3.HG002.delly.summary.json | 156     | 156     | 45  | 1308 | 0.7761194029850741 | 0.10655737704918  | 0.18738738738738703 |

_SVbenchmark_

| Tool  | File                     | TP_base | FP  | TP_comp | FN   | Recall | Precision | F1                 |
| ----- | ------------------------ | ------- | --- | ------- | ---- | ------ | --------- | ------------------ |
| test1 | test1.HG002.manta.report | 150     | 77  | 150     | 1314 | 0.1025 | 0.65      | 0.177015250544662  |
| test2 | test2.HG002.lumpy.report | 42      | 17  | 42      | 1422 | 0.0287 | 0.7069    | 0.0551392891450528 |
| test3 | test3.HG002.delly.report | 98      | 113 | 98      | 1366 | 0.0669 | 0.4593    | 0.116850694918833  |

_Wittyer_

| Tool  | File                   | StatsType | TP_base | TP_comp | FP       | FN     | Precision           | Recall            | F1                  |
| ----- | ---------------------- | --------- | ------- | ------- | -------- | ------ | ------------------- | ----------------- | ------------------- |
| test1 | test1.HG002.manta.json | Event     | 151     | 141     | 79       | 1313   | 0.64090909090909    | 0.103142076502732 | 0.177688571380881   |
| test1 | test1.HG002.manta.json | Base      | 31815   | 31886   | 40428    | 363698 | 0.440938130929004   | 0.080439833836055 | 0.136058646053033   |
| test2 | test2.HG002.lumpy.json | Event     | 44      | 42      | 16       | 1420   | 0.7241379310344821  | 0.030054644808743 | 0.057713928794503   |
| test2 | test2.HG002.lumpy.json | Base      | 28640   | 28640   | 31970    | 366873 | 0.47252928559643603 | 0.072412284804797 | 0.125580161491527   |
| test3 | test3.HG002.delly.json | Event     | 100     | 96      | 113      | 1364   | 0.45933014354066903 | 0.068306010928961 | 0.11892668665295701 |
| test3 | test3.HG002.delly.json | Base      | 63683   | 63683   | 32721459 | 331830 | 0.001942434777314   | 0.161013670852791 | 0.00383856195726    |

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
| test3 | test3.HG002.delly.summary.txt | None      | 52      | 54      | 373  | 406 | 0.1265    | 0.1135 | 0.1197 |
| test4 | test4.HG002.svaba.summary.txt | None      | 22      | 22      | 384  | 436 | 0.0542    | 0.048  | 0.0509 |

Note that RTGtools bndeval is only configured to benchmark SVTYPE=BND. Eventough we decomposed SV tpes to BND, SVs with other types can remain and they might be matched with the tool. That is why both precison and recall values are signififcanly lower than SVBenchmark and Truvari. bndeval should only be used for BND SV analysis with corresponding truths.

## Case 3 : Lifting over one test file for Structural Variant benchmarking

### Config file

- [test config](../conf/tests/lintover_test.config)

### Analysis

- Besides to the 4 test cases, we are adding another test case here whom is generated through GRCh38 reference instead of GRCh37. This test case will be lifted over to GRCh37 by the pipeline.
- Preprocessing includes normalization enabling left alignment of the variants, splitting multi allelic variants, and deduplication.
- Size filtering by >30 bases applied to test files to ensure equal SV size in benchmarking.
- Truvari is used as benchmarking method. Tools spesific parameters are given in the corresponding samplesheet.

### Results

| Tool  | File                            | TP_base | TP_comp | FP  | FN  | Precision | Recall | F1     |
| ----- | ------------------------------- | ------- | ------- | --- | --- | --------- | ------ | ------ |
| test1 | test1.HG002.manta.summary.json  | 118     | 118     | 78  | 238 | 0.6020    | 0.3315 | 0.4275 |
| test2 | test2.HG002.lumpy.summary.json  | 40      | 40      | 49  | 316 | 0.4494    | 0.1124 | 0.1798 |
| test3 | test3.HG002.delly.summary.json  | 34      | 34      | 44  | 322 | 0.4359    | 0.0955 | 0.1567 |
| test4 | test4.HG002.svaba.summary.json  | 50      | 50      | 245 | 306 | 0.1695    | 0.1404 | 0.1536 |
| test5 | test5.HG002.dragen.summary.json | 15      | 15      | 12  | 341 | 0.5556    | 0.0421 | 0.0783 |

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
| test6 | test6.HG002.strelka.summary.txt  | 30.000    | 52163   | 52187   | 405  | 569 | 0.9923    | 0.9892 | 0.9908 |
| test6 | test6.HG002.strelka.summary.txt  | None      | 52281   | 52305   | 708  | 451 | 0.9866    | 0.9914 | 0.989  |
| test7 | test7.HG002.bcftools.summary.txt | None      | 51920   | 51947   | 1270 | 812 | 0.9761    | 0.9846 | 0.9804 |

While hap.py reports the metrics SNP and INDELs seperated, RTGtools provides all merged. As you can see, even for small variant benchmarking there are differences with the results.

## Case 5 : Lifting over truth vcf for Germline Small Variant benchmarking

### Config file

- [test config](../conf/tests/liftover_truth.config)

### Analysis

- We are using the same test cases used in _Case 4_ but in this case, we will try to becnhamark them against a truth file which is generated with GRCh37 using lifting over strategy.
- Both RTGtools vcfeval and hap.py methods are used for benchmarking.

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
| test6 | test6.HG002.strelka.summary.txt  | 24.000    | 51964   | 51989   | 301  | 303 | 0.9942    | 0.9942 | 0.9942 |
| test6 | test6.HG002.strelka.summary.txt  | None      | 52028   | 52053   | 449  | 239 | 0.9914    | 0.9954 | 0.9934 |
| test7 | test7.HG002.bcftools.summary.txt | None      | 51650   | 51675   | 1060 | 617 | 0.9799    | 0.9882 | 0.984  |

As you can see, lifting over files from GRCh37 to GRCh38 worked quite well, there is not too much difference.

## Case 6 : GH4GH Germline Small Variant benchmarking

## Case 7 : Germline Small Variant benchmarking with stratifications

## Case 8 : Somatic SNV Variant benchmarking

## Case 9 : Somatic INDEL Variant benchmarking

## Case 10 : Somatic SV Variant benchmarking

## Case 11 : Somatic CNV Variant benchmarking
