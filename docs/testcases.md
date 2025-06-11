# nf-core/variantbenchmarking: Test cases

This pipeline is able to benchmark various type of analysis. Below, explanations of some use cases given in [tests](../conf/tests/) will be explained.

## Case 1 : Germline Structural Variant benchmarking

### Config file

- [test config](../conf/tests/germline_sv.config)

### Analysis

- We are using 3 different publicly available structural variant calls for HG002 to benchmark against Genome in a Bottle HG002 SV analysis. Only chromosome 21 will be used for the analysis.
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
| test6 | test6.HG002.strelka.summary.txt  | 24.000    | 51964   | 51989   | 301  | 303 | 0.9942    | 0.9942 | 0.9942 |
| test6 | test6.HG002.strelka.summary.txt  | None      | 52028   | 52053   | 449  | 239 | 0.9914    | 0.9954 | 0.9934 |
| test7 | test7.HG002.bcftools.summary.txt | None      | 51650   | 51675   | 1060 | 617 | 0.9799    | 0.9882 | 0.984  |

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

| Tool   | File                            | Type | TP_base | TP  | FN  | TP_call | FP   | UNK | Recall             | Precision          | recall_lower       | recall_upper       | recall2            | precision_lower    | precision_upper    | na  | ambiguous | fp.region.size | F1                 |
| ------ | ------------------------------- | ---- | ------- | --- | --- | ------- | ---- | --- | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | --- | --------- | -------------- | ------------------ |
| test8  | test8.SEQC2.freebayes.stats.csv | SNVs | 473     | 6   | 467 | 1719    | 1713 | 0   | 0.0126849894291754 | 0.0034904013961605 | 0.0053115241334789 | 0.0259581656571561 | 0.0126849894291754 | 0.0014581505501261 | 0.0071804021927946 | 0.0 | 0.0       | 46709983.0     | 36.673102621339    |
| test9  | test9.SEQC2.mutect2.stats.csv   | SNVs | 473     | 6   | 467 | 8       | 2    | 0   | 0.0126849894291754 | 0.75               | 0.0053115241334789 | 0.0259581656571561 | 0.0126849894291754 | 0.4083764894266625 | 0.9440331102753412 | 0.0 | 0.0       | 46709983.0     | 0.042817399441143  |
| test10 | test10.SEQC2.strelka.stats.csv  | SNVs | 473     | 9   | 464 | 5251    | 5242 | 0   | 0.0190274841437632 | 0.0017139592458579 | 0.009455525942495  | 0.0344370723145163 | 0.0190274841437632 | 0.0008484051549055 | 0.0031258365589098 | 0.0 | 0.0       | 46709983.0     | 112.22440393523586 |

_RTGtools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP   | FN  | Precision | Recall | F1     |
| ------ | --------------------------------- | --------- | ------- | ------- | ---- | --- | --------- | ------ | ------ |
| test8  | test8.SEQC2.freebayes.summary.txt | 65.000    | 8       | 7       | 1456 | 551 | 0.0048    | 0.0143 | 0.0072 |
| test8  | test8.SEQC2.freebayes.summary.txt | None      | 8       | 7       | 1949 | 551 | 0.0036    | 0.0143 | 0.0057 |
| test9  | test9.SEQC2.mutect2.summary.txt   | None      | 8       | 7       | 27   | 551 | 0.2059    | 0.0143 | 0.0268 |
| test10 | test10.SEQC2.strelka.summary.txt  | None      | 0       | 0       | 0    | 559 |           | 0.0    |        |

The number of TPs found quite similar with som.py and RTGtools, yet test10 was not able to run properly with rtgtools as GT field is not reported in strelka.

## Case 9 : Somatic INDEL Variant benchmarking

- [test config](../conf/tests/somatic_indel.config)

### Analysis

- In order to demonstare INDEL type of benchmarking for somatic variants, 3 example somatic variant calls will be benchmarked against SEQC2 truth INDEL calls.
- This time we are not applying preprocessing to the test variants.
- Using som.py (a version of hap.py tuned specially for somatic benchmarking) and rtgtools vcfeval with "--squash-ploidy" parameter on.

### Results

_SOM.py_

| Tool   | File                            | Type   | TP_base | TP  | FN  | TP_call | FP  | UNK | Recall             | Precision          | recall_lower | recall_upper       | recall2            | precision_lower | precision_upper    | na  | ambiguous | fp.region.size | F1                 |
| ------ | ------------------------------- | ------ | ------- | --- | --- | ------- | --- | --- | ------------------ | ------------------ | ------------ | ------------------ | ------------------ | --------------- | ------------------ | --- | --------- | -------------- | ------------------ |
| test8  | test8.SEQC2.freebayes.stats.csv | indels | 22      | 0   | 22  | 504     | 504 | 0   | 0.0                | 0.0                | 0.0          | 0.1543725128155745 | 0.0                | 0.0             | 0.0072924851130725 | 0.0 | 0.0       | 46709983.0     | 10.789984659168042 |
| test9  | test9.SEQC2.mutect2.stats.csv   | indels | 22      | 1   | 21  | 1       | 0   | 0   | 0.0454545454545454 | 1.0                | 0.0          | 0.1934373646833098 | 0.0454545454545454 | 0.025           | 1.0                | 0.0 | 0.0       | 46709983.0     | 0.0                |
| test10 | test10.SEQC2.strelka.stats.csv  | indels | 22      | 1   | 21  | 445     | 444 | 0   | 0.0454545454545454 | 0.0022471910112359 | 0.0          | 0.1934373646833098 | 0.0454545454545454 | 0.0             | 0.0104547022441303 | 0.0 | 0.0       | 46709983.0     | 9.505462675933751  |

_RTGtools_

| Tool   | File                              | Threshold | TP_base | TP_call | FP   | FN  | Precision | Recall | F1  |
| ------ | --------------------------------- | --------- | ------- | ------- | ---- | --- | --------- | ------ | --- |
| test8  | test8.SEQC2.freebayes.summary.txt | 0.000     | 0       | 0       | 2353 | 44  | 0.0       | 0.0    | 0.0 |
| test8  | test8.SEQC2.freebayes.summary.txt | None      | 0       | 0       | 2353 | 44  | 0.0       | 0.0    | 0.0 |
| test9  | test9.SEQC2.mutect2.summary.txt   | None      | 0       | 0       | 0    | 44  |           | 0.0    |     |
| test10 | test10.SEQC2.strelka.summary.txt  | None      | 0       | 0       | 0    | 44  |           | 0.0    |     |

As the truth set contains very few variants, the number of matcheds quite small and rtgtools couldnt able to find any.

## Case 10 : Somatic SV Variant benchmarking

- [test config](../conf/tests/somatic_sv.config)

### Analysis

- In this case we are using nf-core/sarek variant calls from manta and tiddit. AWS stored files are being pulled automatically with the pipeline.
- We are filtering out extra contigs, and taking only the variants longer than 30bp.
- We use truvari as benchmarking method and all the method related parameters are in samplesheet.
- Please note that SEQC2 SV truth file used in this case is not validated, it is generated only for testing purposes using the dataset reported in (Talsania et all.)[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02816-6]

### Results

| Tool   | File                             | TP_base | TP_comp | FP  | FN  | Precision           | Recall               | F1                  |
| ------ | -------------------------------- | ------- | ------- | --- | --- | ------------------- | -------------------- | ------------------- |
| test11 | test11.SEQC2.manta.summary.json  | 3       | 3       | 5   | 591 | 0.375               | 0.005050505050505001 | 0.009966777408637   |
| test12 | test12.SEQC2.tiddit.summary.json | 29      | 29      | 84  | 565 | 0.25663716814159204 | 0.048821548821548    | 0.08203677510608201 |

## Case 11 : Somatic CNV Variant benchmarking

- [test config](../conf/tests/somatic_cnv.config)

### Analysis

- We are again using nf-core/sarek variant calls from cnvkit, ascat and controlfreec.
- We are filtering out extra contigs and splitting multi-allelic sites.
- We use truvari, wittyer and _intersect_ methods. _intersect_ (using bedtools) is not a benchmarking method but rather to intersect BED regions given for test and truth files. truvari and wittyer will be only applied if input test is in VCF format reming that CNV results are tend to report in BED similar formats.
- The truth file used here is reported previously through (zenodo)[https://zenodo.org/records/14619054]

### Results

_Truvari_

| Tool   | File                             | TP_base | TP_comp | FP  | FN   | Precision           | Recall            | F1                   |
| ------ | -------------------------------- | ------- | ------- | --- | ---- | ------------------- | ----------------- | -------------------- |
| test14 | test14.SEQC2.cnvkit.summary.json | 26      | 26      | 84  | 1058 | 0.23636363636363603 | 0.023985239852398 | 0.043551088777219006 |

_Witty.er_

| Tool   | File                     | StatsType | TP_base | TP_comp | FP        | FN  | Precision | Recall | F1  |
| ------ | ------------------------ | --------- | ------- | ------- | --------- | --- | --------- | ------ | --- |
| test14 | test14.SEQC2.cnvkit.json | Event     | 0       | 0       | 169       | 0   | 0.0       |        |     |
| test14 | test14.SEQC2.cnvkit.json | Base      | 0       | 0       | 893824341 | 0   | 0.0       |        |     |

_Intersect_

| Tool   | File                                    | TP_base | TP_comp | FN  | FP  | Precision          | Recall             | F1                 |
| ------ | --------------------------------------- | ------- | ------- | --- | --- | ------------------ | ------------------ | ------------------ |
| test13 | test13.SEQC2.controlfreec_stats.csv     | 123     | 108     | 537 | 331 | 0.2460136674259681 | 0.1863636363636363 | 0.2120740439186762 |
| test14 | test14.SEQC2.cnvkit_stats.csv           | 144     | 128     | 516 | 246 | 0.3422459893048128 | 0.2181818181818181 | 0.2664816099930603 |
| test14 | test14.SEQC2.cnvkit.converted_stats.csv | 116     | 107     | 544 | 192 | 0.3578595317725752 | 0.1757575757575757 | 0.2357365342247208 |
| test15 | test15.SEQC2.ascat_stats.csv            | 92      | 81      | 568 | 77  | 0.5126582278481012 | 0.1393939393939394 | 0.21918936408024   |

Only CNVkit was reporting the variants in VCF format that is why truvari and wittyer was running only that test sample. Wittyer was not proprely running as variants are not sequence resolved.
We also have two type of results for CNVkit for intersection analysis. \*\_converted one show the result from the VCF input being converted to BED for ntersection analysis while the other result is using input BED file directly.

As a final note, it is highly encourged to run all those test and investigate comparions, plots and tables created. This pipeline serves for multiple tools enabling efficient benchmark analysis.
