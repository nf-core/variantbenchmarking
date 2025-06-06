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
- Size filtering by >30 bases applied to test files to ensure equal SV size in becnhmarking.
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

## Case 3 : Linting over one test file for Structural Variant benchmarking

### Config file

- [test config](../conf/tests/lintover_test.config)

### Analysis

- Besides to the 4 test cases, we are adding another test case here whom is generated through GRCh38 reference instead of GRCh37. This test case will be lifted over to GRCh37 by the pipeline.
- Preprocessing includes normalization enabling left alignment of the variants, splitting multi allelic variants, and deduplication.
- Size filtering by >30 bases applied to test files to ensure equal SV size in becnhmarking.
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
