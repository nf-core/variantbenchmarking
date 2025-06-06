# nf-core/variantbenchmarking: Test cases

This pipeline is able to benchmark various type of analysis. Below, explanations of some use cases given in [tests](../conf/tests/) will be explained.

## Case 1 : Germline Structural variant benchmarking

### Config file

- [test config](../conf/tests/germline_sv.config)

### Analysis

- 4 different publicly available structural variant calls for HG002 are being benchmarked against Genome in a Bottle HG002 SV analysis. Only chromosome 21 will be used for the analysis.
- Preprocessing includes normalization enabling left alignment of the variants, splitting multi allelic variants, and deduplication. Using svync ensures correct formattings of the files.
- Size filtering by >30 bases applied to test files to ensure equal SV size in becnhmarking.
- Truvari and SVbenchmark used as benchmarking method. Tools spesific parameters are given in the corresponding samplesheet.

### Results

_Truvari_

| Tool  | File                           | TP_base | TP_comp | FP  | FN   | Precision | Recall | F1     |
| ----- | ------------------------------ | ------- | ------- | --- | ---- | --------- | ------ | ------ |
| test1 | test1.HG002.manta.summary.json | 404     | 404     | 56  | 1060 | 0.8783    | 0.2760 | 0.4200 |
| test2 | test2.HG002.lumpy.summary.json | 205     | 205     | 58  | 1259 | 0.7795    | 0.1400 | 0.2374 |
| test3 | test3.HG002.delly.summary.json | 157     | 157     | 171 | 1307 | 0.4787    | 0.1072 | 0.1752 |
| test4 | test4.HG002.svaba.summary.json | 286     | 286     | 120 | 1178 | 0.7044    | 0.1954 | 0.3059 |

_SVbenchmark_

| Tool  | File                     | TP_base | FP  | TP_comp | FN   | Recall | Precision | F1     |
| ----- | ------------------------ | ------- | --- | ------- | ---- | ------ | --------- | ------ |
| test1 | test1.HG002.manta.report | 150     | 317 | 150     | 1314 | 0.1025 | 0.3109    | 0.1541 |
| test2 | test2.HG002.lumpy.report | 42      | 222 | 42      | 1422 | 0.0287 | 0.1559    | 0.0485 |
| test3 | test3.HG002.delly.report | 98      | 240 | 98      | 1366 | 0.0669 | 0.2857    | 0.1085 |
| test4 | test4.HG002.svaba.report | 0       | 406 | 0       | 1464 | 0.0    | 0.0       | NA     |

SVBenchmark is not able to benchmark test4 properly since SVs are not sequence resolved. We made sure benchmarking unresolved SVs benchmark by pctseq = 0 in Truvari. Moreover, the number of TPs found in SVbenchmark is significanly lower than Truvari yet this is not primarly as methodological differences but also because of differences of the parameters defining SV comparions. Therefore, it is highly important to set meaningful parameters before starting to perform benchmarks.

## Case 2 : Germline BND (Breakends representation of structural variants) benchmarking

### Config file

- [test config](../conf/tests/germline_bnd.config)

### Analysis

- The same samples in case 1 are used with the same preprocessing.
- Size filtering by >30 bases applied to test files to ensure equal SV size in becnhmarking.
- svdecompose is used to convert SVTYPEs to "BND"
- RTGtools bndeval is used as benchmarking method.

### Results

| Tool  | Threshold | TP_base | TP_call | FP   | FN  | Precision | Recall | F1     |
| ----- | --------- | ------- | ------- | ---- | --- | --------- | ------ | ------ |
| test1 | None      | 163     | 165     | 1302 | 295 | 0.1125    | 0.3559 | 0.1709 |
| test2 | None      | 48      | 48      | 279  | 410 | 0.1468    | 0.1048 | 0.1223 |
| test3 | None      | 52      | 54      | 373  | 406 | 0.1265    | 0.1135 | 0.1197 |
| test4 | None      | 22      | 22      | 384  | 436 | 0.0542    | 0.048  | 0.0509 |

Note that RTGtools bndeval is only configured to benchmark SVTYPE=BND. Eventough we decomposed SV tpes to BND, SVs with other types can remain and they might be matched with the tool. That is why both precison and recall values are signififcanly lower than SVBenchmark and Truvari. bndeval should only be used for BND SV analysis with corresponding truths.
