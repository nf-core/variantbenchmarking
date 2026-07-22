The pipeline overview metro map is generated from docs/images/variantbenchmarking.svg using nf-metro version 1.1.0. If you add or rename pipeline steps, update the .mmd source and regenerate the images:

```
pip install 'nf-metro>=1.1.0'

nf-metro render docs/images/variantbenchmarking.mmd \
    -o docs/images/variantbenchmarking.svg \
    --theme nfcore \
    --no-chrome-css \
    --logo docs/images/nf-core-variantbenchmarking_logo_dark.png \
    --animate

```
