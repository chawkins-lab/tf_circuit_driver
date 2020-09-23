# tf_circuit_driver

```{sh}

----------
Requirements

python >=v3.5, click >=7.0

----------
Installation

cd /path/to/this/repo
pip3 install --editable .

----------
Usage
tf_circuit_driver --actions /path/to/species.protein.actions.v11.0.txt.gz \
--mappings /path/to/species.protein.info.v11.0.txt.gz \
--targets /path/to/targets_0.0/ \
--dge /path/to/sig_dge_mapped.csv \
--output /path/to/test.csv

Mappings: download string PPI mappings for mouse (https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz)
or human (https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz)

Actions: download PPI actions for mouse (https://stringdb-static.org/download/protein.actions.v11.0/10090.protein.actions.v11.0.txt.gz)
or human (https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz)

Targets: run IS-MARA (https://ismara.unibas.ch/mara/) by uploading your fastq
files and run for mouse (mm10) or human (hg19) as appropriate. Once run is
completed and you receive the results email, download ismara_report.tar.gz
using the link provided and decompress to targets_0.0 (this name is the default,
  so make sure you put separate MARA runs in different directories to avoid
  overwriting existing data)

DGE: run standard differential expression with e.g. DESeq2, edgeR. The csv input
file supplied as input should have 3 columns, in this order, headed "symbol",
"fc", "pvalue", where:
    -symbol: HGNC symbol for mouse or human gene
    -fc: log2FC in treatment/control
    -pvalue: adjusted p-value for treatment/control

----------
Example data

Fastq files for siCTR and siMYC transfected A549 cells were downloaded from
GSE112190 (https://mcr.aacrjournals.org/content/molcanres/early/2019/12/20/1541-7786.MCR-19-0657.full.pdf).
After quality trimming, reads were aligned with STAR to GRCh37 and duplicates
marked with Picard. Counts were generated with ht-seq and differential
expression tested with DESeq to provide input for file myc.dge.txt. Fastq were
also processed with IS-MARA and output file (targets.tar.gz) was decompressed.
To run tf_circuit_driver decompress example_data.tar.gz and run"

cd /path/to/example_data
tf_circuit_driver \
    --actions /path/to/9606.protein.actions.v11.0.txt.gz \
    --mappings /path/to/9606.protein.info.v11.0.txt.gz \
    --targets targets_0.0 \
    --dge example.dge.txt \
    --output example.out.csv

```
