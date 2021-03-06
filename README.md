# rna_seq_pipeline

> This is a customized script intergrating tools for reference-depdendent RNA seq data analysis. Users can modify the config file to suit for their own demands.

### Format

- Input: fastq(clean)
- Output: expression matrix

### Features

- One-stop automatic pipeline for rna-seq data analysis
- Parallel processing

### requirements

- samtools 1.9
- hisat2 2.1.0
- featureCounts 1.6.4 (subread 1.6.4)
- HTSeq 0.11.2
- umi_tools 0.5.5 (python2)
- umi_tools 1.0.0 (python3)

### Usage

```shell
# For single-end
python rna_seq_pipeline.py -c SE.conf
# For pair-end
python rna_seq_pipeline.py -c PE.conf
```

