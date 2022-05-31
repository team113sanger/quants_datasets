# QUANTS example

Our example dataset is raw counts from the BRCA1 saturation genome editing (SGE) experiment by Findlay *et al.*. 

> Findlay, G.M., Daza, R.M., Martin, B. *et al.*   
> **Accurate classification of BRCA1 variants with saturation genome editing.**  
> *Nature 562, 217–222 (2018).*   
> [https://doi.org/10.1038/s41586-018-0461-z](https://doi.org/10.1038/s41586-018-0461-z)

This experiment uses a mapping-based approach for quantification. As QUANTS requires a library, we will use [VaLiAnT](https://github.com/cancerit/VaLiAnT) to generate a library corresponding to BRCA1 exon2:

> Barbon, L., Offord, V., Radford, E.J., *et al.*  
> **Variant Library Annotation Tool (VaLiAnT): an oligonucleotide library design and annotation tool for saturation genome editing and other deep mutational scanning experiments.**   
> *Bioinformatics 38(4), 892–899 (2022).*. 
> https://doi.org/10.1093/bioinformatics/btab776

For this example we will be quantifying multiple samples which correspond to BRCA1 exon 2:

* plasmid
* negative control
* day 5 (2 replicates)
* day 11 (2 replicates)

For more information on how to use QUANTS, please see the [documentation](https://github.com/cancerit/QUANTS/tree/master/docs).

## Dependencies and installations

For this example, we ran QUANTS on our LSF compute cluster with Singularity enabled. [QUANTS](https://github.com/cancerit/QUANTS) can be cloned directly from GitHub. QUANTS is a Nextflow pipeline which has been developed using the [nf-core](https://nf-co.re/) framework. This means it is easily portable. 

We recommend using QUANTS with Singularity or Docker enabled for reproducibility. In theory, Conda or local installations of tools is supported but has not yet been fully tested.

Minimum requirement to run QUANTS is an installaton of [nextflow](https://www.nextflow.io) and either [Docker](https://www.docker.com) or [Singularity](https://sylabs.io).

## QUANTS inputs

All of the inputs required to run QUANTS can be found in [example_input](example_input).

### Sample to sequencing data mapping

QUANTS accepts both FASTQ and CRAM formats for sequencing data. In this example, we downloaded the relevant sample FASTQ files from the SRA project [PRJNA481326](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA481326) which correspond to BRCA1 exon 2 day 5 and day 11 samples, plasmid library and negative control. For more information, see [prepare_quants_inputs.ipynb](prepare_quants_inputs.ipynb).

QUANTS requires a comma-separated samplesheet (`--input`) which contains the relevant information to map sample names to their corresponding sequencing file(s). This was manually generated for this example.

As the example dataset is paired end FASTQ, three columns with predefined headers (`sample`, `fastq_1`, `fastq_2`) are required. Here, we use relative paths but, it is recommended to use absolute paths when using Docker or Singularity.

| sample | fastq_1 | fastq_2 |
| --- | --- | --- |
| X2_lib | example_input/fastq/SRR7525820_1.fastq.gz | example_input/fastq/SRR7525820_2.fastq.gz | 

### Oligo library

QUANTS uses an oligo library for quantification of guide abundance. At Sanger, our SGE libraries are currently generated using  [VaLiAnT](https://github.com/cancerit/VaLiAnT). In this example, we provide the input data and commands to generate a library of oligos which correspond to the Findlay *et al.* BRCA1 X2 (exon 2) experiment, including the necessary PAM protection edits. The details for running VaLiAnT can be found in [example_input/valiant_library/README.md](example_input/valiant_library/README.md).

QUANTS quantifies unique read abundance and oligo abundance using [pyCROQUET](https://github.com/cancerit/pycroquet) which has a strict library format detailed [here](https://github.com/cancerit/pycroquet/wiki/Guide-library-format). Details of the conversion of the VaLiAnT library into pyCROQUET format can be found in [prepare_quants_inputs.ipynb](prepare_quants_inputs.ipynb).

### Configurations and parameters

In this example, we used a configuration file (`-c`) that was specific to running QUANTS on an LSF compute cluster with Singularity and a JSON-formatted parameter file (`-params-file`) which contains the pipeline-specific key value pairs.

* configuration file (`-c`) [example_input/lsf.config](example_input/lsf.config)
* parameter file (`-params-file`) [example_input/params.json](example_input/params.json)

The parameter file sets the following parameters:

* `outdir` - where to write pipeline results
* `single_end` - set to `false` as we have paired end data
* `input_type` - set to `fastq` as our sequencing is in FASTQ format
* `read_transform` - set to `null` as we don't need to reverse complement our reads to match the oligos
* `raw_sequencing_qc`, `read_merging_qc` - set to `true` so [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [SeqKit stats](https://bioinf.shenwei.me/seqkit/usage/) are run to give QC metrics/plots per stage
* `read_trimming`, `read_trimming_qc` - set to `null` and `false` respectively as we don't want to trim adapters as a separate step (to mirror Findlay)
* `read_merging` - set to `seqprep` to use [SeqPrep](https://github.com/jstjohn/SeqPrep) for merging paired reads into a single construct
* `seqprep_options` - mirroring SeqPrep command line options from the Findlay paper
* `read_filtering` - set to `null` as we don't want to do any length filtering to save compute
* `quantification` - set to `pycroquet` to use [pyCROQUET](https://github.com/cancerit/pycroquet) to determine guide abundance from the merged reads

For more information, see the [QUANTS documentation](https://github.com/cancerit/QUANTS/tree/master/docs).

## Running QUANTS with the example dataset

Assuming that you have cloned the example dataset repository, followed the steps in [prepare_quants_inputs.ipynb](prepare_quants_inputs.ipynb) and have Nextflow and Singularity installed:

```
QUANTS \
  -c example_input/lsf.config \
  -params-file example_input/params.json \
  --input example_input/samplesheet.csv \
  --oligo_library example_input/BRCA1_exon2_pyCROQUET.tsv
```

## QUANTS pipeline outputs

QUANTS results can be found in [example_output/results](example_output/results) where output files are organised by software application.

### Quality control

QUANTS can generate quality control metrics and plots by sample/stage or in a collated report.

* `fastqc` - FastQC HTML, tables and plots per sample, per stage (i.e. raw, merged...)
* `multiqc` - Compiled HTML report of all FastQC results in [example_output/results/multiqc/multiqc_report.html](example_output/results/multiqc/multiqc_report.html)
* `seqkit_stats` - SeqKit tabular results per sample, per stage (e.g. average read length, Q30 percentage...)

### Merged reads

In this example, we used SeqPrep to merge paired end reads into a single construct (equally, we could have used Flash2 which is also supported in QUANTS).

* `seqprep` - alignments (`*.merged.aln.gz`) and merged reads per sample (`*.merged.fq.gz`)

### Quantification

QUANTS uses pyCROQUET to determine the abundance of unique read sequences and library oligos. 

* `example_output/results/pycroquet/*_merged.cram*` - read to oligo alignments in CRAM format
* `example_output/results/pycroquet/*_merged.counts.tsv.gz` - oligo counts (generated by pyCROQUET)
* `example_output/results/pycroquet/*_merged.query_counts.tsv.gz` - unique read counts (generated by pyCROQUET)
* `example_output/results/pycroquet/*_merged.query_to_library_counts.tsv.gz` - oligo counts (generated by comparing oligos to unique read sequences and taking their counts)
* `example_output/results/pycroquet/*_merged.stats.json` - summary statistics based on `*_merged.counts.tsv.gz`

## Comparing raw counts

For the results of comparing QUANTS and Findlay *et al.* raw counts, please see [comparing_raw_counts.ipynb](comparing_raw_counts.ipynb).