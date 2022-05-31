# VaLiAnT library for BRCA1 exon2

[VaLiAnT](https://github.com/cancerit/VaLiAnT) version 2.0.0 was used to generate a library of sequential single nucleotide variants across the CDS of BRCA1 exon 2 (ENST00000357654.9 chr17:43124017-43124096).

```
valiant sge \
    parameter_input_files/brca1_nuc_targeton_input.txt \
    reference_input_files/chr17.fa \
    brca1_nuc_output \
    'homo sapiens' \
    'GRCh38' \
    --pam parameter_input_files/brca1_protection_edits.vcf \
    --revcomp-minus-strand \
    --gff reference_input_files/ENST00000357654.9.gtf
```