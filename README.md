# branchwater-web-MAGtest

Associated Manuscript will be added once published.

Zenodo with required files (bins) and pre-computed results: 


In brief:
After search with branchwater-web, metagenome results from five genome queries went through stringent testing for verification of match accuracy by creating metagenome assembled genomes (MAGs) from the publicly available SRA metagenomes. For each query genome, the CSV of full results from branchwater-web were downloaded and filtered to cANI levels denoting high match quality and sufficient matches for testing (0.96-98). Metagenomes were downloaded with nf-core/fetchngs v. 1.10.0 (15) and processed into MAGs with nf-core/mag v. 2.3.2, which only kept MAG bins that meet high quality criteria. For metagenomes with bins, the resulting MAGs were searched against NCBI via both sourmash and BBTools/sendsketch (17) to identify potential matches at the genus or species level. Independent taxonomic classification was also carried out using BAT bin annotation, which searches predicted ORFs against NCBI non-redundant protein database (NR from 12/12/2024).

Here, we include the code and configuration used for each analysis.

1. Configuration for MAG generation can be found in the `nfcore` folder

2. The `bbtools sendsketch` workflow can be found in the `sendsketch` folder

3. The `BAT` annotation workflow is in the main directory and uses the `BAT.snakefile` workflow file. Instructions for running can be found at the top of this file.

4. The `sourmash` workflow and final aggregation across tools uses the `sourmash.snakefile` workflow file. Instructions for running can be found at the top of this file.


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17088951.svg)](https://doi.org/10.5281/zenodo.17088951)
