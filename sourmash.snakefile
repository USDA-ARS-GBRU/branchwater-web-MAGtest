# Snakefile for running sourmash analyses on a set of MAG bins and search genomes
# This is part of the 2025 BR-MAG test analysis, and is designed to be run after the bat.snakefile
# It sketches a set of search genomes, sketches the MAG bins, compares them using sourmash manysearch and multisearch,
# and aggregates the results with expected taxonomy and sendsketch results.
# It also aggregates the results from both the bat.snakefile and this sourmash.snakefile into a single summary file.

# You can download pre-computed results from Zenodo:
# Zenodo:

# Instructions to run:
    # 1. Make sure you're in the `branchwater-web-MAGtest` directory and have the the `branchwater-web-MAGtest` conda environment activated.
    #    If you need to install it, run: `conda env create -f environment.yaml`, then `conda activate branchwater-web-MAGtest`.

    # 2. Download metagenome bins from the associated Zenodo repository:
        # After downloading mulit-bins.tar.gz, untar via `tar -xvzf multi-bins.tar.gz`.

    # 3. Download metagenome signatures from the associated Zenodo repository:
        # After downloading `metagenome-sigs.tar.gz`, untar via `tar -xvzf metagenome-sigs.tar.gz`.
        # The file `br-magtest.wort.mf.csv` (included in the github) will look for the metagenome signatures in the `metagenomes/` directory.

    # 4: Download the srcipted sourmash NCBI 'entire' database: `ncbi-euks-all-2025.01-k51.dna.rocksdb.tar.gz` from https://sourmash.readthedocs.io/en/latest/databases-md/ncbi_euks_2025_01.html. Modify the database path below as needed to reflect your path to this database. 
        # command:
        # `curl -O --no-clobber https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/ncbi-euks-2025.01/ncbi-euks-all-2025.01-k51.dna.rocksdb.tar.gz`

        #For higher resulution searches of eukaryotes, you can also download the NCBI fungi, euk-other, and bilateria-minus-verts databases (see additional paths below).

    # 5. Make sure you've run the bat.snakefile to generate the annotated bins. Alternatively, download and untar the bat outputs (`bat-outputs.tar.gz`; `tar -xvzf bat-outputs.tar.gz`).

    # 6. Run the workflow with:
    #
    #    snakemake -s sourmash.snakefile --use-conda --cores 20 --rerun-incomplete --keep-going --latency-wait 60 --printshellcmds -n
    #
    # Note: the `-n` flag will do a dry-run; remove it to actually execute the workflow.

    # 7. The final aggregated results will be in `output.sourmash/multi.aggregated-results.csv` and `output.sourmash/multi.aggregated-results-bins.csv`.

##################################################
######### MODIFY THESE PATHS AS NEEDED ###########
db_ncbi_entire= "/group/ctbrowngrp5/sourmash-db/entire-2025-01-21/entire-2025-01-21.k51.rocksdb"

# Optional: for higher resolution eukaryote searches, you can also download scaled=1000 databases and modify the paths below.
# You would then need to modify the rule all and rule `aggregate_results` to use these higher resolution eukaryote multisearch results.
# (note: db_bilateria_minus_verts is quite large, ~50G)
db_fungi="/group/ctbrowngrp5/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k51.zip",
db_euk_other="/group/ctbrowngrp5/2025-genbank-eukaryotes/eukaryotes-other.sig.zip",
db_bilateria_minus_verts="/group/ctbrowngrp5/2025-genbank-eukaryotes/bilateria-minus-vertebrates.sig.zip",


###################################################
######## NO NEED TO MODIFY BELOW THIS LINE ########
outdir = "output.sourmash"
metagenome_mf = "br-magtest.wort.mf.csv",
bins_file = "multi_bins.txt"
expected_taxonomy = "multi_mapping.csv"
sendsketch_results = "multispecies_tracking.csv"


search_genomes = {"GCF_000888735.1": "Acanthamoeba polyphaga mimivirus",
                  "GCF_001890705.1": "Aspergillus sydowii",
                  "GCF_003013715.1": "Candida auris",
                  "GCF_000165345.1": "Cryptosporidium parvum ",
                  "GCF_002217175.1": "Folsomia candida"}

rule all:
    input:
        # sketches (search genomes and bins) 
        f"{outdir}/search-genomes.sig.zip",
        f"{outdir}/bins.sig.zip",
        # search genomes x metagenomes (branchwater k31, sc1000)
        f"{outdir}/search-genomes-x-brmetagenomes.manysearch.csv",
        # bins x search genomes
        f"{outdir}/bins-x-search-genomes.multisearch.csv",
        # bins x NCBI entire scaled 10,000
        f"{outdir}/bins-x-ncbi-entire.manysearch.csv",
        # bins x NCBI euks scaled 1,000
        #f"{outdir}/bins-x-ncbi-euks.multisearch.sc1000.csv",
        # aggregated results
        f"{outdir}/multi.aggregated-results.csv"


#### 1. Sketch the search genomes ####
rule write_directsketch_csv:
    output:
        csv=f"{outdir}/search-genomes.gbsketch.csv"
    run:
        with open(output.csv, 'w') as f:
            f.write("accession,name\n")
            for acc, species_name in search_genomes.items():
                name = f"{acc} {species_name}"
                f.write(f"{acc},{name}\n")


rule sourmash_sketch_search_genomes:
    # use scaled=100 to allow higher resolution search vs bins (sc100)
    input:
        csv=f"{outdir}/search-genomes.gbsketch.csv",
    output:
        sig=f"{outdir}/search-genomes.sig.zip",
    threads: 4
    benchmark: f"{outdir}/logs/sketch-search-genomes.benchmark"
    shell:
        """
        sourmash scripts gbsketch -p dna,k=21,k=31,k=51,scaled=100 {input.csv} -o {output.sig}
        """

#### 2. Search the branchwater search genomes against the metagenomes, but with k=31 ###
rule sourmash_manysearch_metagenomes:
    """
    search genomes against a set of metagenomes using sourmash manysearch
    params: k 31, scaled 1000
    """
    input:
        genomes_zip=f"{outdir}/search-genomes.sig.zip",
        metagenomes_mf = "br-magtest.wort.mf.csv",
    output:
        csv=f"{outdir}/search-genomes-x-brmetagenomes.manysearch.csv"
    threads: 4
    log: f"{outdir}/logs/search-genomes-x-brmetagenomes.manysearch.log"
    benchmark: f"{outdir}/logs/search-genomes-x-brmetagenomes.manysearch.benchmark"
    shell:
        """
        sourmash scripts manysearch {input.genomes_zip} {input.metagenomes_mf} -k 31 --scaled 1000 --output {output.csv} > {log} 2>&1
        """

#### 3. Sketch the MAG bins ####
rule write_manysketch_bins_csv:
    input:
        bins="multi_bins.txt",
    output:
        csv=f"{outdir}/bins.manysketch.csv"
    run:
        with open(output.csv, 'w') as f:
            f.write("name,genome_filename,protein_filename\n")
            # iterate over input fasta files and write to csv
            bin_list = [i.strip() for i in open(input.bins).readlines()]
            for bin_file in bin_list:
                sample_bin = os.path.basename(bin_file).replace("MEGAHIT-MetaBAT2-", "").replace(".fa", "")
                f.write(f"{sample_bin},{bin_file},\n")

rule manysketch_bins:
    input:
        csv=f"{outdir}/bins.manysketch.csv"
    output:
        sig=f"{outdir}/bins.sig.zip"
    threads: 4
    benchmark: f"{outdir}/logs/sketch-bins.benchmark"
    shell:
        """
        sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=100 {input.csv} -o {output.sig}
        """

#### 4. Use sourmash comparisons to identify/classify bins ####

#### 4a. Compare bins against the search genomes ####
rule multisearch_bins_x_search_genomes:
    """
    search bins against the search genomes using sourmash manysearch
    """
    input:
        genomes_zip=f"{outdir}/search-genomes.sig.zip",
        bins_sig=f"{outdir}/bins.sig.zip",
    output:
        csv=f"{outdir}/bins-x-search-genomes.multisearch.csv"
    threads: 4
    log: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc100.log"
    benchmark: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc100.benchmark"
    shell:
        """
        sourmash scripts multisearch {input.bins_sig} {input.genomes_zip} -k 31 --scaled 100 --output {output.csv} --ani -m DNA --threshold 0.001 > {log} 2>&1
        """

#### 4b. Search the bins against NCBI database (scaled 10,000) ####
rule sourmash_manysearch_ncbi_entire:
    """
    search bins against the entire NCBI database using sourmash manysearch
    params: k 51, scaled 10000
    database: entire-2025-01-21.k51.rocksdb
    """
    input:
        zip=f"{outdir}/bins.sig.zip",
        db=db_ncbi_entire,
    output:
        manysearch=f"{outdir}/bins-x-ncbi-entire.manysearch.csv"
    log: f"{outdir}/logs/bins-x-ncbi-entire.manysearch.log"
    benchmark: f"{outdir}/logs/bins-x-ncbi-entire.manysearch.benchmark"
    shell:
        """
        sourmash scripts manysearch --scaled 10_000 -k 51 {input.zip} {input.db} -o {output.manysearch} > {log} 2>&1
        """

#### 4c. Search the bins against NCBI eukaryote databases (scaled 1000) ####
rule write_euk_db_file:
    """
    Write a simple text file containing the euk databases.
    This allows us to use multiple databases with the branchwater plugin.
    """
    input:
        fungi=db_fungi,
        euk_other=db_euk_other,
        bilateria_minus_verts=db_bilateria_minus_verts,
    output:
        db_txt="euk-dbs.txt"
    run:
        with open(output.db_txt, 'w') as f:
            f.write(f"{input.fungi}\n")
            f.write(f"{input.euk_other}\n")
            f.write(f"{input.bilateria_minus_verts}\n")


rule sourmash_multisearch_ncbi_euk_zips:
    """
    search bins against the NCBI eukaryote databases using sourmash multisearch
    params: k 51, scaled 1000, ani
    databases: fungi, euk-other, bilateria-minus-verts
    """
    input:
        zip=f"{outdir}/bins.sig.zip",
        fungi=db_fungi,
        euk_other=db_euk_other,
        bilateria_minus_verts=db_bilateria_minus_verts,
        db_txt="euk-dbs.txt"
    output:
        multisearch=f"{outdir}/bins-x-ncbi-euks.multisearch.sc1000.csv"
    log: f"{outdir}/logs/bins-x-ncbi-euks.multisearch.sc1000.log"
    benchmark: f"{outdir}/logs/bins-x-ncbi-euk.multisearch.sc1000.benchmark"
    shell:
        """
        sourmash scripts multisearch --scaled 1_000 -m DNA -k 51 --ani {input.zip} {input.db_txt} -o {output.multisearch} > {log} 2>&1
        """


# 5. aggregate and assess the results
rule aggregate_results:
    """
    aggregate the results sourmash, sendsketch, BAT annotation; compare to expected taxonomy
    """
    input:
        expected_taxonomy=expected_taxonomy,
        bins_file = bins_file,
        # search genomes x metagenomes (branchwater k31, sc1000)
        mgx_k31 = f"{outdir}/search-genomes-x-brmetagenomes.manysearch.csv",
        # bins x search genomes
        bins_x_searchgx = f"{outdir}/bins-x-search-genomes.multisearch.csv",
        # bins x NCBI euks scaled 1,000
        #bins_x_ncbi = f"{outdir}/bins-x-ncbi-euks.multisearch.sc1000.csv",
        # sendsketch results (already aggregated with brachwater-web results and metadata)
        sendsketch=sendsketch_results,
        # bins x NCBI entire scaled 10,000
        bins_x_ncbi = f"{outdir}/bins-x-ncbi-entire.manysearch.csv",
    output:
        summary_csv=f"{outdir}/multi.aggregated-results.csv",
        bin_summary_csv=f"{outdir}/multi.aggregated-results-bins.csv",
    log: f"{outdir}/logs/aggregate-results.log"
    benchmark: f"{outdir}/logs/aggregate-results.benchmark"
    shell:
        """
        python3 scripts/aggregate-results.py --expected-csv {input.expected_taxonomy} \
                                             --bins {input.bins_file} \
                                             --mgx-manysearch-csv {input.mgx_k31} \
                                             --bin-manysearch-ncbi-csv {input.bins_x_ncbi} \
                                             --multisearch-bins {input.bins_x_searchgx} \
                                             --sendsketch-csv {input.sendsketch} \
                                             --bin-summary {output.bin_summary_csv} \
                                             --output {output.summary_csv} > {log} 2>&1
        """
