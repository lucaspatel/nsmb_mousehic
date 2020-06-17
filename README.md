# Haplotype Specific Hi-C Pipeline

A pipeline to perform automated mapping, haplotype-assignment, and Hi-C data generation, as in [Dynamic reorganization of the genome shapes the recombination landscape in meiotic prophase](https://www.nature.com/articles/s41594-019-0187-0).

The pipeline ingests paired-end .fastq files from a Hi-C experiment and runs a complete Hi-C analysis pipeline, specifically for the haplotype analysis. 

## Getting Started

This pipeline is optimized for PBS/TORQUE systems. To get started, you'll need to set up a few things:
* Clone the repository to your working directory.
* Add your index files to the `/index` directory
* Add your SNP files to the `/snp` directory
    * The required format is one file per chromosome named `chr#.snps.txt.gz` with each file containing SNPs in the tab-seperated format `POS REF ALT`
* Add your sequence files to the `/fastq` directory
    * By default, the pipeline expects bz2-compressed .fastq files
    * The compression format can be configured in `config.json`
* Configure the `config.json` file with the paths to the required software packages, .fastq prefixes, assemblies, and chromosome references for the pipeline
* Configure the jobscript.pbs with your respective email and username used in your job submission queue. Also load any modules that may be needed here.

### Prerequisites

This pipeline leverages several bioinformatics tools, namely:
* Samtools
* BWA (MEM)
* [WASP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4626402/)
* Juicer
* Cooler
* Pairtools

## Running

To run the pipeline, simply execute `run.sh`. I recommend backgrounding the process, or submitting it to a job scheduler.

```
./run.sh &
```

## FAQ

> What packages are necessary to run the pipeline?

The main packages are listed above, but a full conda environment can be found in `environment.yml`.

> The pipeline failed or completed before generating the final output files, what happened?

Consult the log files in `/mapping/log` and `/pbslog` for more information.

> I can't restart the pipeline after it failed!

Snakemake locks the working directory while running. You may have to unlock the working directory when restarting the pipeline with `snakemake --unlock`.

## Built With

* [Python](https://www.python.org)
    * The pipeline must be executed in a Python 3 environment, but Python 2.7 is still required for a few legacy scripts.
* [Snakemake](http://www.dropwizard.io/1.0.2/docs/) - A pipeline tool for Python

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Lucas Patel** - *Initial and ongoingn work* 
* **Yunjiang Qiu** - *Initial pipeline structure and Hi-C WASP adaptation* 
* **David Gorkin** - *Hi-C WASP adaptation* 
* **Kevin Corbett** - *Experimental design and scientific direction* 

## Citation

This pipeline was created to facilitate the experimental design of the corresponding publicatio, [Dynamic reorganization of the genome shapes the recombination landscape in meiotic prophase](https://www.nature.com/articles/s41594-019-0187-0). If this pipeline is used in your own analyis, please cite the above approrpriately.


## Acknowledgments

* Thanks to Grant McVicker for early guidance.
* Thanks to Anton Goloborodko for pairing score guidance and code.

* The `homolog_pairing.py` script is adapted from `_____`
