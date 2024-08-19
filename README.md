# NCBI Genomes downloader

Download, extract, and flatten NCBI genome archive.


## Why

The NCBI genome archive allows you to download reference genomes through a web interface or through
REST API in a nice `zip` archive.
There is only one issue, it comes in a deep nested structure with a lot of ancilary files.

For instance, after downloading the reference genome of Rice (_Oryza sativa Japonica Group_) with ID `GCF_034140825.1` through the web interface, selecting only the `gtf` file to minimize the download size, the `ncbi_dataset.zip` has the following structure:

```
ncbi_dataset.zip
├── ncbi_dataset
│   └── data
│       ├── assembly_data_report.jsonl
│       ├── dataset_catalog.json
│       ├── data_summary.tsv
│       └── GCF_034140825.1
│           └── genomic.gtf
└── README.md
```

That is a lot of files and folders despite requesting only the `gtf` file.
Typically, only the genomic `fasta` and `gtf` files are required for most analyses.
On top of that, the `gtf` file doesn't even have a nice file name,
having rather generic `genomic.gtf`.
This makes it quite laborious to integrate this easily into pipelines
since a lot of manual intervention: unpacking, sorting, perhaps even renaming, is required.

Here comes `ncbi.r`, a stand-alone R script with no dependencies that does this work for you.
Only up-to-date R version is required.

## Installation

Just copy paste the script where you need it and run it directly with `Rscript`.

Alternatively, add it to your `PATH` and mark it as executable.
See your OS help pages on how to do that.

For instance, on my Linux machine, my home directory contains a `bin` folder that was added to `PATH` system variable. I can just copy the `ncbi.r` there, type:

```
chmod +x $HOME/bin/ncbi.r
```

and now I can simply type `ncbi.r` to run the script from anywhere.

## Usage

If `ncbi.r` is not in your path, navigate to where you downloaded the file and type:

```
Rscript ncbi.r --help
```

This will display usage instructions.

To download, extract, and flatten the `GCF_034140825.1` genome, type:

```
Rscript ncbi.r GCF_034140825.1
```

This will download the Rice reference genome and create two files in your current path:
`GCF_034140825.1.fna` and `GCF_034140825.1.gtf`.
