# radpipe

A bioinformatics pipeline for double digest RAD-seq analysis with support for
running pipeline stages on a distributed compute cluster.

radpipe is based on the Ruffus pipeline library and predominantly uses the
Stacks software suite for the RAD-seq processing steps.

# Dependencies

- [Python 3](https://www.python.org/downloads/)
- [Ruffus](http://www.ruffus.org.uk/index.html)
- [DRMAA](http://apps.man.poznan.pl/trac/slurm-drmaa)
- [Stacks >= 2.0](http://catchenlab.life.illinois.edu/stacks/)
- [Samtools](http://www.htslib.org/download/)
- [BWA](http://bio-bwa.sourceforge.net/) or
  [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](http://multiqc.info/)

# Installation

It is recommended you install `radpipe` in a virtual environment. For example,

```
cd /place/to/install
virtualenv -p /usr/bin/python3 radpipe_venv
source radpipe_venv/bin/activate
pip install -U git+https://github.com/jessicachung/ruffus
pip install -U git+https://github.com/pearg/radpipe
```

# Usage

### Set DRMAA path

If you're using the pipeline on a cluster with a job scheduling system, you'll
need to tell the cluster where your DRMAA library is with `DRMAA_LIBRARY_PATH`.

```
DRMAA_LIBRARY_PATH=/mnt/galaxy/gvl/software/slurm-drmaa-1.0.7/slurm_drmaa/.libs/libdrmaa.so
export DRMAA_LIBRARY_PATH
```

### Configuration file

Running the pipeline requires a configuration file. You can find an example
config file in [example/pipeline.config](example/pipeline.config).

### Dry run

Print out steps the pipeline will run with `-n`. Specify the target task the
pipeline will complete with `--target_tasks`. Valid target tasks are listed
in the Stages section below.

```
radpipe \
    --config pipeline.config \
    --verbose 3 \
    --target_tasks multiqc_fastqc \
    -n
```

### Running the pipeline

Run the pipeline, allowing multiple jobs to be submitted to the queue with
`--jobs` and `--use_threads`.

```
radpipe \
    --config pipeline.config \
    --verbose 3 \
    --log_file radpipe.log \
    --target_tasks multiqc_fastqc \
    --jobs 4 \
    --use_threads
```

### Help message

```
$ radpipe --help
usage: radpipe [-h] [--verbose [VERBOSE]] [-L FILE] [-T JOBNAME] [-j N]
               [--use_threads] [-n] [--touch_files_only] [--recreate_database]
               [--checksum_file_name FILE] [--flowchart FILE]
               [--key_legend_in_graph] [--draw_graph_horizontally]
               [--flowchart_format FORMAT] [--forced_tasks JOBNAME]
               [--config CONFIG] [--jobscripts JOBSCRIPTS] [--version]

RAD-Seq pipeline

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG       Pipeline configuration file in YAML format, defaults
                        to pipeline.config
  --jobscripts JOBSCRIPTS
                        Directory to store cluster job scripts created by the
                        pipeline, defaults to jobscripts
  --version             show program's version number and exit

Common options:
  --verbose [VERBOSE], -v [VERBOSE]
                        Print more verbose messages for each additional
                        verbose level.
  -L FILE, --log_file FILE
                        Name and path of log file

pipeline arguments:
  -T JOBNAME, --target_tasks JOBNAME
                        Target task(s) of pipeline.
  -j N, --jobs N        Allow N jobs (commands) to run simultaneously.
  --use_threads         Use multiple threads rather than processes. Needs
                        --jobs N with N > 1
  -n, --just_print      Don't actually run any commands; just print the
                        pipeline.
  --touch_files_only    Don't actually run any commands; just 'touch' the
                        output for each task to make them appear up to date.
  --recreate_database   Don't actually run any commands; just recreate the
                        checksum database.
  --checksum_file_name FILE
                        Path of the checksum file.
  --flowchart FILE      Don't run any commands; just print pipeline as a
                        flowchart.
  --key_legend_in_graph
                        Print out legend and key for dependency graph.
  --draw_graph_horizontally
                        Draw horizontal dependency graph.
  --flowchart_format FORMAT
                        format of dependency graph file. Can be 'svg', 'svgz',
                        'png', 'jpg', 'psd', 'tif', 'eps', 'pdf', or 'dot'.
                        Defaults to the file name extension of --flowchart
                        FILE.
  --forced_tasks JOBNAME
                        Task(s) which will be included even if they are up to
                        date.
```

### How to kill a running pipeline

If you want to kill a running pipeline that is running in cluster mode, and
`^C` isn't killing the pipeline process, you will need to cancel all individual
pipeline jobs in the system. On SLURM, you can cancel your jobs with
`scancel <job_id>`. If all your jobs are from radpipe, you can cancel all jobs
submitted with `scancel --user=<your_username>`.


# Stages

- **`fastqc`**: Runs FastQC on raw sequencing files and outputs to `results/qc/fastqc/`.
- **`multiqc_fastqc`**: Runs MultiQC on FastQC outputs and outputs to `results/qc/`.
- **`process_radtags`**: Runs Stacks process_radtags to demux samples and outputs to `results/sample_radtags/`.
- **`build_index`**: Builds index for the reference FASTA file for either bwa or bowtie and stores indices in `results/ref/`.
- **`alignment`**: Aligns FASTQ files with either bwa or bowtie and outputs BAMs to `results/alignments/`.
- **`sort_bam`**: Sorts BAM files by coordinate and generates indices. Outputs to `results/alignments/`.
- **`filter_bam`**: Optional step to filter BAM files using Samtools view. Outputs to `results/alignments/`.
- **`flagstat`**: Runs Samtools flagstat on the final BAM files and outputs to `results/qc/flagstat/`.
- **`multiqc_flagstat`**: Runs MultiQC on the flagstat outputs and outputs to `results/qc/`.
- **`gstacks`**: Runs Stacks gstacks on sorted BAMs and outputs to `results/gstacks/`.
- **`populations`**: Runs Stacks populations and outputs to `results/populations/`.

You can also generate a flowchart image with the `--flowchart` option.

```bash
radpipe \
    --config pipeline.config \
    --target_tasks populations,multiqc_fastqc,multiqc_flagstat \
    --flowchart_format png \
    --flowchart radpipe.png
```
