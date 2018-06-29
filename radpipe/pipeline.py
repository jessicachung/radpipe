'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from radpipe.stages import PipelineStages
from radpipe.pipeline_base.utils import safe_make_dir
from radpipe.libraries import parse_libraries
from radpipe.utils import path_list_join, get_output_paths
from collections import Counter, OrderedDict
import radpipe.error_codes
import logging
import sys
import os

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.DEBUG)

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name="radpipe")

    # Stages are dependent on the state
    stages = PipelineStages(state)

    # Get a list of library objects.
    libraries = parse_libraries(libraries=state.config.get_options("libraries"))

    # Get a list of input files
    input_files = [l.files for l in libraries]
    # input_files = [item for sublist in input_files for item in sublist]
    logging.info("Input files: " + str(input_files))

    # Get a list of all samples for each library
    samples_dict = OrderedDict()
    for l in libraries:
        samples_dict[l.name] = l.samples
    logging.debug("Samples: " + str(samples_dict))

    # Make sure that there are no duplicate samples
    sample_list = [item for sublist in samples_dict.values() for item in sublist]
    sample_counts = Counter(sample_list)
    for sample in sample_counts:
        if sample_counts[sample] > 1:
            logging.error("Sample {} appears {} times in the barcodes files. "
                          "Sample names must be unique".format(sample,
                              sample_counts[sample]))
            sys.exit(radpipe.error_codes.INVALID_INPUT_FILE)

    # Define output directories
    output_dir = get_output_paths(state)
    logging.debug(output_dir)

    # Check if alignment_method is valid
    alignment_method = state.config.get_options("alignment_method").strip().lower()
    if alignment_method not in ["bwa mem", "bowtie"]:
        logging.error("Error: Invalid alignment_method in config file. " \
                      "Valid options are ['bwa mem', 'bowtie'].")
        sys.exit(radpipe.error_codes.INVALID_ARGUMENT)

    # If 'alignment' is in target_tasks, specify which type of alignment job
    # TODO: allow multiple comma-separated tasks
    if "alignment" in state.options.target_tasks:
        if alignment_method == "bwa mem":
            state.options.target_tasks = ["bwa_align"]
        elif alignment_method == "bowtie":
            state.options.target_tasks = ["bowtie_align"]
    logging.debug(state)

    # Population map filenames
    popmap_file = "{output_dir}/{name}_popmap.txt".format(
        output_dir=output_dir["populations"],
        name=state.config.get_options("analysis_id"))
    try:
        config_popmap_file = state.config.get_options("popmap_file")
        if config_popmap_file:
            logging.info("Using popmap file: {}".format(config_popmap_file))
        else:
            raise(Exception)
    except Exception:
        config_popmap_file = None
        logging.info("Creating new popmap file: {}".format(popmap_file))

    # Population r values
    populations_r = state.config.get_options("populations_r")
    assert(isinstance(populations_r, list))

    # Dummy stages. These do nothing except provide a node at the beginning
    # for the pipeline graph, giving the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.do_nothing,
        name="original_fastqs",
        output=input_files)

    pipeline.originate(
        task_func=stages.do_nothing,
        name="reference_genome",
        output=state.config.get_options("reference_genome"))

    # Create a copy of the population map file needed for stacks, or create
    # one denovo using the sample list.
    pipeline.originate(
        task_func=stages.create_popmap_file,
        name="create_popmap_file",
        output=[popmap_file],
        extras=[config_popmap_file, sample_list])

    # Create index for reference genome based on alignment method.
    if alignment_method == "bwa mem":
        align_task_name = "bwa_mem"
        pipeline.transform(
            task_func=stages.bwa_index,
            name="bwa_index",
            input=output_from("reference_genome"),
            filter=formatter(".+/(?P<ref>[^/]+).(fa|fasta)"),
            output=path_list_join(output_dir["reference"],
                       ["reference.fa.bwt", "reference.fa.sa"]),
            extras=[output_dir["reference"]])

    if alignment_method == "bowtie":
        align_task_name = "bowtie"
        pipeline.transform(
            task_func=stages.bowtie_index,
            name="bowtie_index",
            input=output_from("reference_genome"),
            filter=formatter(".+/(?P<ref>[^/]+).(fa|fasta)"),
            output=path_list_join(output_dir["reference"],
                       ["reference.1.ebwt", "reference.rev.1.ebwt"]),
            extras=[output_dir["reference"]])

    # FastQC
    pipeline.transform(
        task_func=stages.fastqc,
        name="fastqc",
        input=output_from("original_fastqs"),
        filter=formatter(".+/(?P<lib>[^/]+)/(?P<fn>[^/]+).(fastq|fq).gz"),
        output="%s/{lib[0]}/{fn[0]}_fastqc.zip" % output_dir["fastqc"],
        extras=[output_dir["fastqc"], "{lib[0]}"])

    # MultiQC
    pipeline.merge(
        task_func=stages.multiqc,
        name="multiqc",
        input=output_from("fastqc"),
        output="%s/multiqc_report.html" % output_dir["fastqc"],
        extras=[output_dir["fastqc"]])

    # Stacks: Process RAD-Tags
    pipeline.transform(
        task_func=stages.process_radtags,
        name="process_radtags",
        input=output_from("original_fastqs"),
        filter=formatter(".+/(?P<lib>[^/]+)/[^/]+"),
        output="%s/{lib[0]}/{lib[0]}.success" % output_dir["process_radtags"],
        extras=[output_dir["process_radtags"], "{lib[0]}",
                state.config.get_options("renz_1"),
                state.config.get_options("renz_2"),
                state.config.get_options("process_radtags_options")]
    ).follows("fastqc")

    # Create a list for alignment with the input fastq files from process_radtags
    process_radtags_outputs = []
    for l in libraries:
        for s in l.samples:
            base = "{dir}/{lib}/{sample}".format(
                dir=output_dir["process_radtags"], lib=l.lib_id, sample=s)
            process_radtags_outputs.append([base + ".1.fq.gz",
                                            base + ".2.fq.gz"])
    # print(process_radtags_outputs)

    # Alignment
    if align_task_name == "bwa_mem":
        (pipeline.transform(
            task_func=stages.bwa_align,
            name=align_task_name,
            input=process_radtags_outputs,
            filter=formatter(".+/(?P<sm>[^/]+).1.fq.gz"),
            output="%s/{sm[0]}.bwa.bam" % output_dir["alignments"],
            extras=[os.path.join(output_dir["reference"], "reference.fa"),
                    "{path[0]}", output_dir["alignments"], "{sm[0]}",
                    state.config.get_options("alignment_options")])
        ).follows("bwa_index")

    if align_task_name == "bowtie":
        (pipeline.transform(
            task_func=stages.bowtie_align,
            name=align_task_name,
            input=process_radtags_outputs,
            filter=formatter(".+/(?P<sm>[^/]+).1.fq.gz"),
            output="%s/{sm[0]}.bowtie.bam" % output_dir["alignments"],
            extras=[os.path.join(output_dir["reference"], "reference"),
                    "{path[0]}", output_dir["alignments"], "{sm[0]}",
                    state.config.get_options("alignment_options")])
        ).follows("bowtie_index")

    # Sort BAM and index
    pipeline.transform(
        task_func=stages.sort_bam,
        name="sort_bam",
        input=output_from(align_task_name),
        filter=suffix(".bam"),
        output=".sorted.bam")

    # TODO: filter with samtools view

    # Stacks: gstacks
    pipeline.merge(
        task_func=stages.gstacks,
        name="gstacks",
        input=output_from(stages.sort_bam),
        output="%s/catalog.fa.gz" % output_dir["gstacks"],
        extras=[output_dir["alignments"],
                output_dir["gstacks"],
                align_task_name,
                sample_list])

    # Define outputs from each run of populations
    populations_outputs = []
    for r in populations_r:
        dir_name = "{pop_dir}/{analysis_name}_r{r}".format(
            pop_dir=output_dir["populations"],
            analysis_name=state.config.get_options("analysis_id"),
            r=r)
        populations_outputs.append(os.path.join(dir_name, "populations.snps.vcf"))
    # print(populations_outputs)

    # Stacks: populations
    pipeline.originate(
        task_func=stages.populations,
        name="popluations",
        output=populations_outputs,
        extras=[output_dir["gstacks"], output_dir["populations"], popmap_file,
                state.config.get_options("populations_options")]
    ).follows("gstacks").follows("create_popmap_file")

    return pipeline
