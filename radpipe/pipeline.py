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
    state.logger.info("Input files: " + str(input_files))

    # Get a list of all samples for each library
    samples_dict = OrderedDict()
    for l in libraries:
        samples_dict[l.name] = l.samples
    state.logger.debug("Samples: " + str(samples_dict))

    # Make sure that there are no duplicate samples
    sample_list = [item for sublist in samples_dict.values() for item in sublist]
    sample_counts = Counter(sample_list)
    for sample in sample_counts:
        if sample_counts[sample] > 1:
            print("Sample {} appears {} times in the barcodes files. "
                  "Sample names must be unique".format(sample,
                      sample_counts[sample]))
            sys.exit(radpipe.error_codes.INVALID_INPUT_FILE)

    # Define output directories
    output_dir = get_output_paths(state)
    state.logger.debug(output_dir)

    # Allow multiple comma-separated tasks
    if len(state.options.target_tasks) == 1:
        state.options.target_tasks = state.options.target_tasks[0].split(",")
    if len(state.options.forced_tasks) == 1:
        state.options.forced_tasks = state.options.forced_tasks[0].split(",")
    state.logger.debug("Target tasks: " + str(state.options.target_tasks))
    state.logger.debug("Forced tasks: " + str(state.options.forced_tasks))

    # Check if alignment_method is valid
    alignment_method = state.config.get_options("alignment_method").strip().lower()
    if alignment_method not in ["bwa mem", "bowtie"]:
        print("Error: Invalid alignment_method in config file. " \
              "Valid options are ['bwa mem', 'bowtie'].")
        sys.exit(radpipe.error_codes.INVALID_ARGUMENT)
    if alignment_method == "bwa mem":
        align_task_name = "bwa_mem"
        index_task_name = "bwa_index"
    else:
        align_task_name = "bowtie"
        index_task_name = "bowtie_index

    # TODO: Refactor this
    # If 'alignment' is in target_tasks or forced_tasks, specify which
    # type of alignment job
    if "alignment" in state.options.target_tasks:
        index = state.options.target_tasks.index("alignment")
        state.options.target_tasks[index] = align_task_name
    if "alignment" in state.options.forced_tasks:
        index = state.options.forced_tasks.index("alignment")
        state.options.forced_tasks[index] = align_task_name

    # If 'build_index' is in target_tasks or forced_tasks, specify which
    # type of index job
    if "build_index" in state.options.target_tasks:
        index = state.options.target_tasks.index("build_index")
        state.options.target_tasks[index] = index_task_name
    if "build_index" in state.options.forced_tasks:
        index = state.options.forced_tasks.index("build_index")
        state.options.forced_tasks[index] = index_task_name
    state.logger.debug(state)

    # Whether to include filter_bam stage or not
    filter_bams = False
    try:
        samtools_view_options = state.config.get_options("samtools_view_options")
        if samtools_view_options:
            filter_bams = True
    except:
        pass
    state.logger.info("Filter bams: {}".format(filter_bams))

    # Population map filenames
    popmap_file = "{output_dir}/{name}_popmap.txt".format(
        output_dir=output_dir["populations"],
        name=state.config.get_options("analysis_id"))
    try:
        config_popmap_file = state.config.get_options("popmap_file")
        if config_popmap_file:
            state.logger.info("Using popmap file: {}".format(config_popmap_file))
        else:
            raise(Exception)
    except Exception:
        config_popmap_file = None
        state.logger.info("Creating new popmap file: {}".format(popmap_file))

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
        pipeline.transform(
            task_func=stages.bwa_index,
            name="bwa_index",
            input=output_from("reference_genome"),
            filter=formatter(".+/(?P<ref>[^/]+).(fa|fasta)"),
            output=path_list_join(output_dir["reference"],
                       ["reference.fa.bwt", "reference.fa.sa"]),
            extras=[output_dir["reference"]])

    if alignment_method == "bowtie":
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

    # MultiQC: FastQC
    pipeline.merge(
        task_func=stages.multiqc_fastqc,
        name="multiqc_fastqc",
        input=output_from("fastqc"),
        output="%s/multiqc_fastqc_report.html" % output_dir["qc"],
        extras=[output_dir["qc"], output_dir["fastqc"]])

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
                state.config.get_options("process_radtags_options")])

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
        ).follows("bwa_index").follows("process_radtags")

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
        ).follows("bowtie_index").follows("process_radtags")

    # Sort BAM and index
    pipeline.transform(
        task_func=stages.sort_bam,
        name="sort_bam",
        input=output_from(align_task_name),
        filter=suffix(".bam"),
        output=".sorted.bam")

    if filter_bams:
        final_bam_task_name = "filter_bam"
        pipeline.transform(
            task_func=stages.filter_bam,
            name="filter_bam",
            input=output_from("sort_bam"),
            filter=suffix(".sorted.bam"),
            output=".sorted.filtered.bam",
            extras=[state.config.get_options("samtools_view_options")])
    else:
        final_bam_task_name = "sort_bam"

    # Samtools flagstat
    pipeline.transform(
        task_func=stages.flagstat,
        name="flagstat",
        input=output_from(final_bam_task_name),
        filter=suffix(".bam"),
        output=".flagstat.txt",
        output_dir=output_dir["flagstat"])

    # MultiQC: flagstat
    pipeline.merge(
        task_func=stages.multiqc_flagstat,
        name="multiqc_flagstat",
        input=output_from("flagstat"),
        output="%s/multiqc_flagstat_report.html" % output_dir["qc"],
        extras=[output_dir["qc"], output_dir["flagstat"]])

    # Stacks: gstacks
    pipeline.merge(
        task_func=stages.gstacks,
        name="gstacks",
        input=output_from(final_bam_task_name),
        output="%s/catalog.fa.gz" % output_dir["gstacks"],
        extras=[output_dir["alignments"],
                output_dir["gstacks"],
                align_task_name,
                final_bam_task_name,
                sample_list,
                state.config.get_options("gstacks_options")])

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
