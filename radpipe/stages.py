'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from radpipe.pipeline_base.utils import safe_make_dir, run_java
from radpipe.pipeline_base.stages import Stages
from radpipe.pipeline_base.runner import run_stage
import os
from math import floor
import shutil
from radpipe.utils import *


# java_tmp = "-Djava.io.tmpdir=$TMPDIR"

class PipelineStages(Stages):
    def __init__(self, *args, **kwargs):
        super(PipelineStages, self).__init__(*args, **kwargs)
        self.reference_genome = self.get_options("reference_genome")

    def do_nothing(self, *args):
        '''Do nothing'''
        pass

    def create_popmap_file(self, output, config_popmap_file, sample_list):
        '''
        Copy population map for populations if provided, else create a
        population map file using the sample list.
        '''
        safe_make_dir(os.path.dirname(output))
        if config_popmap_file:
            # Copy popmap file if user provided
            shutil.copyfile(config_popmap_file, output)
        else:
            # If no provided popmap file, create popmap using sample_list
            popmap = ["{}\t1".format(x) for x in sample_list]
            with open(output, "w") as f:
                f.write("\n".join(popmap))

    def bwa_index(self, input, outputs, reference_dir):
        '''Index reference genome for bwa alignment'''
        safe_make_dir(reference_dir)
        command = "ln -sf {ref_fasta} {ref_symlink} && bwa index {ref_symlink}".format(
                      ref_fasta=os.path.abspath(input),
                      ref_symlink=os.path.join(reference_dir, "reference.fa"))
        run_stage(self.state, "build_index", command)

    def bowtie_index(self, input, outputs, reference_dir):
        '''Index reference genome for bowtie alignment'''
        safe_make_dir(reference_dir)
        command = "ln -sf {ref_fasta} {ref_symlink} && bowtie-build " \
                  "{ref_symlink} {index_base}".format(
                      ref_fasta=os.path.abspath(input),
                      ref_symlink=os.path.join(reference_dir, "reference.fa"),
                      index_base= os.path.join(reference_dir, "reference"))
        run_stage(self.state, "build_index", command)

    def fastqc(self, input, outputs, fastqc_dir, lib):
        '''Run FastQC on fastq files'''
        safe_make_dir(fastqc_dir)
        fastqc_output_dir = os.path.join(fastqc_dir, lib)
        safe_make_dir(fastqc_output_dir)
        assert(isinstance(input, list))
        # Remove barcode file
        fastq_input = input[:-1]
        fastq_input = " ".join(fastq_input)
        command = "fastqc -o {fastqc_output_dir} -f fastq {fastq_input}".format(
                      fastqc_output_dir=fastqc_output_dir,
                      fastq_input=fastq_input)
        run_stage(self.state, "fastqc", command)

    def multiqc_fastqc(self, input, output, qc_dir, fastqc_dir):
        '''Run MultiQC on the FastQC directory'''
        command = "multiqc --module fastqc --outdir {qc_dir} " \
                  "--filename multiqc_fastqc {fastqc_dir}".format(
                      qc_dir=qc_dir, fastqc_dir=fastqc_dir)
        run_stage(self.state, "multiqc", command)

    def process_radtags(self, inputs, output, output_dir, lib, re_1, re_2,
                        extra_options):
        '''Process radtags to separate into separate fastq files'''
        lib_output_dir = os.path.join(output_dir, lib)
        safe_make_dir(lib_output_dir)
        command = "process_radtags -1 {r1} -2 {r2} -b {barcodes_file} " \
                  "-i gzfastq -o {lib_output_dir} --inline_inline --renz_1 {re_1} " \
                  "--renz_2 {re_2} {extra_options} && touch {success_file}".format(
                      r1=inputs[0], r2=inputs[1], barcodes_file=inputs[2],
                      lib_output_dir=lib_output_dir, re_1=re_1, re_2=re_2,
                      extra_options=extra_options, success_file=output)
        run_stage(self.state, "process_radtags", command)

    def bwa_align(self, inputs, output, ref_fasta, input_path, output_path, sm,
                  extra_options):
        '''Align fastq files with BWA mem'''
        safe_make_dir(os.path.dirname(output))
        # cores = min(self.get_stage_options("alignment", "cores") - 1, 1)
        cores = self.get_stage_options("alignment", "cores")
        # If PE fastq inputs, join into a string
        if isinstance(inputs, list):
            fastq_input = " ".join(inputs)
        else:
            fastq_input = inputs
        lb = os.path.basename(input_path)
        log_file = os.path.join(output_path, "{sm}.bwa.log".format(sm=sm))
        read_group_string = "@RG\tID:{lb}\tLB:{lb}\tPU:{lb}\tPL:ILLUMINA\tSM:{sm}" \
                            "".format(lb=lb, sm=sm)
        command = 'set -o pipefail; bwa mem -t {cores} {extra_options} ' \
                  '-R "{rg}" {ref_fasta} {fastq_input} 2> {log_file} | ' \
                  'samtools view -b - > {output}' .format(
                          cores=cores, extra_options=extra_options,
                          rg=read_group_string, ref_fasta=ref_fasta,
                          fastq_input=fastq_input, output=output,
                          log_file=log_file)
        # Use bash instead of sh
        command = "bash -c '{}'".format(command)
        run_stage(self.state, "alignment", command)

    def bowtie_align(self, inputs, output, index_base, input_path, output_path,
                     sm, extra_options):
        '''Align fastq files with Bowtie'''
        safe_make_dir(os.path.dirname(output))
        cores = self.get_stage_options("alignment", "cores")
        # If PE fastq inputs, join into a string
        # Bowtie doesn't accept gzipped files
        if isinstance(inputs, list):
            assert(len(inputs) == 2)
            fastq_input = "-1 <(zcat {r1}) -2 <(zcat {r2})" \
                          .format(r1=inputs[0], r2=inputs[1])
        else:
            fastq_input = "<(zcat {})".format(inputs)
        lb = os.path.basename(input_path)
        log_file = os.path.join(output_path, "{sm}.bowtie.log".format(sm=sm))
        command = "set -o pipefail; bowtie --threads {cores} --sam {extra_options} " \
                  "--sam-RG ID:{lb} --sam-RG LB:{lb} --sam-RG PU:{lb} " \
                  "--sam-RG PL:ILLUMINA --sam-RG SM:{sm} " \
                  "{index_base} {fastq} 2> {log_file} | samtools view -b - " \
                  "> {output}".format(
                          cores=cores, extra_options=extra_options, lb=lb,
                          sm=sm, index_base=index_base, fastq=fastq_input,
                          log_file=log_file, output=output)
        # Use bash instead of sh
        command = 'bash -c "{}"'.format(command)
        run_stage(self.state, "alignment", command)

    def sort_bam(self, input, output):
        '''Sort BAM file by coordinates'''
        cores = self.get_stage_options("sort_bam", "cores")
        mem = max(floor(self.get_stage_options("sort_bam", "mem") / cores) - 1, 1)
        command = "samtools sort -@ {cores} -m {mem}G -o {output} {input} " \
                  "&& samtools index {output}".format(cores=cores, mem=mem,
                          output=output, input=input)
        run_stage(self.state, "sort_bam", command)

    def filter_bam(self, input, output, extra_options):
        '''Filter BAM file with Samtools view'''
        command = "samtools view -b {extra_options} {input} > {output} && " \
                  "samtools index {output}".format(extra_options=extra_options,
                          input=input, output=output)
        run_stage(self.state, "filter_bam", command)

    def flagstat(self, input, output):
        '''Run Samtools flagstat on final BAM files'''
        safe_make_dir(os.path.dirname(output))
        command = "samtools flagstat {input} > {output}".format(
                          input=input, output=output)
        run_stage(self.state, "flagstat", command)

    def multiqc_flagstat(self, input, output, qc_dir, flagstat_dir):
        '''Run MultiQC on the flagstat directory'''
        command = "multiqc --module samtools --outdir {qc_dir} " \
                  "--filename multiqc_flagstat {flagstat_dir}".format(
                      qc_dir=qc_dir, flagstat_dir=flagstat_dir)
        run_stage(self.state, "multiqc", command)

    def gstacks(self, inputs, output, input_dir, output_dir, aligner_name,
                final_bam_name, sample_list):
        '''Run gstacks'''
        safe_make_dir(output_dir)
        # Create popmap file using sample_list
        popmap_filename = os.path.join(output_dir, "popmap.txt")
        gstacks_popmap = ["{}\t1".format(x) for x in sample_list]
        with open(popmap_filename, "w") as f:
            f.write("\n".join(gstacks_popmap))
        cores = self.get_stage_options("gstacks", "cores")
        if aligner_name == "bwa_mem":
            suffix = ".bwa.sorted"
        else:
            suffix = ".{aligner}.sorted".format(aligner=aligner_name)
        if final_bam_name == "filter_bam":
            suffix = suffix + ".filtered.bam"
        else:
            suffix = suffix + ".bam"
        command = "gstacks -t {cores} -M {popmap} -I {input_dir} -S {suffix} " \
                  "-O {output_dir}".format(cores=cores, popmap=popmap_filename,
                          input_dir=input_dir, suffix=suffix,
                          output_dir=output_dir)
        run_stage(self.state, "gstacks", command)

    def populations(self, output, gstacks_dir, populations_dir, popmap_file,
                    populations_options):
        '''Run Stacks populations'''
        output_dir = os.path.dirname(output)
        safe_make_dir(populations_dir)
        safe_make_dir(output_dir)
        cores = self.get_stage_options("populations", "cores")
        # Copy popmap file to directory
        popmap_copy = os.path.join(output_dir, "popmap.txt")
        shutil.copyfile(popmap_file, popmap_copy)
        # Get r value
        r = output_dir.split("_")[-1][1:]
        command = "populations -P {gstacks_dir} -O {output_dir} -t {cores} " \
                  "-M {popmap} -r {r} --vcf".format(gstacks_dir=gstacks_dir,
                          output_dir=output_dir, cores=cores,
                          popmap=popmap_file, r=r)
        run_stage(self.state, "populations", command)
