'''
Various utility functions that don't have a sensible home elsewhere
'''
import os.path

def get_output_paths(state):
    results_dir = state.config.get_options("results_dir")
    output_path = {
        "reference": "reference",
        "fastqc": "fastqc",
        "process_radtags": "sample_radtags",
        "alignments": "alignments",
        "gstacks": "gstacks",
        "populations": "populations"
    }
    output_path = [(a, os.path.join(results_dir, b)) for a,b in output_path.items()]
    output_path = dict(output_path)
    return output_path

def path_list_join(dir, file_list):
    '''Join directory to a list of files'''
    return [os.path.join(dir, x) for x in file_list]

def run_picard(state, stage, args):
    mem = int(state.config.get_stage_options(stage, "mem"))
    return run_java(state, stage, PICARD_JAR, mem, args)

def run_trimmomatic(state, stage, args):
    mem = int(state.config.get_stage_options(stage, "mem"))
    return run_java(state, stage, TRIMMOMATIC_JAR, mem, args)

def create_empty_outputs(outputs):
    '''Create empty dummy files for testing purposes'''
    if isinstance(outputs, str):
        outputs = [outputs]
    for output_filename in outputs:
        with open(output_filename, "w"):
            pass
