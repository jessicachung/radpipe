'''
Various utility functions that don't have a sensible home elsewhere
'''

import os
from radpipe.pipeline_base.runner import run_stage

def safe_make_dir(path):
    '''Make a directory if it does not already exist'''
    if not os.path.exists(path):
        os.makedirs(path)


def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)


def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)
