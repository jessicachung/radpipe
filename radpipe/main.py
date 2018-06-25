'''Run the pipeline as a program.'''

from radpipe.pipeline import make_pipeline
import radpipe.pipeline_base.main
import pkg_resources

PROGRAM_NAME = "radpipe"
PROGRAM_INFO = pkg_resources.require(PROGRAM_NAME)[0]
PROGRAM_VERSION = PROGRAM_INFO.version

def main():
    radpipe.pipeline_base.main.main(PROGRAM_NAME, PROGRAM_VERSION, make_pipeline)
