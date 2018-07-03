'''
A Bioinformatics pipeline based on Ruffus.

Author: Bernie Pope (bjpope@unimelb.edu.au).

Copyright: 2015

See ../README.md for details about how to use the program.

Repository: https://github.com/bjpop/complexo_pipeline
'''

from __future__ import print_function
from ruffus import *
import ruffus.cmdline as cmdline
import sys
from radpipe.pipeline_base.config import Config
from radpipe.pipeline_base.state import State
import radpipe.pipeline_base.error_codes
import logging
# from radpipe.pipeline_base.logger import Logger

# set logging level
LOGGING_LEVEL = logging.INFO
# default place to save cluster job scripts
# (mostly useful for post-mortem debugging)
DEFAULT_JOBSCRIPT_DIR = 'jobscripts'
# default name of the pipeline configuration file
DEFAULT_CONFIG_FILE = 'pipeline.config'


def parse_command_line(version):
    '''Parse the command line arguments of the pipeline'''
    parser = cmdline.get_argparse(description='RAD-Seq pipeline',
        ignored_args = ["version"] )
    parser.add_argument('--config', type=str, default=DEFAULT_CONFIG_FILE,
        help='Pipeline configuration file in YAML format, defaults to {}' \
            .format(DEFAULT_CONFIG_FILE))
    parser.add_argument('--jobscripts', type=str,
        default=DEFAULT_JOBSCRIPT_DIR,
        help='Directory to store cluster job scripts created by the ' \
             'pipeline, defaults to {}'.format(DEFAULT_JOBSCRIPT_DIR))
    parser.add_argument('--version', action='version',
        version='%(prog)s ' + version)
    return parser.parse_args()

def main(program_name, program_version, make_pipeline):
    '''Initialise the pipeline, then run it'''
    # Parse command line arguments
    options = parse_command_line(program_version)
    # Initialise the logger
    # logger = Logger(__name__, options.log_file, options.verbose)
    if options.log_file:
        logging.basicConfig(
            filename=options.log_file,
            level=LOGGING_LEVEL,
            filemode="a",
            format="%(asctime)s %(levelname)s - %(message)s",
            datefmt="%m-%d-%Y %H:%M:%S")
    logger = logging.getLogger(__name__)
    # Log the command line used to run the pipeline
    logger.info("*** radpipe ***")
    logger.info(' '.join(sys.argv))
    drmaa_session = None
    try:
        # Set up the DRMAA session for running cluster jobs
        import drmaa
        drmaa_session = drmaa.Session()
        drmaa_session.initialize()
    except Exception as e:
        print("{progname} error using DRMAA library".format(progname=program_name), file=sys.stdout)
        print("Error message: {msg}".format(msg=e.message, file=sys.stdout))
        exit(error_codes.DRMAA_ERROR)
    # Parse the configuration file, and initialise global state
    config = Config(options.config)
    config.validate()
    state = State(options=options, config=config, logger=logger,
                  drmaa_session=drmaa_session)
    # Build the pipeline workflow
    pipeline = make_pipeline(state)
    # Run (or print) the pipeline
    cmdline.run(options)
    if drmaa_session is not None:
        # Shut down the DRMAA session
        drmaa_session.exit()
