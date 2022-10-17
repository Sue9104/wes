import re
import os
import sys
import argparse
import time
import luigi
import logging
import logging.config
import yaml
import subprocess
from workflow import DoWES

# Command line parameters
today = time.strftime("%Y%m%d", time.localtime())
parser = argparse.ArgumentParser( description='Human Whole Exome Analysis')
parser.add_argument('infile', type=str,
                    help='input csv file, header must be sample,lane,r1,r2')
parser.add_argument('--outdir', type=str,
                    help='output directory')
parser.add_argument('--project', type=str, default="ngs",
                    help='project name')
parser.add_argument('--workers', type=int, default=20,
                    help='how many workers use to run WES, default is 20')
parser.add_argument('--local-scheduler', action='store_true', default=False,
                    help="use local scheduler instead of central scheduler")

args = parser.parse_args()
infile = os.path.abspath(args.infile)
outdir = args.outdir
if not outdir:
    outdir = "{}-CallVariants-{}".format(args.project, today)
outdir = os.path.abspath(outdir)
os.makedirs(outdir, exist_ok=True)
os.chdir(outdir)

# logging
script_path = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(script_path,'logging.yaml'), 'r') as stream:
    config = yaml.load(stream, Loader=yaml.FullLoader)
logging.config.dictConfig(config)
rootlogger = logging.getLogger('root')
cmdlogger = logging.getLogger('cmd')
timelogger = logging.getLogger('ptime')
## write command to logging
commit_sha = subprocess.check_output(
    'git rev-parse HEAD', encoding='utf-8', shell=True,
    cwd=os.path.dirname(os.path.abspath(__file__)))
cmd = 'python {} {} --outdir {} --workers {}'.format(
    __file__, infile, outdir, args.workers
)
rootlogger.info(cmd)
rootlogger.info("Commit sha: {}".format(commit_sha))
cmdlogger.info(cmd)

# Luigi run WES
luigi.build(
    [DoWES(infile=infile,outdir=outdir)],
    workers=args.workers,
    local_scheduler=args.local_scheduler
)
