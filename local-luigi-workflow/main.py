import os
import sys
import argparse
import time
import luigi
import logging
import subprocess
#from workflow import DoFreebayes, DoMpileup, DoGATK
from tasks import DoWES

# Command line parameters
today = time.strftime("%Y-%m-%d", time.localtime())
parser = argparse.ArgumentParser( description='Human Whole Exome Analysis')
parser.add_argument('infile', type=str,
                    help='input csv file, header must be #SAMPLE,LANE,R1,R2, sample must be internal-id_bgi-id')
parser.add_argument('--outdir', type=str, default='wes_' + today,
                    help='output directory')
parser.add_argument('--workers', type=int, default=8,
                    help='how many workers use to run WES, default is 6')
parser.add_argument('--local-scheduler', action='store_true', default=False,
                    help="use local scheduler instead of central scheduler")

args = parser.parse_args()
infile = os.path.abspath(args.infile)
outdir = os.path.abspath(args.outdir)
os.makedirs(outdir, exist_ok=True)

# Logging
os.makedirs(outdir, exist_ok=True)
logfile_handler = logging.FileHandler("{}/wes_{}.log".format(outdir, today))
standout_handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(
    level=logging.INFO,
    handlers = [standout_handler, logfile_handler],
    format='%(asctime)s - [line:%(lineno)d] - %(levelname)s: %(message)s'
)
commit_sha = subprocess.check_output('git rev-parse HEAD', shell=True, encoding='utf-8',
                                     cwd=os.path.dirname(os.path.abspath(__file__)).rstrip())
cmd = 'python {} {} --outdir {} --workers {}'.format(__file__, infile, outdir, args.workers)
logging.info(cmd)
logging.info("Commit sha: {}".format(commit_sha))

# Luigi run WES
#luigi.build([DoFreebayes(infile=infile,outdir=outdir),
#             DoMpileup(infile=infile,outdir=outdir),
#             DoGATK(infile=infile,outdir=outdir)],
#            workers=args.workers, local_scheduler=args.local_scheduler)
luigi.build(
    [
        DoWES(infile=infile,outdir=outdir)
    ],
    workers=args.workers, local_scheduler=args.local_scheduler
)
