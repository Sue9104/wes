import re
import os
import glob
import pandas as pd
import argparse
parser = argparse.ArgumentParser( description='Generate Fastq Path for WES Analysis')
parser.add_argument( 'project', type=str, help='project id, as D210001')
parser.add_argument( '--indir', type=str, help='fastq parent folder',
                    default='/home/junyu/exome')
args = parser.parse_args()

data_path = "{}/{}".format(os.path.abspath(args.indir), args.project)
samples = sorted([
    d for d in os.listdir("{}/raw-data".format(data_path))
    if os.path.isdir("{}/raw-data/{}".format(data_path, d))
])
rows = []
for sample in samples:
    print("Sample: " + sample)
    sample_path = "{}/raw-data/{}".format(data_path, sample)
    for lane in sorted(os.listdir(sample_path)):
        print("Lane: " + lane)
        lane_id = re.search("_(?P<lane>L[0-9]+)", lane).group('lane')
        fastqs = sorted(glob.glob(
            '{}/raw-data/{}/{}/*q.gz'.format(data_path,sample,lane)
        ))
        index = ['#SAMPLE','LANE','R1','R2']
        internal_id = '{}_{}'.format(args.project, sample)
        row = pd.Series(
            [internal_id, lane_id] + fastqs,
            index=index
        )
        rows.append(row)
outfile = "{}/{}.fastq_path.csv".format(data_path, args.project)
print("OUTFILE: " + outfile)
pd.DataFrame(rows).to_csv(outfile, index=False)

