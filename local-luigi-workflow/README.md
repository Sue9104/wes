# Whole Exome Analysis

## Environment

conda env create -n wes --file ./wes_environment.yaml

## Run Pipelines

1. run luigid server in backgroud

```
   sudo makedir /etc/luigi && cp files/luigid/luigi.cfg /etc/luigi/luigi.cfg
   luigid --port 8082 --pidfile pid --logdir log --background
```


2. copy and edit "files/cfg/luigi.cfg" to running path


3. set sample_fastq.csv

- first line must be "#SAMPLE,LANE,R1,R2"
- automatically generate sample_fastq.csv, directory structure must be \
"sample_name|info_Laneid|*fq.gz"

```
   python tools/generate_fastq_path_for_specified_folder.py D210038
```


4. run main pipeline


```
   python main.py --workers 20 sample_fastq.csv
```

5. run submodule

```
   PYTHONPATH=~/Pipelines/wes/local-luigi-workflow/ luigi --module tasks Annovar --vcf vcf_name
```
