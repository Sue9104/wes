import luigi
import tasks
import os
import pandas as pd

class DoMpileup(luigi.WrapperTask):
    """Bcftools Mpileup calling variants

    Flow: Trimmomatic -> Mapping -> RemoveDuplicates -> CallVariants

    Attributes:
        samples (list): list of sample id
        outdir (str): the parent directory of output

    Todo:
        * remove the time limitation after all mapping is done.
    """
    infile = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        os.makedirs(self.outdir, exist_ok=True)
        sample_infos = pd.read_csv(self.infile, index_col=0)
        # prepare fastq
        yield [ PrepareFastq(sample = sample, outdir = self.outdir,
                             fastq = sample_infos.loc[sample].to_dict())
               for sample in sample_infos.index]
        # quality control
        yield [FastQC(sample = sample, outdir = self.outdir)
               for sample in sample_infos.index]
        # call variants
        yield [Mpileup(sample = sample, outdir = self.outdir)
               for sample in sample_infos.index]



class PrepareFastq(tasks.PrepareFastq):
    """link raw fastq to outdir/raw-data

    Attributes:
        sample (str): sample name
        outdir (str): output directory
        fastq (dict): fastq path, such as '{"R1":"test_R1.fq", "R2": "test_R2.fq"}'
    """

    def requires(self):
        return []

class FastQC(tasks.FastQC):
    """Quality control for sequencing data using FastQC

    Website: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """
    def requires(self):
        return PrepareFastq(sample = self.sample, outdir = self.outdir, fastq = '{"":""}')

class MultiQC(tasks.MultiQC):
    """Aggregate results across many samples into a single report

    Website: https://multiqc.info/

    Attributes:
        indir (str): input directory
        outdir (str): output directory,
    """

    def requires(self):
        return []


class Trimmomatic(tasks.Trimmomatic):
    """Cut adapter and remove low quality bases

    Website: http://www.usadellab.org/cms/?page=trimmomatic

    """
    def requires(self):
        return PrepareFastq(sample = self.sample, outdir = self.outdir, fastq = '{"":""}')


class BWA(tasks.BWA):
    """Mapping to genome using BWA"""
    def requires(self):
        return Trimmomatic(sample=self.sample, outdir = self.outdir)

class RemoveDuplicates(tasks.RemoveDuplicates):
    """Remove pcr duplicates using samtools rmdup"""
    def requires(self):
        return BWA(sample=self.sample, outdir = self.outdir)

class Mpileup(tasks.Mpileup):
    """Call variants using freebayes"""
    def requires(self):
        return RemoveDuplicates(sample=self.sample, outdir = self.outdir)
