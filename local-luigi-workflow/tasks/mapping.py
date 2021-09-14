import luigi
import os
import glob
import logging
import subprocess
from .reference import Reference

class BWA(luigi.Task):
    resources = {"cpu": 10, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.sorted.bam".format(
            outdir = self.outdir, sample = self.sample))

    def run(self):
        outdir_mapping = os.path.join(self.outdir, 'mapping')
        os.makedirs(outdir_mapping, exist_ok=True)
        fastq = glob.glob("{outdir}/trim-adapter/{sample}_*P.fq.gz".format(
            outdir = self.outdir, sample = self.sample
        ))
        joined_fastq = " ".join(fastq)
        cmd = """
sample={sample} && bwa mem -t 20 \
  -R "@RG\\tID:$sample\\tSM:$sample\\tLB:WES\\tPL:Illumina" {genome} {fq} \
  | samtools sort -@ 2 -o {output} -
        """.format(
            sample = self.sample, fq = joined_fastq,
            genome = Reference().genome, output = self.output().path
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

