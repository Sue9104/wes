import luigi
import subprocess
import logging
import os
import glob

class PrepareFastq(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    fastq = luigi.DictParameter(significant=False)
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/raw-data/{sample}_R1.fq.gz".format(
            outdir = self.outdir, sample = self.sample))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'raw-data'), exist_ok=True)
        cmd = "ln -s {R1} {outdir}/raw-data/{sample}_R1.fq.gz\n".format(
            R1=self.fastq["R1"], sample = self.sample, outdir=self.outdir)
        if "R2" in self.fastq:
            cmd += "ln -s {R2} {outdir}/raw-data/{sample}_R2.fq.gz\n".format(
                R2=self.fastq["R2"], sample = self.sample, outdir=self.outdir)
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class FastQC(luigi.Task):
    resources = {"cpu": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget('{}/quality-control/{}_R1_fastqc.html'.format(
            self.outdir, self.sample))

    def run(self):
        outdir_qc = os.path.join(self.outdir, 'quality-control')
        os.makedirs(outdir_qc, exist_ok=True)
        fq1 = "{outdir}/raw-data/{sample}_R1.fq.gz".format(outdir = self.outdir, sample = self.sample)
        # multiple thread parameter: -t
        cmd = "fastqc -t 2 -f fastq -o {outdir} {fq}".format(outdir = outdir_qc, fq = fq1)
        fq2 = fq1.replace('_R1', '_R2')
        if os.path.exists(fq2):
            cmd += " {fq}\n".format(fq = fq2)
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class MultiQC(luigi.Task):
    resources = {"cpu": 6, "memory": 1}
    infile = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget('{outdir}/quality-control/multiqc_report.html'.format(
            outdir = self.outdir))

    def run(self):
        cmd = "multiqc -o {outdir}/quality-control {outdir}".format(outdir = self.outdir)
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class Trimmomatic(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    illumina_adapter = luigi.Parameter()

    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        trimmed_1 = "{}/trim-adapter/{}_1P.fq.gz".format(self.outdir, self.sample)
        trimmed_2 = "{}/trim-adapter/{}_2P.fq.gz".format(self.outdir, self.sample)
        return [luigi.LocalTarget(trimmed_1), luigi.LocalTarget(trimmed_2)]

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'trim-adapter'), exist_ok=True)
        fastq = glob.glob("{}/raw-data/{}_R*.fq.gz".format(self.outdir, self.sample))
        joined_fastq = " ".join(fastq)
        cmd = """
trimmomatic PE -threads 2 -phred33 \
    -trimlog {outdir}/trim-adapter/{sample}.trimmomatic.log \
    -summary {outdir}/trim-adapter/{sample}.trimmomatic.summary \
    {fastq} -baseout {outdir}/trim-adapter/{sample}.fq.gz \
    ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    2> {outdir}/trim-adapter/{sample}_trimmomatic.stdErr
        """.format(
            fastq = joined_fastq, adapter = self.illumina_adapter,
            sample = self.sample, outdir = self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class Qualimap(luigi.Task):
    """Evaluate alignment data: bamqc"""
    resources = {"cpu": 4, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}/qualimapReport.html".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = "qualimap bamqc -bam {outdir}/mapping/{sample}.sorted.bam -nt 8 \
            -outdir {outdir}/mapping/{sample} -outformat PDF:HTML".format(
                outdir = self.outdir, sample = self.sample
        )
        logging.info(cmd)
        subprocess.run(cmd, shell = True)
