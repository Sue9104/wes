import luigi
import subprocess
import logging
import os
from .reference import Reference

class RemoveDuplicates(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.deduped.bam".format(
            sample = self.sample, outdir = self.outdir))

    def run(self):
        cmd = "samtools rmdup {outdir}/mapping/{sample}.sorted.bam {outdir}/mapping/{sample}.deduped.bam".format(
            sample = self.sample, outdir = self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class MarkDuplicate(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        output_bam = "{}/mapping/{}.mark_dedup.bam".format(self.outdir, self.sample)
        output_metrics = "{}/mapping/{}.mark_dedup.metrics.txt".format(self.outdir, self.sample)
        return [luigi.LocalTarget(output_bam), luigi.LocalTarget(output_metrics)]

    def run(self):
        cmd = """
gatk MarkDuplicatesSpark -conf 'spark.executor.cores=3' \
    --optical-duplicate-pixel-distance 2500 \
    --input {outdir}/mapping/{sample}.sorted.bam \
    --output {outdir}/mapping/{sample}.mark_dedup.bam \
    --metrics-file {outdir}/mapping/{sample}.mark_dedup.metrics.txt
    """.format(
            sample = self.sample, outdir = self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class BaseScoreQualityRecalibrator(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.bqsr_recal_data.table".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = """
gatk BaseRecalibrator --reference {genome} --intervals {interval} \
    --known-sites {dbsnp} --known-sites {mills_gold_standard} \
    --use-original-qualities \
    --input {outdir}/mapping/{sample}.mark_dedup.bam \
    --output {output}
        """.format(
            genome = Reference().genome, interval = Reference().interval,
            dbsnp = Reference().dbsnp, mills_gold_standard = Reference().mills_gold_standard,
           outdir = self.outdir, sample = self.sample, output = self.output().path
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class ApplyBQSR(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.bqsr.bam".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = "gatk ApplyBQSR --java-options '-Xms3g' \
               --use-original-qualities --reference {genome} \
               --input {outdir}/mapping/{sample}.mark_dedup.bam \
               --bqsr-recal-file {outdir}/mapping/{sample}.bqsr_recal_data.table \
               --output {output}".format(
                   genome = Reference().genome, outdir = self.outdir,
                   sample = self.sample, output = self.output().path
               )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class AnalyzeCovariates(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    before_recal_table = luigi.Parameter()
    after_recal_table = luigi.Parameter()
    output_plot = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget(self.output_plot)

    def run(self):
        cmd = "gatk AnalyzeCovariates \
            -before {before_recal} -after {after_recal} -plots {output}".format(
            before_recal = self.before_recal_table,
            after_recal = self.after_recal_table,
            output = self.output_plot,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


