import os
import luigi
import pandas as pd
import tasks


class DoGATK(luigi.WrapperTask):
    """GATK Calling Variants

    ClassFlow: Trimmomatic -> Mapping -> MarkDuplicate -> BaseScoreQualityRecalibrator -> ApplyBQSR -> HaplotypeCaller -> HardFilter or CNN
    Process: Remove adapter and filt low quality bases -> Mapping to genome -> Mark duplicate -> Recalibrate base score quality -> Call variants

    Attributes:
        infile (str): input sample fastq
        outdir (str): output directory

    Todo:
        * Compare hard filter and CNN
    """
    resources = {"cpu": 1, "memory": 1}
    infile = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        os.makedirs(self.outdir, exist_ok=True)
        sample_infos = pd.read_csv(self.infile, index_col=0)
        # workflow for single sample
        for sample in sample_infos.index:
            # prepare fastq
            yield PrepareFastq(sample = sample, outdir = self.outdir,
                               fastq = sample_infos.loc[sample].to_dict())
            # quality control
            yield FastQC(sample = sample, outdir = self.outdir)
            yield Qualimap(sample = sample, outdir = self.outdir)
            # call variants
            yield HaplotypeCaller(sample = sample, outdir = self.outdir)
            vcf = "{outdir}/call-variants/{sample}.gatk.raw.vcf".format(
                outdir = self.outdir, sample = sample
            )
            # annovar
            anno_prefix = "{outdir}/annotation/{sample}.gatk.annovar".format(
                outdir = self.outdir, sample = sample
            )
            #yield Annovar(sample = sample, outdir = self.outdir,
            #              input_vcf = vcf, output_prefix = anno_prefix)
            # snpeff
            #yield SnpEff(sample = sample, outdir = self.outdir,
            #             input_vcf = vcf, output_prefix = "gatk")
        yield MultiQC(infile = self.infile, outdir = self.outdir)

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
        infile (str): input csv file
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

class MarkDuplicate(tasks.MarkDuplicate):
    """Remove pcr duplicates using samtools rmdup"""
    def requires(self):
        return BWA(sample=self.sample, outdir = self.outdir)

class BaseScoreQualityRecalibrator(tasks.BaseScoreQualityRecalibrator):
    """Generate bqsr recal table"""
    def requires(self):
        return MarkDuplicate(sample=self.sample, outdir = self.outdir)

class ApplyBQSR(tasks.ApplyBQSR):
    """Generate bqsr recal table"""
    def requires(self):
        return BaseScoreQualityRecalibrator(sample=self.sample, outdir = self.outdir)

class HaplotypeCaller(tasks.HaplotypeCaller):
    """Generate raw vcf using HaplotypeCaller"""
    def requires(self):
        return ApplyBQSR(sample=self.sample, outdir = self.outdir)

class HardFilter(tasks.HardFilter):
    """Filt raw vcf """
    def requires(self):
        return HaplotypeCaller(sample=self.sample, outdir = self.outdir, gvcf = False)

class Annovar(tasks.Annovar):
    """Annotation filtered vcf"""
    def requires(self):
        return HaplotypeCaller(sample=self.sample, outdir = self.outdir, gvcf = False)

class SnpEff(tasks.SnpEff):
    """Annotation vcf using snpEff"""
    def requires(self):
        return HaplotypeCaller(sample=self.sample, outdir = self.outdir, gvcf = False)

class Qualimap(tasks.Qualimap):
    """qc for bam"""
    def requires(self):
        return BWA(sample=self.sample, outdir = self.outdir)

class MultiQC(tasks.MultiQC):
    """integrate quality control"""
    def requires(self):
        sample_infos = pd.read_csv(self.infile, index_col=0)
        qc_snpeff = [ SnpEff(sample = sample, outdir = self.outdir,
                             input_vcf = "", output_prefix = "gatk")
                     for sample in sample_infos.index]
        qc_fastqc = [ FastQC(sample = sample, outdir = self.outdir)
                     for sample in sample_infos.index]
        qc_qualimap = [ Qualimap(sample = sample, outdir = self.outdir)
                     for sample in sample_infos.index]
        return qc_snpeff + qc_fastqc + qc_qualimap
