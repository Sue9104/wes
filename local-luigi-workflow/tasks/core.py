import luigi
from luigi.util import inherits
import pandas as pd
import numpy as np
import subprocess
import os
import logging
import time
from .reference import Reference


##################################################
## Callbacks for execution
##################################################
@luigi.Task.event_handler(luigi.Event.START)
def start_time(task):
    logging.warning('Task:{task}\t StartTime:{time}'.format(
        task = task.__class__.__name__,
        time = time.strftime('%Y-%m-%d %H:%M', time.localtime())
    ))

@luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
def execution_time(task, processing_time):
    logging.warning('Task:{task}\t ProcessingTime:{time:.0f} min'.format(
        task = task.__class__.__name__,
        time = processing_time / 60
    ))


##################################################
## Class Definition
##################################################
class Input(object):
    """Input Parameter for analysis

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    Properties:
        samples (list): a list of sample name from infile
        info (dict): a detail dict for sample information, for example:
            {"sample_name": named_turple("LANE", "R1", "R2")}

    """

    def __init__(self, infile, outdir):
        """object initialize"""
        self.infile = infile
        self.outdir = outdir

    @property
    def samples(self):
        sample_info_data = pd.read_csv(self.infile)
        samples = list( set(sample_info_data["#SAMPLE"].values) )
        return samples

    @property
    def info(self):
        info = {}
        sample_info_data = pd.read_csv(self.infile)
        for row in sample_info_data.itertuples(index=False, name="Info"):
            info[row.sample] = \
                info[row.sample] + [row] if row.sample in info else [row]
        return info


@inherits(Input)
class QualityControl(luigi.WrapperTask):
    """Using FastQC to quality control

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    Properties:
        outdir_qc (str): {outdir}/quality-control

    Raise:
        NotImplementedError: class need the requires function!!!
    """

    resources = {"cpu": 6}

    @property
    def outdir_qc(self):
        outdir_qc = os.path.join(self.outdir, 'quality-control')
        os.makedirs(outdir_qc, exist_ok=True)
        return outdir_qc

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget('{}/multiqc_report.html'.format(self.outdir_qc))

    def run(self):
        cmd = "multiqc -o {0} {0}".format(self.outdir_qc)
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class CheckExists(luigi.ExternalTask):
    """Check file if exists or not"""
    resources = {"cpu": 1}
    file = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.file)


class FastQC(luigi.Task):
    """quality control for fastq file

    Attributes:
        fastq (str): fastq file path
        outdir(str): output directory

    Yield:
        - _fastqc.zip
        - _fastqc.html

    Raise:
        NotImplementedError: class need the requires function!!!
    """

    resources = {"cpu": 1}
    fastq = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        out_prefix = os.path.basename(self.fastq).split('.')[0]
        return luigi.LocalTarget('{}/{}_fastqc.html'.format(self.outdir, out_prefix))

    def run(self):
        cmd = "fastqc -f fastq -o {} {}".format(self.outdir_qc, self.fastq)
        logging.info(cmd)
        subprocess.run(cmd, shell=True)



@inherits(Input)
class FiltLowQuality(luigi.WrapperTask):
    resources = {"cpu": 6, "memory": 1}

    def requires(self):
        outdir_trim_adapter = os.path.join(self.outdir, 'trim-adapter')
        os.makedirs(outdir_trim_adapter, exist_ok=True)
        sample_infos = pd.read_csv(self.infile, index_col=[0,1])
        return [Trimmomatic(sample=sample, lane = lane,
                           R1=sample_infos.loc[sample]["R1"],
                           R2=sample_infos.loc[sample]["R2"],
                           outdir_trim_adapter=outdir_trim_adapter)
               for sample,lane in sample_infos.index]

class Trimmomatic(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    illumina_adapter = luigi.Parameter()

    sample = luigi.Parameter()
    lane = luigi.Parameter()
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    outdir_trim_adapter = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        trimmed_1 = "{}/{}_{}_1P.fq.gz".format(self.outdir_trim_adapter, self.sample, self.lane)
        trimmed_2 = "{}/{}_{}_2P.fq.gz".format(self.outdir_trim_adapter, self.sample, self.lane)
        return [luigi.LocalTarget(trimmed_1), luigi.LocalTarget(trimmed_2)]

    def run(self):
        os.makedirs(self.outdir_trim_adapter, exist_ok=True)
        prefix = "{}/{}_{}.fq.gz".format(self.outdir_trim_adapter, self.sample, self.lane)
        cmd = "trimmomatic PE -threads 2 -phred33 -trimlog {outdir}/trimmomatic.log -summary {outdir}/trimmomatic.summary {R1} {R2} -baseout {prefix} ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
            R1 = self.R1, R2=self.R2, prefix = prefix, adapter = self.illumina_adapter, outdir = self.outdir_trim_adapter
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

@inherits(Input)
class Mapping(luigi.WrapperTask):
    resources = {"cpu": 6, "memory": 1}

    def requires(self):
        sample_infos = pd.read_csv(self.infile, index_col=[0,1])
        return [ MergeSampleBWA(infile=self.infile, outdir=self.outdir, sample=sample)
               for sample, value in sample_infos.groupby('#SAMPLE') ]

@inherits(Input)
class MergeSampleBWA(luigi.Task):
    resources = {"cpu": 2, "memory": 1}
    sample = luigi.Parameter()
    bamdst = luigi.Parameter()

    def requires(self):
        sample_infos = pd.read_csv(self.infile, index_col=[0,1])
        return [BWA(infile=self.infile, sample=self.sample, lane=index[1], outdir=self.outdir)
               for index in sample_infos.groupby('#SAMPLE').get_group(self.sample).index]

    def output(self):
        return luigi.LocalTarget("{}/mapping/{}.merged.bam".format(self.outdir, self.sample))

    def run(self):
        sample_infos = pd.read_csv(self.infile, index_col=[0,1])
        lanes =[ index[1] for index in sample_infos.groupby('#SAMPLE').get_group(self.sample).index]
        os.makedirs( os.path.join(self.outdir, 'quality-control', self.sample), exist_ok=True)
        cmd = 'samtools merge - {bams} | tee {outdir}/mapping/{sample}.merged.bam | samtools index - {outdir}/mapping/{sample}.merged.bam.bai && {bamdst} -p {bed} -o {outdir}/quality-control/{sample}/ {outdir}/mapping/{sample}.merged.bam'.format(
            bams = ' '.join([ '{}/mapping/{}_{}.sorted.bam'.format(self.outdir, self.sample, lane) for lane in lanes]),
            outdir = self.outdir,
            sample = self.sample,
            bamdst = self.bamdst,
            bed = Reference().bed
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

@inherits(Input)
class BWA(luigi.Task):
    resources = {"cpu": 15, "memory": 10}
    bamdst = luigi.Parameter()

    sample = luigi.Parameter()
    lane = luigi.Parameter()

    def requires(self):
        sample_infos = pd.read_csv(self.infile, index_col=[0,1])
        return Trimmomatic(sample=self.sample, lane=self.lane, R1=sample_infos.loc[self.sample, self.lane]["R1"], R2=sample_infos.loc[self.sample, self.lane]["R2"],
                           outdir_trim_adapter=os.path.join(self.outdir,'trim-adapter'))

    def output(self):
        return luigi.LocalTarget("{}/mapping/{}_{}.sorted.bam".format(self.outdir, self.sample, self.lane))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'mapping'), exist_ok=True)
        cmd = 'bwa mem -t 20 -R "@RG\\tID:{lane}\\tSM:{sample}\\tLB:WES\\tPL:Illumina" {genome} {outdir}/trim-adapter/{sample}_{lane}_1P.fq.gz {outdir}/trim-adapter/{sample}_{lane}_2P.fq.gz | samtools sort -@ 2 -o {outdir}/mapping/{sample}_{lane}.sorted.bam - && \mkdir {outdir}/quality-control/{sample}_{lane}/ && {bamdst} -p {bed} -o {outdir}/quality-control/{sample}_{lane} {outdir}/mapping/{sample}_{lane}.sorted.bam'.format(
            sample=self.sample, lane=self.lane, genome=Reference().genome, outdir=self.outdir,
            bamdst=self.bamdst, bed=Reference().bed
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class RemoveDuplicates(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{}/mapping/{}.deduped.bam".format(self.outdir, self.sample))

    def run(self):
        cmd = "samtools rmdup {outdir}/mapping/{sample}.merged.bam {outdir}/mapping/{sample}.deduped.bam".format(
            outdir=self.outdir, sample=self.sample
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

@inherits(Input)
class CallVariants(luigi.Task):
    resources = {"cpu": 1, "memory": 1}

    def requires(self):
        return Mapping(infile=self.infile, outdir=self.outdir)

    def output(self):
        sample_infos = pd.read_csv(self.infile)
        sample_list = list( set(sample_infos["#SAMPLE"].values) )
        return [
            luigi.LocalTarget(
                "{}/call-variants/{}.{}.gatk.raw.vcf.gz".format(
                    self.outdir, sample, Reference().genome_version
                )
            )
            for sample in sample_list
        ]

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        sample_infos = pd.read_csv(self.infile)
        sample_list = list( set(sample_infos["#SAMPLE"].values) )
        #yield Freebayes(samples=sample_list, outdir=self.outdir)
        #yield Mpileup(samples=sample_list, outdir=self.outdir)
        yield GATK(samples=sample_list, outdir=self.outdir)


@inherits(Input)
class DataProcessingBeforeCalling(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    samples = luigi.ListParameter()

    def requires(self):
        return [ApplyBQSR(sample=sample, infile=self.infile, outdir=self.outdir)
               for sample in self.samples]

    def output(self):
        return [luigi.LocalTarget("{}/mapping/{}.bqsr.bam".format(self.outdir, sample))
                for sample in self.samples]


class MarkDuplicate(luigi.Task):
    resources = {"cpu": 20, "memory": 18}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.mark_dedup.bam".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd ="""
gatk --java-options '-Xmx16G' MarkDuplicatesSpark -conf 'spark.executor.cores=8' \
  --input {outdir}/mapping/{sample}.merged.bam \
  --output {outdir}/mapping/{sample}.mark_dedup.bam \
  --metrics-file {outdir}/mapping/{sample}.mark_dedup.metrics.txt
        """.format(
            genome = Reference().genome,
            outdir = self.outdir, sample = self.sample,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class BaseScoreQualityRecalibrator(luigi.Task):
    resources = {"cpu": 4, "memory": 4}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return MarkDuplicate(sample=self.sample, outdir=self.outdir)

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.bqsr_recal_data.table".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = "gatk BaseRecalibrator --java-options '-Xmx16G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal' \
               --reference {genome} \
               --intervals {interval} \
               --known-sites {mills_gold_standard} \
               --known-sites {snp_1000g} \
               --known-sites {indel_1000g} \
               --known-sites {hapmap} \
               --known-sites {omni} \
               --input {outdir}/mapping/{sample}.mark_dedup.bam \
               --output {output}".format(
                   genome = Reference().genome,
                   interval = Reference().interval,
                   dbsnp = Reference().dbsnp,
                   mills_gold_standard = Reference().mills_gold_standard,
                   snp_1000g = Reference().snp_1000g,
                   indel_1000g = Reference().indel_1000g,
                   hapmap = Reference().hapmap,
                   omni = Reference().omni,
                   outdir = self.outdir,
                   sample = self.sample,
                   output = self.output().path
               )
        dbsnp = Reference().dbsnp
        phase3_1000g = Reference().phase3_1000g
        if dbsnp != "None":
            cmd += " --known-sites " + dbsnp
        if phase3_1000g != "None":
            cmd += " --known-sites " + phase3_1000g

        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class ApplyBQSR(luigi.Task):
    resources = {"cpu": 2, "memory": 20}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return BaseScoreQualityRecalibrator(sample=self.sample, outdir=self.outdir)

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.bqsr.bam".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = "gatk ApplyBQSR --java-options '-Xmx16g -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10' \
               --reference {genome} \
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
    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.bqsr_analyze_covariates.pdf".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = "gatk AnalyzeCovariates \
               -before {outdir}/mapping/{sample}.bqsr_recal_data.table \
               -after {outdir}/mapping/{sample}.bqsr_recal_data.after.table \
               -plots {outdir}/mapping/{sample}.bqsr_analyze_covariates.pdf".format(
                   outdir = self.outdir, sample = self.sample
               )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class GATK(luigi.WrapperTask):
    """GATK Calling Variants

    ClassFlow: Trimmomatic -> Mapping -> MarkDuplicate -> BaseScoreQualityRecalibrator -> ApplyBQSR -> HaplotypeCaller
    Process: Remove adapter and filt low quality bases -> Mapping to genome -> Mark duplicate -> Recalibrate base score quality -> Call variants
    Filter Default: None
    Filter Options:
          HaplotypeCaller -> CNNFilt -> MergeVcfs
          HaplotypeCaller -> HardFilt -> MergeVcfs

    Attributes:
        samples (list): List of sample id
        outdir (int, optional): Parent directory of output

    Todo:
        * Compare hard filter and CNN
    """
    resources = {"cpu": 1, "memory": 16}
    samples = luigi.ListParameter()
    outdir = luigi.Parameter()

    def requires(self):
        #tasks = [HardFilt(sample=sample, outdir=self.outdir) for sample in self.samples]
        #tasks += [CNNFilt(sample=sample, outdir=self.outdir) for sample in self.samples]
        tasks = [HaplotypeCaller(sample=sample, outdir=self.outdir) for sample in self.samples]
        return tasks

    def output(self):
        return [
            luigi.LocalTarget(
                "{outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz".format(
                    outdir = self.outdir,
                    sample = sample,
                    version = Reference().genome_version
                )
            ) for sample in self.samples
        ]


class MergeVcfs(luigi.Task):
    vcfs = luigi.ListParameter()
    out = luigi.Parameter()
    def requires(self):
        return []
    def output(self):
        return luigi.LocalTarget(self.out)
    def run(self):
        cmd = "gatk MergeVcfs {vcfs} -O {out}".format(
            vcfs = ' '.join([" -I {} ".format(vcf) for vcf in self.vcfs]),
            out = self.out
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class HaplotypeCaller(luigi.Task):
    resources = {"cpu": 10, "memory": 16}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return ApplyBQSR(sample=self.sample, outdir=self.outdir)

    def output(self):
        return luigi.LocalTarget(
            "{outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz".format(
                outdir = self.outdir,
                sample = self.sample,
                version = Reference().genome_version
            )
        )

    def run(self):
        cmd = """
gatk --java-options '-Xmx16g' HaplotypeCaller -G StandardAnnotation -G StandardHCAnnotation \
    --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads 10 \
    -R {genome} --dbsnp {dbsnp} \
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
    -I {outdir}/mapping/{sample}.bqsr.bam \
    -O {outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz \
    -bamout {outdir}/call-variants/{sample}.{version}.haplotypecaller.bam
gatk CollectHsMetrics -R {genome} \
    --BAIT_INTERVALS {interval} --TARGET_INTERVALS {interval} \
    -I {outdir}/mapping/{sample}.bqsr.bam \
    -O {outdir}/mapping/{sample}.{version}.hs.metrics.txt \
    """.format(
            outdir = self.outdir, sample = self.sample,
            genome = Reference().genome,
            version = Reference().genome_version,
            dbsnp = Reference().dbsnp,
            interval = Reference().interval
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class CNNFilt(luigi.Task):
    resources = {"cpu": 4, "memory": 16}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return HaplotypeCaller(sample=self.sample, outdir=self.outdir)

    def output(self):
        return luigi.LocalTarget(
            "{outdir}/call-variants/{sample}.{version}.gatk.cnn_filt.vcf.gz".format(
                outdir = self.outdir, sample = self.sample,
                version = Reference().genome_version
            )
        )

    def run(self):
        cmd = """
gatk --java-options '-Xmx16G' CNNScoreVariants -tensor-type read_tensor -R {genome} \
    -I {outdir}/mapping/{sample}.bqsr.bam \
    -V {outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.cnn_mark.vcf.gz
gatk FilterVariantTranches --info-key CNN_2D \
    --resource {hapmap} --resource {mills_gold_standard} \
    --snp-tranche 99.95 --indel-tranche 99.4 --invalidate-previous-filters \
    -V {outdir}/call-variants/{sample}.{version}.gatk.cnn_mark.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.cnn_filt.vcf.gz
        """.format(
            outdir = self.outdir, sample = self.sample,
            genome = Reference().genome,
            version = Reference().genome_version,
            hapmap = Reference().hapmap, mills_gold_standard = Reference().mills_gold_standard,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class HardFilt(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return HaplotypeCaller(sample=self.sample, outdir=self.outdir)

    def output(self):
        return luigi.LocalTarget(
            "{outdir}/call-variants/{sample}.{version}.gatk.hard_filt.vcf.gz".format(
                outdir = self.outdir, sample = self.sample,
                version = Reference().genome_version
            )
        )

    def run(self):
        cmd = """
gatk SelectVariants --select-type-to-include SNP -R {genome} \
    -V {outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.raw.snp.vcf.gz
gatk VariantFiltration -R {genome} \
    --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || QUAL < 30.0 || SOR > 3.0' \
    --filter-name 'hard_filter' \
    -V {outdir}/call-variants/{sample}.{version}.gatk.raw.snp.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.raw.snp.hard_filt.vcf.gz
gatk SelectVariants --select-type-to-include INDEL -R {genome} \
    -V {outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.raw.indel.vcf.gz
gatk VariantFiltration -R {genome} \
    --filter-expression 'QD < 2.0 || FS > 200.0 || QUAL < 30.0' \
    --filter-name 'hard_filter' \
    -V {outdir}/call-variants/{sample}.{version}.gatk.raw.indel.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.raw.indel.hard_filt.vcf.gz
gatk MergeVcfs \
    -I {outdir}/call-variants/{sample}.{version}.gatk.raw.snp.hard_filt.vcf.gz \
    -I {outdir}/call-variants/{sample}.{version}.gatk.raw.indel.hard_filt.vcf.gz \
    -O {outdir}/call-variants/{sample}.{version}.gatk.hard_filt.vcf.gz
rm {outdir}/call-variants/{sample}.{version}.gatk.raw.snp.vcf.gz* \
   {outdir}/call-variants/{sample}.{version}.gatk.raw.indel.hard_filt.vcf.gz*
        """.format(
            outdir = self.outdir, sample = self.sample,
            genome = Reference().genome,
            version = Reference().genome_version
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class Freebayes(luigi.Task):
    """Freebayes calling variants

    Flow: Trimmomatic -> Mapping -> RemoveDuplicates -> CallVariants

    Attributes:
        samples (list): list of sample id
        outdir (str): the parent directory of output

    Todo:
        * remove the time limitation after all mapping is done.
    """
    resources = {"cpu": 1, "memory": 1}
    samples = luigi.ListParameter()
    outdir = luigi.Parameter()

    def requires(self):
        return [RemoveDuplicates(sample=sample, outdir=self.outdir)
                for sample in self.samples]

    def output(self):
        return luigi.LocalTarget("{}/call-variants/freebayes.raw.vcf".format(self.outdir))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        joined_dedup_bam = ' '.join(["{}/mapping/{}.deduped.bam".format(self.outdir, sample)
                                     for sample in self.samples])
        cmd = "freebayes -f {genome} {bam} > {outdir}/call-variants/freebayes.raw.vcf".format(
            genome=Reference().genome, bam=joined_dedup_bam, outdir=self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class Mpileup(luigi.Task):
    """Bcftools Mpileup calling variants

    Flow: Trimmomatic -> Mapping -> RemoveDuplicates -> CallVariants

    Attributes:
        samples (list): list of sample id
        outdir (str): the parent directory of output

    Todo:
        * remove the time limitation after all mapping is done.
    """
    resources = {"cpu": 1, "memory": 1}
    samples = luigi.ListParameter()
    outdir = luigi.Parameter()

    def requires(self):
        return [RemoveDuplicates(sample=sample, outdir=self.outdir)
                for sample in self.samples]

    def output(self):
        return luigi.LocalTarget("{}/call-variants/mpileup.raw.vcf".format(self.outdir))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        joined_dedup_bam = ' '.join(["{}/mapping/{}.deduped.bam".format(self.outdir, sample)
                                     for sample in self.samples])
        cmd = "bcftools mpileup -f {genome} {bam}| bcftools call -mv --ploidy {ploidy} -o {outdir}/call-variants/mpileup.raw.vcf".format(
            genome=Reference().genome, ploidy=Reference.genome_version, bam=joined_dedup_bam, outdir=self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class Annovar(luigi.Task):
    """VCF Annotation using Annovar

    Flow: vcf -> avinput_1(begin_with_chr1) -> avinput_2(begin_with_1) -> anotated_csv

    Attributes:
        vcf (str): vcf
        outfile (str): default vcf.hg19_multianno.csv

    """
    resources = {"cpu": 2, "memory": 2}
    vcf = luigi.Parameter()
    annovar_software_dir = luigi.Parameter()
    annovar_database_dir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{}.hg19_multianno.final.csv".format(self.vcf))

    def run(self):
        cmd = """
perl {annovar_software_dir}/convert2annovar.pl --format vcf4  {vcf} --outfile {vcf}.avinput --withzyg
sed -i 's/^chr//g' {vcf}.avinput
perl {annovar_software_dir}/table_annovar.add_spliceai_exome.pl -buildver hg19 -protocol refGene,ensGene,cytoBand,phastConsElements100way,tfbsConsSites,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,avsnp150,snp138NonFlagged,popfreq_all_20150413,popfreq_max_20150413,gnomad211_exome,dbnsfp41a,gerp++gt2,gerp++elem,dbscsnv11,clinvar_20210501,ipmch592,fudan75,ClinPred,hgmd,intervar_20180118,spliceai_exome -operation g,g,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -csvout -polish -remove -out {vcf} {vcf}.avinput {annovar_database_dir}
xsv fmt -t "\t" {vcf}.hg19_multianno.csv | awk 'BEGIN{{FS="\t";OFS="\t"}}{{if (($122=="Benign") || ((($27>0.05) || ($33>0.05) || ($46>0.05)) && ($122="."))){{print $0}}}}' | xsv select Chr,Start,End,Ref,Alt,1000G_ALL,ExAC_ALL,AF -d "\t" > {vcf}.hg19_multianno.dropped_sites.csv
xsv join --left 1-5 {vcf}.hg19_multianno.csv 1-5 {vcf}.hg19_multianno.dropped_sites.csv | xsv fmt -t "\t" | awk 'BEGIN{{FS"\t";OFS="\t"}}{{if ($176==""){{print $0}}}}' | xsv select -d "\\t" 1-175 | uniq  > {vcf}.hg19_multianno.final.csv.tmp
head -1 {vcf}.hg19_multianno.csv | cat - {vcf}.hg19_multianno.final.csv.tmp > {vcf}.hg19_multianno.final.csv
rm {vcf}.hg19_multianno.final.csv.tmp
        """.format(
            annovar_software_dir = self.annovar_software_dir,
            annovar_database_dir = self.annovar_database_dir,
            vcf = self.vcf
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


@inherits(Input)
class Annotation(luigi.Task):
    resources = {"cpu": 2, "memory": 2}

    def requires(self):
        return CallVariants(infile=self.infile, outdir=self.outdir)

    def output(self):
        sample_infos = pd.read_csv(self.infile)
        sample_list = list( set(sample_infos["#SAMPLE"].values) )
        return [
            luigi.LocalTarget(
                "{}/call-variants/{}.{}.gatk.raw.vcf.gz.hg19_multianno.final.csv".format(
                    self.outdir, sample, Reference().genome_version
                )
            )
            for sample in sample_list
        ]

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        sample_infos = pd.read_csv(self.infile)
        sample_list = list( set(sample_infos["#SAMPLE"].values) )
        yield [
            Annovar(
                vcf = "{}/call-variants/{}.{}.gatk.raw.vcf.gz".format(
                    self.outdir, sample, Reference().genome_version
                )
            ) for sample in sample_list
        ]