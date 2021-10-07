import luigi
from luigi.util import inherits
import pandas as pd
import numpy as np
import subprocess
import os
import logging
import time
from .reference import Reference
from .tools import Tools


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
        return luigi.LocalTarget(
            '{}/{}_fastqc.html'.format(self.outdir, out_prefix)
        )

    def run(self):
        cmd = "fastqc -f fastq -o {} {}".format(self.outdir_qc, self.fastq)
        logging.info(cmd)
        subprocess.run(cmd, shell=True)



@inherits(Input)
class FiltLowQuality(luigi.WrapperTask):
    """Filt Adapter and Low Quality by Sample

    This is a wrapper task for one sample with multiple lanes

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"
        sample (str): sample name

    Properties:
        outdir_trim_adapter (str): "{outdir}/trim-adapter"

    """
    resources = {"cpu": 6, "memory": 1}
    sample = luigi.Parameter()

    @property
    def outdir_trim_adapter(self):
        outdir_trim_adapter = os.path.join(self.outdir, 'trim-adapter')
        os.makedirs(outdir_trim_adapter, exist_ok=True)
        return outdir_trim_adapter

    def requires(self):
        return [
            Trimmomatic(sample=self.sample, lane = data.lane,
                r1=data.r1, r2=data.r2,
                outdir=self.outdir_trim_adapter)
            for data in self.info[self.sample]
        ]


class Trimmomatic(luigi.Task):
    """Trimmomatic to filt adapter and low quality bases

    Attributes:
        sample (str): sample name
        lane (str): lane name
        r1 (str): lane R1 fastq path
        r2 (str): lane R2 fastq path

    Output:
        * {outdir}/{sample}_{lane}_1P.fq.gz
        * {outdir}/{sample}_{lane}_2P.fq.gz

    """
    resources = {"cpu": 1, "memory": 1}

    sample = luigi.Parameter()
    lane = luigi.Parameter()
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        trimmed_1 = "{}/{}_{}_1P.fq.gz".format(
            self.outdir, self.sample, self.lane)
        trimmed_2 = "{}/{}_{}_2P.fq.gz".format(
            self.outdir, self.sample, self.lane)
        return [luigi.LocalTarget(trimmed_1), luigi.LocalTarget(trimmed_2)]

    def run(self):
        prefix = "{}/{}_{}.fq.gz".format(self.outdir, self.sample, self.lane)
        cmd = """trimmomatic PE -threads 2 -phred33 \
-trimlog {outdir}/trimmomatic.log -summary {outdir}/trimmomatic.summary \
{R1} {R2} -baseout {prefix} ILLUMINACLIP:{adapter}:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36""".format(
            R1 = self.R1, R2=self.R2, prefix = prefix,
            adapter = Reference.illumina_adapter, outdir = self.outdir,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

@inherits(Input)
class Mapping(luigi.WrapperTask):
    """Mapping and Coverage Statistics by Sample

    Workflow: BWA -> MergeSampleBWA -> MappingStatistics

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"
        sample (str): sample name

    Properties:
        outdir_mapping (str): "{outdir}/mapping"

    Output:
        - {outdir}/mapping/{sample}.merged.bam
        - {outdir}/quality-control/{sample}/coverage.report

    """
    resources = {"cpu": 6, "memory": 1}
    sample = luigi.Parameter()

    @property
    def outdir_mapping(self):
        outdir_mapping = os.path.join(self.outdir, 'mapping')
        os.makedirs(outdir_mapping, exist_ok=True)
        return outdir_mapping

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return [
            luigi.LocalTarget(
                "{}/{}.merged.bam".format(self.outdir_mapping, self.sample),
                "{}/quality-control/{}/coverage.report".format(
                    self.outdir, self.sample
                ),
        ]

    def runs(self):
        # mapping by lanes
        yield [
            BWA(
                sample=data.sample, lane=data.lane, outdir=self.outdir_mapping,
                fq1="{outdir}/trim-adapter/{sample}_{lane}_1P.fq.gz".format(
                    sample=data.sample, lane=data.lane,
                ),
                fq2="{outdir}/trim-adapter/{sample}_{lane}_2P.fq.gz".format(
                    sample=data.sample, lane=data.lane,
                ),
            )
            for data in self.info[self.sample]
        ]
        # merge multiple bams
        bams = [
            '{}/{}_{}.sorted.bam'.format(
                self.outdir_mapping, self.sample, data.lane
            ) for data in self.info[self.sample]
        ]
        yield MergeSampleBams(
            sample=self.sample, outdir=self.outdir_mapping,
            bams = bams,
        )
        merged_bam = "{}/{}.merged.bam".format(self.outdir_mapping, self.sample)
        # mapping statistics
        yield MappingStatistics(
            bam = merged_bam,
            outdir = "{}/quality-control/{}".format(self.outdir, self.sample),
        )
        # mark PCR dupliates
        markdup_bam = "{}/{}.markdup.bam".format(self.outdir_mapping,self.sample)
        yield MarkDuplicate(outfile=markdup_bam, bam=merged_bam)
        # recalibrate base quality
        bqsr_bam = "{}/{}.bqsr.bam".format(self.outdir_mapping, self.sample)
        yield BaseScoreQualityRecalibrator(bam=markdup_bam, outfile=bqsr_bam)


@inherits(Input)
class PostMapping(luigi.WrapperTask):
    """Mapping Postprocessing before GATK CallVariants

    Workflow: MarkDuplicateSpark -> BaseScoreQualityRecalibrator

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"
        sample (str): sample name

    Properties:
        outdir_mapping (str): {outdir}/mapping

    Output:
        - {outdir}/mapping/{sample}.bqsr.bam

    """
    resources = {"cpu": 6, "memory": 1}
    sample = luigi.Parameter()

    @property
    def outdir_mapping(self):
        outdir_mapping = os.path.join(self.outdir, 'mapping')
        os.makedirs(outdir_mapping, exist_ok=True)
        return outdir_mapping

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget(
            "{}/{}.bqsr.bam".format(self.outdir_mapping, self.sample),
        )

    def runs(self):
        merged_bam = "{}/{}.merged.bam".format(self.outdir_mapping, self.sample)
        # mark PCR dupliates
        markdup_bam = "{}/{}.markdup.bam".format(self.outdir_mapping,self.sample)
        yield MarkDuplicate(outfile=markdup_bam, bam=merged_bam)
        # recalibrate base quality
        bqsr_bam = "{}/{}.bqsr.bam".format(self.outdir_mapping, self.sample)
        yield BaseScoreQualityRecalibrator(bam=markdup_bam, outfile=bqsr_bam)


class MergeSampleBams(luigi.Task):
    """Merge Multiple Bam Files for One Sample and Coverage Statistics

    Attributes:
        bams (list): a list of bam files
        outfile (str): output bam filename

    Output:
        - {outdir}/mapping/{sample}.merged.bam
    """
    resources = {"cpu": 2, "memory": 1}
    outfile = luigi.Parameter()
    bams = luigi.ListParameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self):
        cmd = """samtools merge - {bams} | tee {outfile} \
| samtools index - {outfile}.bai """.format(
            bams = " ".join(self.bams), outfile = self.outfile,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class MappingStatistics(luigi.Task):
    """Mapping Coverage and Depth Statistics

    Attributes:
        bam (str): bam file
        outdir (str): output directory

    Output:
        - {outdir}/coverage.report
    """
    resources = {"cpu": 2, "memory": 1}
    outdir = luigi.Parameter()
    bam = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(
            "{}/coverage.report".format(
                self.outdir, self.sample
            )
        )

    def run(self):
        os.makedirs(self.outdir, exist_ok=True)
        cmd = "{bamdst} -p {bed} -o {outdir}/ {bam}".format(
            outdir=self.outdir, bam=self.bam,
            bamdst=Tools.bamdst, bed=Reference.bed,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


@inherits(Input)
class BWA(luigi.Task):
    """BWA Mapping

    Attributes:
        sample (str): sample name
        lane (str): lane name
        fq1 (str): fastq1 path
        fq2 (str): fastq2 path
        outdir (str): output directory, existence is required.

    Output:
        * {outdir}/{sample}_{lane}.sorted.bam

    """
    resources = {"cpu": 15, "memory": 10}

    sample = luigi.Parameter()
    lane = luigi.Parameter()
    fq1 = luigi.Parameter()
    fq2 = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(
            "{}/{}_{}.sorted.bam".format(self.outdir, self.sample, self.lane)
        )

    def run(self):
        cmd = """bwa mem -t 20 \
-R "@RG\\tID:{lane}\\tSM:{sample}\\tLB:WES\\tPL:Illumina" {genome} {fq1} {fq2} \
| samtools sort -@ 2 -o {outdir}/{sample}_{lane}.sorted.bam - """.format(
            sample=self.sample, lane=self.lane, fq1=self.fq1, fq2=self.fq2,
            outdir=self.outdir, genome=Reference.genome,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class MarkDuplicate(luigi.Task):
    """MarkDupliate after Mapping

    Attributes:
        bam (str): bam file
        outfile (str): output bam name

    Output:
        - {outdir}/{sample}.mark_dedup.bam
    """

    resources = {"cpu": 20, "memory": 18}
    bam = luigi.Parameter()
    outfile = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self):
        cmd ="""gatk --java-options '-Xmx16G' MarkDuplicatesSpark \
-conf 'spark.executor.cores=8' --input {bam} \
--output {outfile} --metrics-file {outfile}.metrics.txt""".format(
            outfile=self.outdir, bam=self.bam,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class BaseScoreQualityRecalibrator(luigi.Task):
    """BaseScoreQualityRecalibrator after MarkDuplicate

    Attributes:
        bam (str): bam file
        outfile (str): output bam

    """
    resources = {"cpu": 4, "memory": 4}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(
            "{outdir}/{sample}.bqsr_recal_data.table".format(
                outdir = self.outdir, sample = self.sample)
        )

    def run(self):
        cmd = """gatk BaseRecalibrator \
--java-options '-Xmx16G -XX:GCTimeLimit=50 \
-XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal' \
--reference {genome} --intervals {interval} \
--known-sites {mills_gold_standard} --known-sites {snp_1000g} \
--known-sites {dbsnp} --known-sites {indel_1000g} --known-sites {hapmap} \
--known-sites {omni} \
--input {bam} --output {outfile}.recal_table

gatk ApplyBQSR --java-options '-Xmx16g \
-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
-XX:+PrintGCDetails -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10' \
--reference {genome} --input {bam} \
--bqsr-recal-file {outfile}.recal_table --output {outfile}""".format(
            genome=Reference.genome, interval=Reference.interval,
            mills_gold_standard = Reference.mills_gold_standard,
            snp_1000g=Reference.snp_1000g, indel_1000g=Reference.indel_1000g,
            hapmap=Reference.hapmap, omni=Reference.omni, dbsnp=Reference.dbsnp,
            bam=self.bam, outfile=self.outfile,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class AnalyzeCovariates(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{outdir}/mapping/{sample}.bqsr_cov.pdf".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = """gatk AnalyzeCovariates \
-before {outdir}/mapping/{sample}.bqsr_recal_data.table \
-after {outdir}/mapping/{sample}.bqsr_recal_data.after.table \
-plots {outdir}/mapping/{sample}.bqsr_analyze_covariates.pdf""".format(
                   outdir = self.outdir, sample = self.sample
               )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class HaplotypeCaller(luigi.Task):
    """GATK HaplotypeCaller

    Attributes:
        bam (str): bam file
        outfile (str): output vcf name

    Output:
        - {outdir}/call-variants/{sample}.{version}.gatk.raw.vcf.gz
    """

    resources = {"cpu": 10, "memory": 16}
    bam = luigi.Parameter()
    outfile = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self):
        cmd = """gatk --java-options '-Xmx16g' HaplotypeCaller \
-G StandardAnnotation -G StandardHCAnnotation \
--pair-hmm-implementation AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads 10 \
-R {genome} --dbsnp {dbsnp} \
-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
-I {bam} -O {outfile} -bamout {outfile}.bam

gatk CollectHsMetrics -R {genome} \
--BAIT_INTERVALS {interval} --TARGET_INTERVALS {interval} \
-I {outfile} -O {outfile}.hs-metrics.txt """.format(
            genome=Reference.genome, version=Reference.genome_version,
            dbsnp=Reference.dbsnp, interval=Reference.interval,
            bam=self.bam, outfile=self.outfile,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class CNNFilt(luigi.Task):
    """CNNFilt for Raw VCF

    Attributes:
        bam (str): bam file
        raw (str): raw vcf
        outfile (str): output vcf name

    Output:
        - {outdir}/call-variants/{sample}.{version}.gatk.cnn_filt.vcf.gz
    """

    resources = {"cpu": 4, "memory": 16}
    bam = luigi.Parameter()
    raw = luigi.Parameter()
    outfile = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self):
        cmd = """gatk --java-options '-Xmx16G' CNNScoreVariants \
-tensor-type read_tensor -R {genome} \
-I {bam} -V {raw} -O {outfile}.cnn_mark.gz

gatk FilterVariantTranches --info-key CNN_2D \
--resource {hapmap} --resource {mills_gold_standard} \
--snp-tranche 99.95 --indel-tranche 99.4 --invalidate-previous-filters \
-V {outfile}.cnn_mark.gz -O {outfile}""".format(
            genome = Reference.genome, version = Reference.genome_version,
            hapmap = Reference.hapmap,
            mills_gold_standard = Reference.mills_gold_standard,
            bam = self.bam, raw = self.raw, outfile = self.outfile
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class HardFilt(luigi.Task):
    """HardFilt for Raw VCF

    Attributes:
        bam (str): bam file
        raw (str): raw vcf
        outfile (str): output vcf name

    Output:
        - {outdir}/call-variants/{sample}.{version}.gatk.cnn_filt.vcf.gz
    """
    resources = {"cpu": 1, "memory": 1}
    bam = luigi.Parameter()
    raw = luigi.Parameter()
    outfile = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self):
        cmd = """gatk SelectVariants --select-type-to-include SNP -R {genome} \
-V {raw} -O {outfile}.snp.gz
gatk VariantFiltration -R {genome} \
--filter-expression 'QD<2.0 || FS>60.0 || MQ<40.0 || QUAL<30.0 || SOR>3.0' \
--filter-name 'hard_filter' \
-V {outfile}.snp.gz -O {outfile}.snp.hard-filt.gz

gatk SelectVariants --select-type-to-include INDEL -R {genome} \
-V {raw} -O {outfile}.indel.gz
gatk VariantFiltration -R {genome} \
--filter-expression 'QD < 2.0 || FS > 200.0 || QUAL < 30.0' \
--filter-name 'hard_filter' \
-V {outfile}.indel.gz -O {outfile}.indel.hard-filt.gz
gatk MergeVcfs -I {outfile}.snp.hard-filt.gz -I {outfile}.indel.hard-filt.gz \
-O {outfile}
rm {outfile}.snp* {outfile}.indel* """.format(
            genome = Reference.genome,
            bam=self.bam, raw=self.raw, outfile=self.outfile,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class Annovar(luigi.Task):
    """VCF Annotation using Annovar

    Flow: vcf->avinput1(begin_with_chr1)->avinput2(begin_with_1)->anotated_csv

    Attributes:
        vcf (str): vcf
        outfile (str): default vcf.hg19_multianno.csv

    """
    resources = {"cpu": 2, "memory": 2}
    vcf = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{}.hg19_multianno.final.csv".format(self.vcf))

    def run(self):
        cmd = """
perl {annovar_software_dir}/convert2annovar.pl --format vcf4  {vcf} \
--outfile {vcf}.avinput --withzyg
sed -i 's/^chr//g' {vcf}.avinput

perl {annovar_software_dir}/table_annovar.add_spliceai_exome.pl \
-buildver hg19 -protocol refGene,ensGene,cytoBand,phastConsElements100way,tfbsConsSites,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,avsnp150,snp138NonFlagged,popfreq_all_20150413,popfreq_max_20150413,gnomad211_exome,dbnsfp41a,gerp++gt2,gerp++elem,dbscsnv11,clinvar_20210501,ipmch592,fudan75,ClinPred,hgmd,intervar_20180118,spliceai_exome \
-operation g,g,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . \
-csvout -polish -remove -out {vcf} {vcf}.avinput {annovar_database_dir}

xsv fmt -t "\t" {vcf}.hg19_multianno.csv \
| awk 'BEGIN{{FS="\t";OFS="\t"}}{{if (($122=="Benign") || ((($27>0.05) || ($33>0.05) || ($46>0.05)) && ($122="."))){{print $0}}}}' \
| xsv select Chr,Start,End,Ref,Alt,1000G_ALL,ExAC_ALL,AF -d "\t" > {vcf}.hg19_multianno.dropped_sites.csv

xsv join --left 1-5 {vcf}.hg19_multianno.csv \
1-5 {vcf}.hg19_multianno.dropped_sites.csv \
| xsv fmt -t "\t" | awk 'BEGIN{{FS"\t";OFS="\t"}}{{if ($176==""){{print $0}}}}' \
| xsv select -d "\\t" 1-175 | uniq  > {vcf}.hg19_multianno.final.csv.tmp

head -1 {vcf}.hg19_multianno.csv \
| cat - {vcf}.hg19_multianno.final.csv.tmp > {vcf}.hg19_multianno.final.csv

rm {vcf}.hg19_multianno.final.csv.tmp
        """.format(
            annovar_software_dir = Reference.annovar_software_dir,
            annovar_database_dir = Reference.annovar_database_dir,
            vcf = self.vcf
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


@inherits(Input)
class RemoveDuplicates(luigi.Task):
    """Remove Duplicates Before Calling Variants: Mpileup and Freebayes

    Attributes:
        sample (str): sample name
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    """
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()

    @property
    def inbam(self):
        return "{}/mapping/{}.merged.bam".format(self.outdir, self.sample)

    @property
    def outfile(self):
        return "{}/mapping/{}.dedup.bam".format(self.outdir, self.sample)

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self):
        cmd = "samtools rmdup {bam} {outfile}".format(
            bam=self.inbam, outfile=self.outfile,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class Freebayes(luigi.Task):
    """Freebayes Calling Variants

    Attributes:
        inbam (file): input dedup bam file
        outvcf (str): output vcf file

    """
    resources = {"cpu": 1, "memory": 1}
    inbam = luigi.ListParameter()
    outvcf = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outvcf)

    def run(self):
        cmd = "freebayes -f {genome} {bam} > {outfile}".format(
            genome=Reference.genome, bam=self.inbam, outfile=self.outvcf
        )
        # bcftools filter -e 'QUAL < 20' -s LOWQUAL {rawvcf} {outvcf}
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class Mpileup(luigi.Task):
    """Bcftools Mpileup calling variants

    Attributes:
        inbam (file): input dedup bam file
        outvcf (str): output vcf file

    """
    resources = {"cpu": 1, "memory": 1}
    inbam = luigi.ListParameter()
    outvcf = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outvcf)

    def run(self):
        cmd = """bcftools mpileup -f {genome} {bam} \
| bcftools call -mv --ploidy {ploidy} -o {outvcf}""".format(
            bam=self.inbam, vcf=self.outvcf,
            genome=Reference.genome, ploidy=Reference.genome_version
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


@inherits(Input)
class GATKSingle(luigi.Task):
    """GATK Calling Variants

    Preparation: FiltLowQuality -> Mapping -> PostMapping
    Current Workflow: GATKhaplotypeCaller -> HardFilt|CNNFilt

    Attributes:
        sample (str): sample name
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    """
    resources = {"cpu": 1, "memory": 16}
    sample = luigi.Parameter()

    @property
    def outdir_variants(self):
        outdir_variants = "{}/call-variants".format(self.outdir)
        os.makedirs(outdir_variants, exist_ok=True)
        return outdir_variants

    @property
    def vcf(self):
        return "{}/{}.{}.raw.vcf.gz".format(
            self.outdir_variant, self.sample, Reference.genome_version,
        )

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget("{}.hg19_multianno.final.csv".format(self.vcf))

    def run(self):
        bam = "{}/mapping/{}.bqsr.bam".format(self.outdir, self.sample)
        yield HaplotypeCaller(bam=bam, outfile=self.vcf)
        # hard filt
        #yield HardFilt(
        #    bam=bam, raw=raw_vcf,
        #    outfile="{}.hard-filt.vcf.gz".format(self.prefix)
        #)
        yield Annovar(vcf=self.vcf)


@inherits(Input)
class FreebayesSingle(luigi.Task):
    """Freebayes Calling Variants for Single Sample

    Attributes:
        sample (str): sample name
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    """
    @property
    def outdir_variants(self):
        outdir_variants = "{}/call-variants".format(self.outdir)
        os.makedirs(outdir_variants, exist_ok=True)
        return outdir_variants

    @property
    def vcf(self):
        return "{}/{}.{}.freebayes.raw.vcf.gz".format(
            self.outdir_variants, self.sample, Reference.genome_version,
        )

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget("{}.hg19_multianno.final.csv".format(self.vcf))

    def run(self):
        inbam = "{}/mapping/{}.dedup.bam".format(self.outdir, self.sample)
        yield Freebayes(inbam=inbam, outvcf=self.vcf)
        yield Annovar(vcf=self.vcf)


@inherits(Input)
class MpileupSingle(luigi.Task):
    """Mpileup Calling Variants for Single Sample

    Attributes:
        sample (str): sample name
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    """
    @property
    def outdir_variants(self):
        outdir_variants = "{}/call-variants".format(self.outdir)
        os.makedirs(outdir_variants, exist_ok=True)
        return outdir_variants

    @property
    def vcf(self):
        return "{}/{}.{}.mpileup.raw.vcf.gz".format(
            self.outdir_variants, self.sample, Reference.genome_version,
        )

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget("{}.hg19_multianno.final.csv".format(self.vcf))

    def run(self):
        inbam = "{}/mapping/{}.dedup.bam".format(self.outdir, self.sample)
        yield Mpileup(inbam=inbam, outvcf=self.vcf)
        yield Annovar(vcf=self.vcf)


@inherits(Input)
class DoWES(luigi.WrapperTask):
    """WES Analysis

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    Todo:
        * freebayes
        * mpileup

    """

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")


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
