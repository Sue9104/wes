import luigi
from luigi.util import inherits
import pandas as pd
import numpy as np
import subprocess
import os
import glob
import logging
from datetime import datetime, timedelta
from .reference import Reference

# logging
rootlogger = logging.getLogger("root")
cmdlogger = logging.getLogger("cmd")
timelogger = logging.getLogger("ptime")

##################################################
## Callbacks for execution
##################################################
@luigi.Task.event_handler(luigi.Event.START)
def start_time(task):
    text = 'Task:{task}\tStartTime:{time}'.format(
        task = task.__class__.__name__,
        time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    )
    rootlogger.info(text)

@luigi.Task.event_handler(luigi.Event.SUCCESS)
def end_time(task):
    text = 'Task:{task}\tEndTime:{time}'.format(
        task = task.__class__.__name__,
        time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    )
    rootlogger.info(text)

@luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
def execution_time(task, processing_time):
    end = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start = (datetime.now() - timedelta(seconds = processing_time)).strftime('%Y-%m-%d %H:%M:%S')
    text = "Task:{task}\t{ptime:.1f}min\tStarted:{start}\tEnded:{end}".format(
        task = task.__class__.__name__,
        ptime = processing_time / 60,
        start = start,
        end = end
    )
    rootlogger.info(text)
    timelogger.info(text)


##################################################
## Class Definition
##################################################
class Input(luigi.Config):
    """Input Parameter for analysis

    Attributes:
        infile (str): input csv file, header must be "sample,lane,r1,r2"
        outdir (str, optional): output result directory, default is "wes-output"

    Properties:
        samples (list): a list of sample name from infile
        info (dict): a detail dict for sample information, for example:
            {"sample_name": named_turple("LANE", "R1", "R2")}

    """
    infile = luigi.Parameter()
    outdir = luigi.Parameter()

    @property
    def samples(self):
        os.makedirs(self.outdir, exist_ok=True)
        sample_info_data = pd.read_csv(self.infile)
        samples = list( set(sample_info_data["sample"].values) )
        return samples

    @property
    def info(self):
        info = {}
        sample_info_data = pd.read_csv(self.infile)
        for row in sample_info_data.itertuples(index=False, name="Info"):
            info[row.sample] = \
                info[row.sample] + [row] if row.sample in info else [row]
        return info


class QualityControl(Input, luigi.WrapperTask):
#class QualityControl(Input, luigi.Task):
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

    #def output(self):
    #    return luigi.LocalTarget('{}/multiqc_report.html'.format(self.outdir_qc))

    #def run(self):
    #    cmd = "multiqc -o {0} {0}".format(self.outdir_qc)
    #    rootlogger.info(cmd)
    #    cmdlogger.info(cmd)
    #    subprocess.run(cmd, shell=True)


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

    """

    resources = {"cpu": 1}
    fastq = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        out_prefix = os.path.basename(self.fastq).split('.')[0]
        return luigi.LocalTarget(
            '{}/{}_fastqc.html'.format(self.outdir, out_prefix)
        )

    def run(self):
        cmd = "fastqc -f fastq -o {} {}".format(self.outdir, self.fastq)
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)



class FiltLowQuality(Input, luigi.WrapperTask):
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
                R1=data.r1, R2=data.r2,
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
            adapter = Reference().illumina_adapter, outdir = self.outdir,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)

class Mapping(Input, luigi.Task):
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
    resources = {"cpu": 10, "memory": 1}
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
            "{}/quality-control/{}/coverage.report".format(
                self.outdir, self.sample,
            )
        )

    def run(self):
        # mapping by lanes
        sample_bams = yield [
            BWA(
                sample=data.sample, lane=data.lane, outdir=self.outdir_mapping,
                fq1="{outdir}/trim-adapter/{sample}_{lane}_1P.fq.gz".format(
                    outdir = self.outdir, sample=data.sample, lane=data.lane,
                ),
                fq2="{outdir}/trim-adapter/{sample}_{lane}_2P.fq.gz".format(
                    outdir = self.outdir, sample=data.sample, lane=data.lane,
                ),
            )
            for data in self.info[self.sample]
        ]
        # merge multiple bams
        merged_bam = "{}/{}.merged.bam".format(self.outdir_mapping, self.sample)
        yield MergeSampleBams(
            outbam=merged_bam,
            inbams = [bam.path for bam in sample_bams],
        )
        # mapping statistics
        yield MappingStatistics(
            inbam = merged_bam,
            outdir = "{}/quality-control/{}".format(self.outdir, self.sample),
        )


class PostMapping(Input, luigi.Task):
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

    def run(self):
        merged_bam = "{}/{}.merged.bam".format(self.outdir_mapping, self.sample)
        markdup_bam ="{}/{}.markdup.bam".format(self.outdir_mapping, self.sample)
        # mark PCR dupliates
        yield MarkDuplicate(inbam=merged_bam, outbam=markdup_bam)
        # recalibrate base quality
        yield BaseScoreQualityRecalibrator(
            inbam=markdup_bam, outbam=self.output().path,
        )


class MergeSampleBams(luigi.Task):
    """Merge Multiple Bam Files for One Sample and Coverage Statistics

    Attributes:
        inbams (list): a list of bam files
        outbam (str): output bam filename

    Output:
        - {outdir}/mapping/{sample}.merged.bam
    """
    resources = {"cpu": 2, "memory": 1}
    outbam = luigi.Parameter()
    inbams = luigi.ListParameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outbam)

    def run(self):
        cmd = """samtools merge - {bams} | tee {outfile} \
| samtools index - {outfile}.bai """.format(
            bams = " ".join(self.inbams), outfile = self.outbam,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)


class MappingStatistics(luigi.Task):
    """Mapping Coverage and Depth Statistics

    Attributes:
        inbam (str): bam file
        outdir (str): output directory

    Output:
        - {outdir}/{sample}/coverage.report
    """
    resources = {"cpu": 2, "memory": 1}
    outdir = luigi.Parameter()
    inbam = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("{}/coverage.report".format(self.outdir))

    def run(self):
        os.makedirs(self.outdir, exist_ok=True)
        cmd = "bamdst -p {bed} -o {outdir}/ {bam}".format(
            outdir=self.outdir, bam=self.inbam,
            bed=Reference().bed,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)


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
    resources = {"cpu": 30, "memory": 10}

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
        cmd = """bwa mem -t 60 \
-R "@RG\\tID:{lane}\\tSM:{sample}\\tLB:WES\\tPL:Illumina" {genome} {fq1} {fq2} \
| samtools sort -@ 2 -o {outdir}/{sample}_{lane}.sorted.bam - """.format(
            sample=self.sample, lane=self.lane, fq1=self.fq1, fq2=self.fq2,
            outdir=self.outdir, genome=Reference().genome,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)


class MarkDuplicate(luigi.Task):
    """MarkDupliate after Mapping

    Attributes:
        inbam (str): input bam file
        outbam (str): output bam name

    Output:
        - {outdir}/{sample}.mark_dedup.bam
    """

    resources = {"cpu": 20, "memory": 256}
    inbam = luigi.Parameter()
    outbam = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outbam)

    def run(self):
        cmd ="""gatk --java-options '-Xmx256G' MarkDuplicatesSpark \
-conf 'spark.executor.cores=40' --input {bam} \
--output {outfile} --metrics-file {outfile}.metrics.txt""".format(
            outfile=self.outbam, bam=self.inbam,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)

class BaseScoreQualityRecalibrator(luigi.Task):
    """BaseScoreQualityRecalibrator after MarkDuplicate

    Attributes:
        inbam (str): bam file
        outbam (str): output bam

    """
    resources = {"cpu": 4, "memory": 256}
    inbam = luigi.Parameter()
    outbam = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outbam)

    def run(self):
        cmd = """gatk BaseRecalibrator \
--java-options '-Xmx256G -XX:GCTimeLimit=50 \
-XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal' \
--reference {genome} --intervals {interval} \
--known-sites {mills_gold_standard} --known-sites {snp_1000g} \
--known-sites {dbsnp} --known-sites {indel_1000g} --known-sites {hapmap} \
--known-sites {omni} \
--input {bam} --output {outfile}.recal_table

gatk ApplyBQSR --java-options '-Xmx256g \
-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
-XX:+PrintGCDetails -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10' \
--reference {genome} --input {bam} \
--bqsr-recal-file {outfile}.recal_table --output {outfile}""".format(
            genome=Reference().genome, interval=Reference().interval,
            mills_gold_standard = Reference().mills_gold_standard,
            snp_1000g=Reference().snp_1000g, indel_1000g=Reference().indel_1000g,
            hapmap=Reference().hapmap, omni=Reference().omni, dbsnp=Reference().dbsnp,
            bam=self.inbam, outfile=self.outbam,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
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
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)


class HaplotypeCaller(luigi.Task):
    """GATK HaplotypeCaller

    CNNFilt and HardFilt is available. Default no filter is used.

    Attributes:
        bam (str): bam file
        outvcf (str): output vcf name

    """

    resources = {"cpu": 10, "memory": 16}
    bam = luigi.Parameter()
    outvcf = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outvcf)

    def run(self):
        cmd = """gatk --java-options '-Xmx256g' HaplotypeCaller \
-G StandardAnnotation -G StandardHCAnnotation \
--pair-hmm-implementation AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads 10 \
-R {genome} --dbsnp {dbsnp} \
-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
-I {bam} -O {outfile} -bamout {outfile}.bam

gatk CollectHsMetrics -R {genome} \
--BAIT_INTERVALS {interval} --TARGET_INTERVALS {interval} \
-I {outfile} -O {outfile}.hs-metrics.txt """.format(
            genome=Reference().genome, version=Reference().genome_version,
            dbsnp=Reference().dbsnp, interval=Reference().interval,
            bam=self.bam, outfile=self.outvcf,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)

class CNNFilt(luigi.Task):
    """CNNFilt for Raw VCF

    Attributes:
        bam (str): bam file
        rawvcf (str): raw vcf
        outvcf (str): output vcf name

    """

    resources = {"cpu": 4, "memory": 16}
    bam = luigi.Parameter()
    rawvcf = luigi.Parameter()
    outvcf = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outvcf)

    def run(self):
        cmd = """gatk --java-options '-Xmx16G' CNNScoreVariants \
-tensor-type read_tensor -R {genome} \
-I {bam} -V {raw} -O {outfile}.cnn_mark.gz

gatk FilterVariantTranches --info-key CNN_2D \
--resource {hapmap} --resource {mills_gold_standard} \
--snp-tranche 99.95 --indel-tranche 99.4 --invalidate-previous-filters \
-V {outfile}.cnn_mark.gz -O {outfile}""".format(
            genome = Reference().genome, version = Reference().genome_version,
            hapmap = Reference().hapmap,
            mills_gold_standard = Reference().mills_gold_standard,
            bam = self.bam, raw = self.rawvcf, outfile = self.outvcf,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)

class HardFilt(luigi.Task):
    """HardFilt for Raw VCF

    Attributes:
        bam (str): bam file
        rawvcf (str): raw vcf
        outvcf (str): output vcf name

    """
    resources = {"cpu": 1, "memory": 1}
    bam = luigi.Parameter()
    rawvcf = luigi.Parameter()
    outvcf = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(self.outvcf)

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
            genome = Reference().genome,
            bam=self.bam, raw=self.rawvcf, outfile=self.outvcf,
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
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
            annovar_software_dir = Reference().annovar_software_dir,
            annovar_database_dir = Reference().annovar_database_dir,
            vcf = self.vcf
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)


class RemoveDuplicates(Input, luigi.Task):
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
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
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
            genome=Reference().genome, bam=self.inbam, outfile=self.outvcf
        )
        # bcftools filter -e 'QUAL < 20' -s LOWQUAL {rawvcf} {outvcf}
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
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
            genome=Reference().genome, ploidy=Reference().genome_version
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell=True)


class GATKSingle(Input, luigi.Task):
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

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget(
            "{}/backup/{}/backup.finished".format(self.outdir, self.sample)
        )

    def run(self):
        vcf = "{}/{}.{}.gatk.raw.vcf.gz".format(
            self.outdir_variants, self.sample, Reference().genome_version,
        )
        yield HaplotypeCaller(bam=self.input().path, outvcf=vcf)
        # hard filt
        #yield HardFilt(
        #    bam=bam, rawvcf=raw_vcf,
        #    outvcf="{}.hard-filt.vcf.gz".format(self.prefix)
        #)
        #yield Annovar(vcf=vcf)
        yield Backup(sample=self.sample, indir=self.outdir)


class FreebayesSingle(Input, luigi.Task):
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

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget(
            "{}/backup/{}/backup.finished".format(self.outdir, self.sample)
        )

    def run(self):
        vcf = "{}/{}.{}.gatk.raw.vcf.gz".format(
            self.outdir_variants, self.sample, Reference().genome_version,
        )
        yield Freebayes(inbam=self.input().path, outvcf=vcf)
        #yield Annovar(vcf=vcf)
        yield Backup(sample=self.sample, indir=self.outdir)


class MpileupSingle(Input, luigi.Task):
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

    def requires(self):
        raise NotImplementedError("Task Need to implement requires function")

    def output(self):
        return luigi.LocalTarget(
            "{}/backup/{}/backup.finished".format(self.outdir, self.sample)
        )

    def run(self):
        vcf = "{}/{}.{}.mpileup.raw.vcf.gz".format(
            self.outdir_variants, self.sample, Reference().genome_version,
        )
        yield Mpileup(inbam=self.input().path, outvcf=vcf)
        #yield Annovar(vcf=vcf)
        yield Backup(sample=self.sample, indir=self.outdir)


class DoWES(Input, luigi.WrapperTask):
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


class Backup(luigi.Task):
    """Backup for Important Files

    Attributes:
        indir: analysis directory
        sample: sample name

    Output:
        - bam
        - mapping statistics
        - vcf
        - annotated csv

    """
    resources = {"cpu": 1}
    sample = luigi.Parameter()
    indir = luigi.Parameter()

    @property
    def outdir_backup(self):
        outdir_backup = "{}/backup/{}".format(self.indir, self.sample)
        os.makedirs(outdir_backup, exist_ok=True)
        return outdir_backup

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget('{}/backup.finished'.format(self.outdir_backup))

    def run(self):
        vcfs = glob.glob("{indir}/call-variants/{sample}*.raw.vcf.gz".format(
            indir = self.indir, sample=self.sample,
        ))
        csvs = glob.glob(
            "{indir}/call-variants/{sample}*hg19_multianno*csv".format(
                indir = self.indir, sample=self.sample)
        )
        cmd = """cp -r {indir}/quality-control/{sample} {backup}/quality-control\
 && cp {vcf} {csv} {backup}/ """.format(
            indir=self.indir, sample=self.sample,
            vcf = " ".join(vcfs), csv=" ".join(csvs),
            backup=self.outdir_backup,
        )
        # for gatk
        bqsr_files = glob.glob(
            "{}/mapping/{}*bqsr.ba*".format(self.indir,self.sample)
        )
        if len(bqsr_files) > 0:
            cmd += " && cp {} {}/".format(" ".join(bqsr_files), self.outdir_backup)
        # for
        dedup_files = glob.glob(
            "{}/mapping/{}*dedup.ba*".format(self.indir,self.sample)
        )
        if len(dedup_files) > 0:
            cmd += " && cp {} {}/".format(" ".join(dedup_files), self.outdir_backup)
        cmd += " && touch {}/backup.finished".format(self.outdir_backup)
        rootlogger.info(cmd)
        subprocess.run(cmd, shell=True)


class Qualimap(luigi.Task):
    """Evaluate alignment data: bamqc"""
    resources = {"cpu": 10, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget(
            "{outdir}/mapping/{sample}/qualimapReport.html".format(
                outdir = self.outdir, sample = self.sample
            )
        )

    def run(self):
        cmd = "qualimap bamqc -bam {outdir}/mapping/{sample}.sorted.bam -nt 20 \
            -outdir {outdir}/mapping/{sample} -outformat PDF:HTML".format(
                outdir = self.outdir, sample = self.sample
        )
        rootlogger.info(cmd)
        cmdlogger.info(cmd)
        subprocess.run(cmd, shell = True)

