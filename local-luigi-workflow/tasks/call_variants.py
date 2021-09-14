import luigi
import os
import subprocess
import logging
from .reference import Reference

class HaplotypeCaller(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()
    gvcf = luigi.BoolParameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/call-variants/{sample}.gatk.raw.vcf".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        cmd = "gatk --java-options '-Xmx4g' HaplotypeCaller \
        -G StandardAnnotation -G StandardHCAnnotation \
        -R {genome} -L {interval} \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        -I {outdir}/mapping/{sample}.bqsr.bam \
        -O {outdir}/call-variants/{sample}.gatk.raw.vcf \
        -bamout {outdir}/call-variants/{sample}.gatk.raw.bam".format(
            genome = Reference().genome, interval = Reference().interval,
            sample = self.sample, outdir = self.outdir
        )
        if self.gvcf:
            cmd += " -ERC GVCF"
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class CNNFilter(luigi.Task):
    """Filter using deep learning, experimental"""
    resources = {"cpu": 1, "memory": 1}
    sample_bam = luigi.Parameter()
    sample_vcf = luigi.Parameter()
    output_vcf = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget(self.output_vcf)

    def run(self):
        cmd = """
gatk CNNScoreVariants -tensor-type read_tensor -R {genome} \
    -I {input_bam} -V {input_vcf} -O {temp_vcf}
gatk FilterVariantTranches --info-key CNN_2D \
    --resource {hapmap} --resource {mills_gold_standard} \
    --snp-tranche 99.95 --indel-tranche 99.4 --invalidate-previous-filters \
    -V {temp_vcf} -O {output_vcf}
        """.format(
            genome = Reference().genome, hapmap = Reference().hapmap,
            mills_gold_standard = Reference().mills_gold_standard,
            input_bam = self.sample_bam, input_vcf = self.sample_vcf,
            temp_vcf = self.sample_vcf.replace('vcf', 'cnn.vcf'),
            output_vcf = self.output_vcf
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


class HardFilter(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/call-variants/{sample}.gatk.final.vcf".format(
            outdir = self.outdir, sample = self.sample
        ))

    def run(self):
        cmd = """
gatk SelectVariants -R {genome} -select-type SNP \
    -V {outdir}/call-variants/{sample}.gatk.raw.vcf \
    -O {outdir}/call-variants/{sample}.gatk.raw.snp.vcf
gatk --java-options '-Xmx4g' VariantFiltration -R {genome} \
    -filter 'QD < 2.0' --filter-name 'QD2' \
    -filter 'FS > 60.0' --filter-name 'FS60' \
    -filter 'MQ < 40.0' --filter-name 'MQ40' \
    -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
    -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8.0' \
    -V {outdir}/call-variants/{sample}.gatk.raw.snp.vcf \
    -O {outdir}/call-variants/{sample}.gatk.final.snp.vcf
gatk SelectVariants -R {genome} -select-type INDEL \
    -V {outdir}/call-variants/{sample}.gatk.raw.vcf \
    -O {outdir}/call-variants/{sample}.gatk.raw.indel.vcf
gatk --java-options '-Xmx4g' VariantFiltration -R {genome} \
    -filter 'QD < 2.0' --filter-name 'QD2' \
    -filter 'FS > 200.0' --filter-name 'FS200' \
    -filter 'SOR > 10.0' --filter-name 'SOR10' \
    -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20.0' \
    -V {outdir}/call-variants/{sample}.gatk.raw.indel.vcf \
    -O {outdir}/call-variants/{sample}.gatk.final.indel.vcf
gatk MergeVcfs -I {outdir}/call-variants/{sample}.gatk.raw.snp.vcf \
    -I {outdir}/call-variants/{sample}.gatk.raw.indel.vcf \
    -O {outdir}/call-variants/{sample}.gatk.final.vcf
rm {outdir}/call-variants/{sample}.gatk.raw.*.vcf
    """.format(
            genome = Reference().genome,
            outdir = self.outdir, sample = self.sample
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class CombineGVCF(luigi.Task):
    resources = {"cpu": 4, "memory": 1}
    vcfs = luigi.ListParameter()
    output_vcf = luigi.ListParameter()
    db = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget(self.output)

    def run(self):
        joined_vcf = " -V ".join(self.vcfs)
        cmd = """
gatk --java-options '-Xmx4g -Xms4g' GenomicsDBImport \
    --reader-threads 4 --batch-size 50 -L {interval} \
    {joined_vcf} --genomicsdb-workspace-path {db}
gatk --java-options '-Xmx4g' GenotypeGVCFs -R {genome} \
    -D {dbsnp} -G StandardAnnotation -L {interval} -V gendb://{db} \
    -O {output_vcf}
        """.format(
            interval = Reference().interval, genome = Reference().genome,
            dbsnp = Reference().dbsnp, db = self.db,
            joined_vcf = joined_vcf, output_vcf = self.output_vcf
        )
        logging.info(cmd)
        subprocess.run(cmd)


class VQSR(luigi.Task):
    """Require: 1 whole genome or 30 exomes"""
    resources = {"cpu": 1, "memory": 1}

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget()

    def run(self):
        cmd = """
gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
    --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet \
    -V {cohort.vcf.gz} -O {cohort_excesshet.vcf.gz}  \
gatk MakeSitesOnlyVcf -I cohort_excesshet.vcf.gz -O cohort_sitesonly.vcf.gz
gatk --java-options '-Xmx10g -Xms10g' VariantRecalibrator --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL \
    --max-gaussians 4 \
    -resource mills,known=false,training=true,truth=true,prior=12:Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource axiomPoly,known=false,training=true,truth=false,prior=10:Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    -resource dbsnp,known=true,training=false,truth=false,prior=2:Homo_sapiens_assembly38.dbsnp138.vcf \
    -V cohort_sitesonly.vcf.gz -O cohort_indels.recal \
    --tranches-file cohort_indels.tranches
gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP \
    --max-gaussians 6 \
    -resource hapmap,known=false,training=true,truth=true,prior=15:hapmap_3.3.hg38.vcf.gz \
    -resource omni,known=false,training=true,truth=true,prior=12:1000G_omni2.5.hg38.vcf.gz \
    -resource 1000G,known=false,training=true,truth=false,prior=10:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource dbsnp,known=true,training=false,truth=false,prior=7:Homo_sapiens_assembly38.dbsnp138.vcf \
    -V cohort_sitesonly.vcf.gz -O cohort_snps.recal \
    --tranches-file cohort_snps.tranches
gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
    --recal-file cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -V cohort_excesshet.vcf.gz -O indel.recalibrated.vcf.gz
gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
    --recal-file ${snps_recalibration} \
    --tranches-file ${snps_tranches} \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode SNP \
    -V indel.recalibrated.vcf.gz -O snp.recalibrated.vcf.gz \
        """.format(
            genome = Reference().genome, dbsnp = Reference().dbsnp,
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class CollectHsMetrics(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample_bam = luigi.Parameter()
    output_metric = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget(self.output_metric)

    def run(self):
        cmd = "gatk CollectHsMetrics -R {genome} -I {input_bam} -O {output} \
               --BAIT_INTERVALS {interval} --TARGET_INTERVALS {interval}".format(
                   genome = Reference().genome, interval = Reference().interval,
                   input_bam = self.sample_bam, output = self.output_metric
               )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class CollectVariantCallingMetrics(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    sample_bam = luigi.Parameter()
    output_metric = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget()

    def run(self):
        cmd = "gatk --java-options '-Xmx2g' CollectVariantCallingMetrics \
               --SEQUENCE_DICTIONARY {genome_dict} --TARGET_INTERVALS {interval} \
               --DBSNP {dbsnp} -I {input_vcf} -O {output} ".format(
                   genome_dict = Reference().genome_dict, dbsnp = Reference().dbsnp,
                   interval = Reference().interval,
                   input_vcf = self.sample_vcf, output = self.output_metric
               )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class Freebayes(luigi.Task):
    """Freebayes calling variants

    Flow: Trimmomatic -> Mapping -> RemoveDuplicates -> CallVariants

    Attributes:
        sample (list): list of sample id
        outdir (str): the parent directory of output

    Todo:
        * remove the time limitation after all mapping is done.
    """
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/call-variants/{sample}.freebayes.raw.vcf".format(
            outdir = self.outdir, sample = self.sample))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        cmd = """
freebayes -f {genome} {outdir}/mapping/{sample}.deduped.bam > {outdir}/call-variants/{sample}.freebayes.raw.vcf
bcftools filter -e 'QUAL < 20' -s LOWQUAL {outdir}/call-variants/{sample}.freebayes.raw.vcf > {outdir}/call-variants/{sample}.freebayes.final.vcf
        """.format(
            genome=Reference().genome, sample = self.sample, outdir=self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)

class Mpileup(luigi.Task):
    """Bcftools Mpileup calling variants

    Flow: Trimmomatic -> Mapping -> RemoveDuplicates -> CallVariants

    Attributes:
        sample (list): list of sample id
        outdir (str): the parent directory of output

    Todo:
        * remove the time limitation after all mapping is done.
    """
    resources = {"cpu": 1, "memory": 1}
    sample = luigi.Parameter()
    outdir = luigi.Parameter()

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/call-variants/{sample}.mpileup.raw.vcf".format(
            outdir = self.outdir, sample = self.sample))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'call-variants'), exist_ok=True)
        cmd = """
bcftools mpileup -f {genome} {outdir}/mapping/{sample}.deduped.bam \
    | bcftools call -mv --ploidy {ploidy} -o {outdir}/call-variants/{sample}.mpileup.raw.vcf
        """.format(
            genome = Reference().genome, ploidy = Reference().ploidy,
            sample = self.sample, outdir = self.outdir
        )
        logging.info(cmd)
        subprocess.run(cmd, shell=True)


