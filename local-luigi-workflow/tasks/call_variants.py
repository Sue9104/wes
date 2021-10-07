import luigi
import os
import subprocess
import logging
from .reference import Reference


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

