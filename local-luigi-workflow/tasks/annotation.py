import luigi
import subprocess
import logging
import os

class Annovar(luigi.Task):
    resources = {"cpu": 1, "memory": 1}
    annovar = luigi.Parameter()
    annovar_db = luigi.Parameter()
    annovar_db_version = luigi.Parameter()

    sample = luigi.Parameter()
    outdir = luigi.Parameter()
    input_vcf = luigi.Parameter(significant = False)
    output_prefix = luigi.Parameter(significant = False)

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{prefix}.{version}_multianno.vcf".format(
            prefix = self.output_prefix, version = self.annovar_db_version))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'annotation'), exist_ok=True)
        cmd = "perl {annovar} {input_vcf} -out {prefix} \
            {annovar_db} -buildver {annovar_db_version} \
            -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp35c,gnomad211_exome,esp6500siv2_all,clinvar_20190305 \
            -operation g,r,f,f,f,f,f,f \
            -remove -nastring . -vcfinput -polish".format(
                input_vcf = self.input_vcf, prefix = self.output_prefix,
                annovar = self.annovar, annovar_db = self.annovar_db, annovar_db_version = self.annovar_db_version
        )
        logging.info(cmd)
        subprocess.run(cmd, shell = True)

class SnpEff(luigi.Task):
    resources = {"cpu": 1, "memory": 8}
    snpeff = luigi.Parameter()
    snpeff_db_version = luigi.Parameter()

    sample = luigi.Parameter()
    outdir = luigi.Parameter()
    input_vcf = luigi.Parameter(significant = False)
    output_prefix = luigi.Parameter(significant = False)

    def requires(self):
        raise NotImplemetedError("Need to be implemented!")

    def output(self):
        return luigi.LocalTarget("{outdir}/annotation/{sample}.{prefix}.snpEff.summary.csv".format(
            outdir = self.outdir, sample = self.sample, prefix = self.output_prefix
        ))

    def run(self):
        os.makedirs(os.path.join(self.outdir, 'annotation'), exist_ok=True)
        cmd = """
java -Xmx8g -jar {snpeff} {version} {input_vcf} \
    -csvStats {outdir}/annotation/{sample}.{output_prefix}.snpEff.summary.csv \
    > {outdir}/annotation/{sample}.{output_prefix}.snpEff.vcf
        """.format(
            snpeff = self.snpeff, version = self.snpeff_db_version,
            outdir = self.outdir, sample = self.sample,
            input_vcf = self.input_vcf, output_prefix = self.output_prefix
        )
        logging.info(cmd)
        subprocess.run(cmd, shell = True)

