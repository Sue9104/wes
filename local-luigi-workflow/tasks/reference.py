import luigi

class Reference(luigi.Config):
    """Genome Reference for Exome Sequencing Analysis

    Using gatk resource bundle: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle.
    It can be downloaded in Google cloud from gs://gatk-legacy-bundles.
    Please make sure following attributes are in the luigi.cfg of your running path.

    Attributes:
        genome (str): human hg19 genome fasta path
        genome_version (str, optional): default is hg19
        interval(str): interval list file using by gatk, suffix must be ".interval_list"
        dbsnp(str): latest dbsnp vcf file (hg19)
        mills_gold_standard(str): Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
        snp_1000g(str): 1000G_phase1.snps.high_confidence.hg19.sites.vcf
        phase3_1000g(str): Not available for hg19, default is None
        indel_1000g(str): 1000G_phase1.indels.hg19.sites.vcf
        omni(str): 1000G_omni2.5.hg19.sites.vcf
        hapmap(str): hapmap_3.3.hg19.sites.vcf
    """
    genome = luigi.Parameter()
    genome_version = luigi.Parameter(default="hg19")
    interval = luigi.Parameter()
    bed = luigi.Parameter()
    dbsnp = luigi.Parameter()
    mills_gold_standard = luigi.Parameter()
    omni = luigi.Parameter()
    hapmap = luigi.Parameter()
    snp_1000g = luigi.Parameter()
    phase3_1000g = luigi.Parameter(default=None)
    indel_1000g = luigi.Parameter()

