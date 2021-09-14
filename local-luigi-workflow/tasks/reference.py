import luigi

class Reference(luigi.Config):
    """Reference for Whole Exome Sequencing Analysis. """
    genome = luigi.Parameter()
    genome_dict = luigi.Parameter()
    ploidy = luigi.Parameter()
    interval = luigi.Parameter()
    bed = luigi.Parameter()
    dbsnp = luigi.Parameter()
    mills_gold_standard = luigi.Parameter()
    omni = luigi.Parameter()
    hapmap = luigi.Parameter()
    snp_1000g = luigi.Parameter()
