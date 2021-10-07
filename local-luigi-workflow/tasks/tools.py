import luigi

class Tools(luigi.Config):
    """Software Settings

    Attributes:
        bamdst (str): mapping statistics
        snpeff (str): annotation
    """
    bamdst = luigi.Parameter()
    snpeff = luigi.Parameter(default="hg19")
