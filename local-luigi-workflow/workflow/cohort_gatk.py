import luigi

class GATK(luigi.Task):
    """GATK Calling Variants

    ClassFlow: Trimmomatic -> Mapping -> MarkDuplicate -> BaseScoreQualityRecalibrator -> ApplyBQSR -> HaplotypeCaller -> HardFilter or CNN
    Process: Remove adapter and filt low quality bases -> Mapping to genome -> Mark duplicate -> Recalibrate base score quality -> Call variants

    Attributes:
        samples (list): List of sample id
        outdir (int, optional): Parent directory of output

    Todo:
        * Compare hard filter and CNN
    """
    resources = {"cpu": 1, "memory": 1}
    samples = luigi.ListParameter()
