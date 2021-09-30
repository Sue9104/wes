from tasks import *
import luigi

class DoWES(luigi.WrapperTask):
    infile = luigi.Parameter()
    outdir = luigi.Parameter()
    def requires(self):
        return [
            QualityControl(infile=self.infile, outdir=self.outdir),
            CallVariants(infile=self.infile, outdir=self.outdir),
            Annotation(infile=self.infile, outdir=self.outdir),
            Backup(infile=self.infile, outdir=self.outdir)
        ]

class QualityControl(tasks.QualityControl):
    def requires(self):
        fastqs = []
        for _, lanes in self.info:
            for lane in lanes:
                fastqs += [lane.r1, lane.r2]
        return [FastQC(fastq=fastq, outdir=self.outdir_qc)for fastq in fasts]

class FastQC(tasks.FastQC):
    def requires(self):
        return CheckExists(file=self.fastq)
