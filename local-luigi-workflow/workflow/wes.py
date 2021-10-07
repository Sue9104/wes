from tasks import *
import luigi


class QualityControl(tasks.QualityControl):
    def requires(self):
        fastqs = []
        for _, lanes in self.info:
            for lane in lanes:
                fastqs += [lane.r1, lane.r2]
        return [FastQC(fastq=fastq, outdir=self.outdir_qc)for fastq in fastqs]

class FastQC(tasks.FastQC):
    def requires(self):
        return CheckExists(file=self.fastq)

class Mapping(tasks.Mapping):
    def requires(self):
        return FiltLowQuality(
            infile=self.infile, outdir=self.outdir, sample=self.sample,
        )

class PostMapping(tasks.PostMapping):
    def requires(self):
        return Mapping(
            infile=self.infile, outdir=self.outdir, sample=self.sample,
        )

class RemoveDuplicates(tasks.RemoveDuplicates):
    def requires(self):
        return Mapping(
            infile=self.infile, outdir=self.outdir, sample=self.sample,
        )

class GATKSingle(tasks.GATKSingle):
    def requires(self):
        return PostMapping(
            infile=self.infile, outdir=self.outdir, sample=self.sample,
        )

class FreebayesSingle(tasks.FreebayesSingle):
    def requires(self):
        return RemoveDuplicates(
            infile=self.infile, outdir=self.outdir, sample=self.sample,
        )

class MpileupSingle(tasks.MpileupSingle):
    def requires(self):
        return RemoveDuplicates(
            infile=self.infile, outdir=self.outdir, sample=self.sample,
        )

class DoWES(tasks.DoWES):
    def requires(self):
        tasks = [QualityControl(infile=self.infile, outdir=self.outdir)]
        tasks += [
            GATKSingle(
                infile=self.infile, outdir=self.outdir, sample=sample
            ) for sample in self.samples
        ]

