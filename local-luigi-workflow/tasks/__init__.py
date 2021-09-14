from .do_wes import DoWES, Annovar
#from .annotation import Annovar, SnpEff
from .reference import Reference
from .quality_control import PrepareFastq, MultiQC, FastQC, Trimmomatic, Qualimap
from .mapping import BWA
from .data_process_before_calling import RemoveDuplicates, MarkDuplicate, BaseScoreQualityRecalibrator, ApplyBQSR
from .call_variants import Freebayes, Mpileup, HaplotypeCaller, HardFilter
