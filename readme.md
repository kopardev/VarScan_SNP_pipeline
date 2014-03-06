VarScan_SNP_pipeline
********************
Dependent on:
1. [Samtools](http://samtools.sourceforge.net)
2. [VarScan](http://varscan.sourceforge.net/index.html)
3. [Parallel](http://www.gnu.org/software/parallel)

Make SNP calls using mpileup followed by VarScan
INPUT ... txt file containing list of bam files, generally one per sample
OUTPUT ... vcf file
