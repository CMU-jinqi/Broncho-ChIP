Usage: bash chipseq_pipeline.sh [options]
必需参数：
  --control SRR1,SRR2   对照组SRR列表（逗号分隔）
  --treatment SRR3,SRR4 处理组SRR列表（逗号分隔）
  --genome GENOME.fa    参考基因组FASTA文件
  --outdir DIR          输出目录路径
可选参数：
  --rawdata DIR         原始数据根目录（默认：./rawdata）
  --macs3option 'OPT'   MACS3额外参数（需用引号包裹）
  --threads N           并行线程数（默认：4）

======================================================================
#output:
--bam  #sorted bam and bai files
--bigwig #bigwig files
--fastq #fastq files
--index #genome index
--merged_bam #merged bam files
--peaks #macs3 output
--sra #sra files
--effective_genome_size.txt #genome size
======================================================================
#Data preparation
The reference genome sequence needs to be downloaded and saved to your local disk.

#1. Installation:
conda env create -f chipseq_v2.yml

#2. Activate environment:
conda activate chipseq-pipe2

#3. Example
bash chipseq_pipeline.sh --control SRR23063365,SRR23063366 --treatment SRR23063364,SRR23063363 --genome /mnt/Ddisk/project/genome/genome_mouse_mm10/mm10.fa  --outdir H327K_out

#4. Workflow
The script will use sratools to download the sra files and convert them to fastq files. The reference genome will be indexed using Bowtie2 and the fastq files are then be aligned to the reference genome to obtain sam files which are converted to bam files and then sorted bam files. The sorted bam files are then used by MACS3 to call peaks. The sorted bam files will be converted to bigwig files used for igv visulization.
