# ChIP-seq Analysis Pipeline

## Usage

```bash
bash chipseq_pipeline.sh [options]
必需参数：
--control SRR1,SRR2：对照组 SRR 列表（逗号分隔）

--treatment SRR3,SRR4：处理组 SRR 列表（逗号分隔）

--genome GENOME.fa：参考基因组 FASTA 文件

--outdir DIR：输出目录路径

可选参数：
--rawdata DIR：原始数据根目录（默认：./rawdata）

--macs3option 'OPT'：MACS3 额外参数（需用引号包裹）

--threads N：并行线程数（默认：4）

Output
*.bam：sorted BAM 和 BAI 文件

*.bigwig：bigWig 可视化文件

*.fastq：fastq 格式的测序文件

genome index：Bowtie2 生成的索引

merged_bam：合并后的 BAM 文件

peaks：MACS3 输出的 peaks

*.sra：原始测序文件

effective_genome_size.txt：MACS3 所需的基因组大小文件

Data Preparation
确保参考基因组（FASTA 文件）已下载并保存至本地。

Installation
# Step 1: 创建环境
conda env create -f chipseq_v2.yml

# Step 2: 激活环境
conda activate chipseq-pipe2

Example
bash chipseq_pipeline.sh \
  --control SRR23063365,SRR23063366 \
  --treatment SRR23063364,SRR23063363 \
  --genome /mnt/Ddisk/project/genome/genome_mouse_mm10/mm10.fa \
  --outdir H327K_out
Workflow
使用 SRAtools 下载 .sra 文件并转换为 .fastq。

使用 Bowtie2 对参考基因组进行索引。

将 .fastq 文件比对至参考基因组，生成 .sam 文件。

将 .sam 文件转换为 .bam 文件并进行排序。

使用 MACS3 对排序后的 BAM 文件进行 peak calling。

将 BAM 文件转换为 .bigwig 文件，用于 IGV 可视化。

