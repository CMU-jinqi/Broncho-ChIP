#!/bin/bash

# ChIP-seq双端数据分析完整流程

set -euo pipefail

### 初始化变量 ###
CONTROL_SRR=""
TREATMENT_SRR=""
GENOME=""
OUTDIR=""
RAWDATA="./rawdata"
MACS3_OPTIONS=""
THREADS=4

### 参数解析模块 ###
show_help() {
    echo "Usage: $0 [options]"
    echo "必需参数："
    echo "  --control SRR1,SRR2   对照组SRR列表（逗号分隔）"
    echo "  --treatment SRR3,SRR4 处理组SRR列表（逗号分隔）"
    echo "  --genome GENOME.fa    参考基因组FASTA文件"
    echo "  --outdir DIR          输出目录路径"
    echo "可选参数："
    echo "  --rawdata DIR         原始数据根目录（默认：./rawdata）"
    echo "  --macs3option 'OPT'   MACS3额外参数（需用引号包裹）"
    echo "  --threads N           并行线程数（默认：4）"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --control)      CONTROL_SRR="$2"; shift 2 ;;
        --treatment)    TREATMENT_SRR="$2"; shift 2 ;;
        --genome)       GENOME="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --rawdata)      RAWDATA="$2"; shift 2 ;;
        --macs3option)  MACS3_OPTIONS="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        -h|--help)      show_help ;;
        *) echo "错误：未知参数 $1"; exit 1 ;;
    esac
done

### 核心函数模块 ###
check_dependency() {
    local missing=()
    for cmd in prefetch fasterq-dump bowtie2 samtools macs3 bedtools bamCoverage; do
        if ! command -v $cmd &>/dev/null; then
            missing+=("$cmd")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "错误：缺少必需命令 - ${missing[*]}"
        exit 1
    fi
}

find_sra_file() {
    local sra="$1"
    # 扩展搜索路径，兼容新旧版prefetch结构
    local search_paths=(
        "$RAWDATA/sra/${sra}.sra"                 # 传统路径
        "$RAWDATA/sra/${sra}/${sra}.sra"          # 新版嵌套结构
        "/mnt/Ddisk/project/Chipseq-pipe/rawdata_sra/${sra}.sra"
        "$OUTDIR/sra/${sra}.sra"                  # 当前流程旧版结构
        "$OUTDIR/sra/${sra}/${sra}.sra"           # 当前流程新版结构
    )
    
    for path in "${search_paths[@]}"; do
        if [[ -f "$path" ]]; then
            echo "找到SRA文件：$path"
            echo "$path"
            return 0
        fi
    done
    return 1
}

process_sra_download() {
    local sra="$1"
    local target_dir="$OUTDIR/sra"
    mkdir -p "$target_dir"
    
    echo "下载SRA数据：$sra"
    if ! prefetch --output-directory "$target_dir" --type sra "$sra"; then
        echo "错误：SRA下载失败"
        # 清理可能存在的空目录
        rm -rf "${target_dir}/${sra}"
        exit 1
    fi

    # 处理新版嵌套结构
    local nested_path="${target_dir}/${sra}/${sra}.sra"
    if [[ -f "$nested_path" ]]; then
        echo "检测到嵌套结构，移动文件..."
        mv -v "$nested_path" "$target_dir/"
        rmdir "${target_dir}/${sra}" 2>/dev/null || true
    fi

    # 最终路径验证
    local sra_path="${target_dir}/${sra}.sra"
    if [[ ! -f "$sra_path" ]]; then
        echo "错误：SRA文件未出现在预期路径"
        exit 1
    fi
    echo "SRA文件已就绪：$sra_path"
}

convert_fastq() {
    local sra="$1"
    echo "转换FASTQ文件：$sra"
    
    mkdir -p "$OUTDIR/fastq"
    (
        cd "$OUTDIR/fastq"
        if ! fasterq-dump --split-3 --threads $THREADS "../sra/${sra}.sra"; then
            echo "错误：FASTQ转换失败"
            exit 1
        fi
        
        # 验证文件完整性
        if [[ ! -f "${sra}_1.fastq" || ! -f "${sra}_2.fastq" ]]; then
            echo "错误：生成不完整的FASTQ文件"
            exit 1
        fi
        
        # 并行压缩
        if ! (pigz -p $THREADS "${sra}_1.fastq" && pigz -p $THREADS "${sra}_2.fastq"); then
            echo "错误：压缩失败"
            exit 1
        fi
        
        # 清理中间文件
        rm -f "${sra}_1.fastq" "${sra}_2.fastq"
    ) || exit 1
}

process_sra() {
    local sra="$1"
    echo "===== 处理样本 $sra ====="
    
    # 检查 FASTQ 文件是否已存在
    if [[ -f "$OUTDIR/fastq/${sra}_1.fastq.gz" && -f "$OUTDIR/fastq/${sra}_2.fastq.gz" ]]; then
        echo "FASTQ文件已存在，跳过处理"
        return
    fi
    
    # 获取SRA文件
    if ! src_sra=$(find_sra_file "$sra"); then
        process_sra_download "$sra"
    else
        echo "使用现有SRA文件：$src_sra"
        mkdir -p "$OUTDIR/sra"
        cp -v "$src_sra" "$OUTDIR/sra/${sra}.sra"
    fi
    
    # 转换FASTQ
    convert_fastq "$sra"
}

align_reads() {
    local sra="$1"
    local index_prefix="$OUTDIR/index/$(basename "${GENOME%.*}")"
    
    # 检查 BAM 文件是否已存在
    if [[ -f "$OUTDIR/bam/${sra}.sorted.bam" && -f "$OUTDIR/bam/${sra}.sorted.bam.bai" ]]; then
        echo "BAM文件已存在，跳过比对"
        return
    fi
    
    echo "比对样本：$sra"
    bowtie2 -p $THREADS -x "$index_prefix" \
        -1 "$OUTDIR/fastq/${sra}_1.fastq.gz" \
        -2 "$OUTDIR/fastq/${sra}_2.fastq.gz" \
        -S "$OUTDIR/bam/${sra}.sam" || {
            echo "错误：比对失败"
            exit 1
        }
    
    # 处理BAM文件
    samtools view -@ $THREADS -Sb "$OUTDIR/bam/${sra}.sam" > "$OUTDIR/bam/${sra}.bam"
    samtools sort -@ $THREADS -m 2G -o "$OUTDIR/bam/${sra}.sorted.bam" "$OUTDIR/bam/${sra}.bam"
    samtools index -@ $THREADS "$OUTDIR/bam/${sra}.sorted.bam"
    rm -f "$OUTDIR/bam/${sra}.sam" "$OUTDIR/bam/${sra}.bam"
}

calculate_effective_genome_size() {
    echo "计算有效基因组大小..."
    if [[ -f "$OUTDIR/effective_genome_size.txt" ]]; then
        echo "有效基因组大小已存在，跳过计算"
        effective_genome_size=$(cat "$OUTDIR/effective_genome_size.txt")
        return
    fi

    # 使用 samtools faidx 统计有效碱基数
    samtools faidx "$GENOME"
    awk '{sum += $2} END {print sum}' "$GENOME.fai" > "$OUTDIR/effective_genome_size.txt"
    effective_genome_size=$(cat "$OUTDIR/effective_genome_size.txt")
    echo "有效基因组大小：$effective_genome_size"
}

### 主流程 ###
main() {
    # 初始化检查
    check_dependency
    [[ -z "$CONTROL_SRR" || -z "$TREATMENT_SRR" ]] && { 
        echo "错误：必须指定--control和--treatment参数"
        show_help
        exit 1
    }
    mkdir -p "$OUTDIR"/{fastq,sra,bam,index,merged_bam,peaks,bigwig}
    
    # 处理所有样本
    for sra in ${CONTROL_SRR//,/ } ${TREATMENT_SRR//,/ }; do
        process_sra "$sra"
    done
    
    # 构建索引
    local index_prefix="$OUTDIR/index/$(basename "${GENOME%.*}")"
    if [[ ! -f "${index_prefix}.1.bt2" ]]; then
        echo "构建Bowtie2索引..."
        bowtie2-build --threads $THREADS "$GENOME" "$index_prefix" || {
            echo "错误：索引构建失败"
            exit 1
        }
    fi
    
    # 比对步骤
    for sra in ${CONTROL_SRR//,/ } ${TREATMENT_SRR//,/ }; do
        align_reads "$sra"
    done
    
    # 合并BAM文件
    echo "合并对照组..."
    control_bams=($(printf "$OUTDIR/bam/%s.sorted.bam " ${CONTROL_SRR//,/ }))
    if [[ ! -f "$OUTDIR/merged_bam/control_merged.bam" ]]; then
        samtools merge -@ $THREADS -o "$OUTDIR/merged_bam/control_merged.bam" "${control_bams[@]}"
    else
        echo "对照组 BAM 文件已存在，跳过合并"
    fi
    
    echo "合并处理组..."
    treatment_bams=($(printf "$OUTDIR/bam/%s.sorted.bam " ${TREATMENT_SRR//,/ }))
    if [[ ! -f "$OUTDIR/merged_bam/treatment_merged.bam" ]]; then
        samtools merge -@ $THREADS -o "$OUTDIR/merged_bam/treatment_merged.bam" "${treatment_bams[@]}"
    else
        echo "处理组 BAM 文件已存在，跳过合并"
    fi
    
    # 峰检测
    echo "运行MACS3..."
    if [[ ! -f "$OUTDIR/peaks/chipseq_peaks.narrowPeak" ]]; then
        macs3 callpeak \
            -t "$OUTDIR/merged_bam/treatment_merged.bam" \
            -c "$OUTDIR/merged_bam/control_merged.bam" \
            -n chipseq \
            --outdir "$OUTDIR/peaks" \
            $MACS3_OPTIONS || {
                echo "错误：MACS3运行失败"
                exit 1
            }
    else
        echo "峰值文件已存在，跳过 MACS3"
    fi
    
    # 计算有效基因组大小
    calculate_effective_genome_size
    
    # 使用 deeptools 将 sorted bam 文件转换为 bigwig 文件
    echo "将 sorted bam 文件转换为 bigwig 文件..."
    for bam_file in "$OUTDIR/bam"/*.sorted.bam; do
        base_name=$(basename "$bam_file" .sorted.bam)
        if [[ ! -f "$OUTDIR/bigwig/${base_name}.bw" ]]; then
            bamCoverage -b "$bam_file" -o "$OUTDIR/bigwig/${base_name}.bw" \
                --binSize 10 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize "$effective_genome_size" \
                --extendReads 200
        else
            echo "BigWig 文件已存在，跳过：$base_name"
        fi
    done
    
    echo "===== 分析成功完成 ====="
    echo "输出目录：$OUTDIR"
    echo "峰值文件：ls -l $OUTDIR/peaks/chip*"
    echo "BigWig 文件：ls -l $OUTDIR/bigwig/*.bw"
}

### 执行入口 ###
main "$@"
