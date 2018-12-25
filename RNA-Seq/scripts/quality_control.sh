#!/usr/bin/env bash

set -e

function print_usage {
    cat <<EOF
usage:
$ export trimmo_jar="/path/to/trimmomatic.jar"
$ export adapter_fa="/path/to/adapter.fa"
$ quality_control RAW_DIR CLEAN_DIR QC_REPORT_DIR
        [-t/--threads THREADS]
        [--pbs]

INPUT: raw data: fastq.gz file.

OUTPUT: 1. Fastqc report of raw data.
        2. Trimmed fastqc.gz file.
        3. Fastqc on trimmed data.

NOTE: 1. fastq should named as _R1.fastq.gz _R2.fastq.gz
      2. This pipeline only for pair-end library
      3. log stored in qc_trim_qc.log
EOF
}


## PARSE ARG

if [ -z $trimmo_jar ]; then
    echo "must define variable trimmo_jar firstly, use command:"
    echo '$ export trimmo_jar="/path/to/trimmomatic.jar"'
    exit 1
fi
if [ -z $adapter_fa ]; then
    echo "must define variable adapter_fa firstly, use command:"
    echo '$ export adapter_fa="/path/to/adapter.fa"'
    exit 1
fi

# default arguments
threads=4
pbs=0 # `1` for use pbs system

if ! options=$(getopt -o t: --long threads:,pbs -- "$@"); then
    exit 1
fi

set -- $options

while true ; do
    case "$1" in
        -t|--threads)
            case "$2" in
                "") shift 2 ;;
                *) threads=$2 ; shift 2 ;;
            esac ;;
        --pbs) pbs=1 ; shift ;;
        --) shift ; break ;;
        *) print_usage ; exit 1 ;;
    esac
done

if [ $# -ne 3 ]; then
    print_usage
    exit 1
fi

raw_dir=$1
clean_dir=$2
report_dir=$3

# print parameters
echo "start pipeline with parameters:"
echo "    RAW_DIR: $raw_dir"
echo "    CLEAN_DIR: $clean_dir"
echo "    QC_REPORT_DIR: $report_dir"
echo "    threads: $threads"
echo "    pbs: $pbs"

export raw_dir=${raw_dir//\'/}
export clean_dir=${clean_dir//\'/}
export report_dir=${report_dir//\'/}
export threads=${threads//\'/}
export pbs=${pbs//\'/}

## END PARSE ARG


if [ ! -d $raw_dir ]; then
    echo "[error] raw data dir not exist."
    exit 1
fi

if [ ! -d $report_dir/raw ]; then
    mkdir -p $report_dir/raw
fi
if [ ! -d $report_dir/clean ]; then
    mkdir -p $report_dir/clean
fi

if [ ! -d $clean_dir ]; then
    mkdir $clean_dir
fi


function qc_trim_qc {
    # run related commands on single sample

    # usage: qc ID

    id=$1

    r_r1=$raw_dir/"$id"_R1.fastq.gz
    r_r2=$raw_dir/"$id"_R2.fastq.gz

    # fastqc on raw data
    fastqc -t $threads $r_r1 -o $report_dir/raw
    fastqc -t $threads $r_r2 -o $report_dir/raw

    c_r1=$clean_dir/"$id"_R1.fastq.gz
    c_r2=$clean_dir/"$id"_R2.fastq.gz

    # trimm
    java -jar $trimmo_jar PE $r_r1 $r_r2 \
        $c_r1 $clean_dir/"$id"_R1.unpair.fastq.gz \
        $c_r2 $clean_dir/"$id"_R2.unpair.fastq.gz \
        ILLUMINACLIP:$adapter_fa:3:30:10:1:TRUE LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        -threads $threads

    rm $clean_dir/"$id"_R1.unpair.fastq.gz
    rm $clean_dir/"$id"_R2.unpair.fastq.gz

    # fastqc on clean data
    fastqc -t $threads $c_r1 -o $report_dir/clean
    fastqc -t $threads $c_r2 -o $report_dir/clean
}

export -f qc_trim_qc

function main {
    # deploy all tasks of all samples
    # 

    for id in `ls $raw_dir | grep fastq.gz | sed 's/.fastq.gz//g' | sed 's/_R1//' | sed 's/_R2//' | sort -u`;
    do
        if [ $pbs -eq 1 ]; then
            echo "run on pbs cluster"
            echo "qc_trim_qc $id" | qsub -d $PWD -V -l nodes=1:ppn="$threads" -N QC_"$id"
        else
            echo "run on front end" 
            qc_trim_qc $id 2>> qc_trim_qc.log
        fi
    done

}

main
