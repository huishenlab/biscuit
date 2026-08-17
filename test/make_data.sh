# FIXME: Doesn't quite setup data as expected - will need to create the `data/` directory and then run everything in there

# Sherman rarely needs to be rerun, but it's left here on the off-chance that it is needed (set `run_sherman="true"`)
run_sherman="false"

if [ "${run_sherman}" = "true" ]; then
    # Check Sherman is in PATH
    if [[ `which Sherman 2>&1 > /dev/null` ]]; then
        >&2 echo "Cannot find Sherman in PATH. Add Sherman to PATH and rerun."
        exit 1
    else
        echo Using Sherman found at `which Sherman` | tee make_data.log
    fi

    echo Simulating Paired-end WGBS data
    Sherman \
        --length 150 \
        --number_of_seqs 1000 \
        --genome_folder ref \
        --paired_end \
        --CG_conversion 30 \
        --CH_conversion 99 > sim_pe_wgbs_output.txt 2>&1

    mv simulated_1.fastq sim_pe_wgbs_1.fastq
    mv simulated_2.fastq sim_pe_wgbs_2.fastq
    gzip -f sim_pe_wgbs_1.fastq
    gzip -f sim_pe_wgbs_2.fastq

    echo Simulating Single-end WGBS data
    Sherman \
        --length 150 \
        --number_of_seqs 1000 \
        --genome_folder ref \
        --CG_conversion 30 \
        --CH_conversion 99 > sim_se_wgbs_output.txt 2>&1

    mv simulated.fastq sim_se_wgbs.fastq
    gzip -f sim_se_wgbs.fastq

    echo Simulating Paired-end WGS data
    Sherman \
        --length 150 \
        --number_of_seqs 1000 \
        --genome_folder ref \
        --paired_end \
        --conversion_rate 0 > sim_pe_wgs_output.txt 2>&1

    mv simulated_1.fastq sim_pe_wgs_1.fastq
    mv simulated_2.fastq sim_pe_wgs_2.fastq
    gzip -f sim_pe_wgs_1.fastq
    gzip -f sim_pe_wgs_2.fastq

    echo Simulating Single-end WGS data
    Sherman \
        --length 150 \
        --number_of_seqs 1000 \
        --genome_folder ref \
        --conversion_rate 0 > sim_se_wgs_output.txt 2>&1

    mv simulated.fastq sim_se_wgs.fastq
    gzip -f sim_se_wgs.fastq
else
    echo "Skipping Sherman run and moving on to test data creation" | tee make_data.log
fi

function create_data {
    # Inputs
    log=${1}
    tag=${2}
    ise=${3}
    fqs=${4}

    # Create output directories
    adir=dynamic/${tag}_00_align
    pdir=dynamic/${tag}_01_pileup
    vdir=dynamic/${tag}_02_vcf2bed
    mdir=dynamic/${tag}_03_mergecg
    bdir=dynamic/${tag}_04_bsconv
    sdir=dynamic/${tag}_05_bsstrand
    cdir=dynamic/${tag}_06_cinread
    qdir=dynamic/${tag}_07_qc
    rdir=dynamic/${tag}_08_qc_coverage

    mkdir -p ${adir}
    mkdir -p ${pdir}
    mkdir -p ${vdir}
    mkdir -p ${mdir}
    mkdir -p ${bdir}
    mkdir -p ${sdir}
    mkdir -p ${cdir}
    mkdir -p ${qdir}
    mkdir -p ${rdir}

    # Prep single-end command line options
    single_end=
    if [ ${ise} = "true" ]; then
        single_end="-s"
    fi

    echo "Making ${tag} files" | tee --append ${log}
    echo -e "\tbiscuit align" | tee --append ${log}
    biscuit align ${IDX} ${fqs} > ${adir}/current.sam
    biscuit align -v 4 ${IDX} ${fqs} > ${adir}/current.debug
    samtools sort -o ${adir}/current.bam -O BAM --write-index ${adir}/current.sam

    echo -e "\tbiscuit pileup" | tee --append ${log}
    biscuit pileup -o ${pdir}/current.vcf ${REF} ${adir}/current.bam

    echo -e "\tbiscuit vcf2bed" | tee --append ${log}
    biscuit vcf2bed ${pdir}/current.vcf > ${vdir}/current.bed

    echo -e "\tbiscuit mergecg" | tee --append ${log}
    biscuit mergecg ${REF} ${vdir}/current.bed > ${mdir}/current.mergecg.bed

    echo -e "\tbiscuit bsconv" | tee --append ${log}
    biscuit bsconv ${REF} ${adir}/current.bam > ${bdir}/current.bsconv

    echo -e "\tbiscuit bsstrand" | tee --append ${log}
    biscuit bsstrand ${REF} ${adir}/current.bam ${sdir}/current.bss.bam 2> ${sdir}/current.bss
    samtools view -h -o ${sdir}/current.bss.sam -O SAM ${sdir}/current.bss.bam

    echo -e "\tbiscuit cinread" | tee --append ${log}
    biscuit cinread ${REF} ${adir}/current.bam > ${cdir}/current.cinread

    echo -e "\tbiscuit qc" | tee --append ${log}
    biscuit qc ${single_end} ${REF} ${adir}/current.bam ${qdir}/current

    echo -e "\tbiscuit qc_coverage" | tee --append ${log}
    biscuit qc_coverage \
        -P ${rdir}/current \
        -B ref/chr22.gc.bot10p.bed.gz \
        -T ref/chr22.gc.top10p.bed.gz \
        ${REF} \
        ref/chr22.cpg.bed.gz \
        ${adir}/current.bam
}

echo "Test data made with biscuit version:" | tee --append make_data.log
biscuit version 2>&1 | tee --append make_data.log # Need to run separately as this outputs to stderr

mkdir -p dynamic

# BISCUIT index (done once)
echo "Making biscuit index files" | tee --append make_data.log
IDX=dynamic/000_index/current
REF=ref/chr22.fa.gz
mkdir -p dynamic/000_index
biscuit index -p ${IDX} ${REF}

# WGBS Paired-end files
create_data make_data.log wgbs_pe "false" "sim_pe_wgbs_1.fastq.gz sim_pe_wgbs_2.fastq.gz"

# WGS Paired-end files
create_data make_data.log wgs_pe "false" "sim_pe_wgs_1.fastq.gz sim_pe_wgs_2.fastq.gz"

# WGBS Single-end files
create_data make_data.log wgbs_se "true" "sim_se_wgbs.fastq.gz"

# WGS Single-end files
create_data make_data.log wgs_se "true" "sim_se_wgs.fastq.gz"
