#!/usr/bin/env bash
set -euo pipefail

# Calculate callable sites and the de novo SNP mutation rate.
# Denominator = 5 * sum(callable bases shared by queen + 5 offspring)

PROJECT="/localdata/qinti/Project/Bee"
REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
PRIMARY_BED="${PROJECT}/00_rawdata/reference/primary_chromosomes.bed"
MAPPING_SUMMARY="${PROJECT}/02_mapping/stats/mapping_summary.tsv"
BAM_DIR="${PROJECT}/02_mapping/bam"
CONFIRMED="${PROJECT}/06_mutations/snp_candidates/all_families.confirmed_denovo_SNP_candidates.tsv.gz"
OUT="${PROJECT}/06_mutations/callable_sites"

FAMILIES=(C1 C2 C4)
MIN_DP=10
MIN_BQ=10
MIN_MQ=20
MAX_DP_MULTIPLIER=2.5
N_OFFSPRING=5

# Override on the command line if needed:
# N_MUTATIONS=10 DETECTION_SENSITIVITY=0.998 bash script.sh
N_MUTATIONS="${N_MUTATIONS:-10}"
DETECTION_SENSITIVITY="${DETECTION_SENSITIVITY:-1.0}"

mkdir -p "${OUT}/family_masks" "${OUT}/logs"

for tool in samtools awk gzip; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: ${tool} is not in PATH" >&2
        exit 1
    }
done

for file in "${REF}" "${REF}.fai" "${PRIMARY_BED}" "${MAPPING_SUMMARY}"; do
    [[ -s "${file}" ]] || {
        echo "ERROR: missing file: ${file}" >&2
        exit 1
    }
done

awk -v x="${DETECTION_SENSITIVITY}" 'BEGIN {exit !(x > 0 && x <= 1)}' || {
    echo "ERROR: DETECTION_SENSITIVITY must be >0 and <=1" >&2
    exit 1
}

# ------------------------------------------------------------
# 1. BED containing only unambiguous A/C/G/T reference bases
# ------------------------------------------------------------

CONTIG_LIST="${OUT}/primary_chromosomes.txt"
ACGT_BED="${OUT}/primary_chromosomes.ACGT_only.bed"

cut -f1 "${PRIMARY_BED}" | awk 'NF && !seen[$1]++ {print $1}' > "${CONTIG_LIST}"
mapfile -t contigs < "${CONTIG_LIST}"

echo "Creating A/C/G/T-only primary chromosome BED..."
samtools faidx "${REF}" "${contigs[@]}" |
awk 'BEGIN {OFS="\t"}
function flush() {
    if (active) {
        print chrom, start, pos
        active=0
    }
}
/^>/ {
    flush()
    header=substr($0,2)
    split(header,a,/[ \t]/)
    chrom=a[1]
    pos=0
    next
}
{
    seq=toupper($0)
    for (i=1; i<=length(seq); i++) {
        base=substr(seq,i,1)
        if (base ~ /^[ACGT]$/) {
            if (!active) {
                start=pos
                active=1
            }
        } else {
            flush()
        }
        pos++
    }
}
END {flush()}' > "${ACGT_BED}"

PRIMARY_ACGT_BASES=$(
    awk '{n += $3-$2} END {printf "%.0f",n}' "${ACGT_BED}"
)
echo "Primary A/C/G/T bases: ${PRIMARY_ACGT_BASES}"

# ------------------------------------------------------------
# 2. Sample-specific maximum depth:
#    max(30, ceil(mean_depth * 2.5)); missing depth -> 100
# ------------------------------------------------------------

MAX_DP_TABLE="${OUT}/sample_max_depth_thresholds.tsv"
awk -F '\t' -v m="${MAX_DP_MULTIPLIER}" '
BEGIN {OFS="\t"}
NR==1 {
    for (i=1;i<=NF;i++) {
        if ($i=="sample_id") s=i
        if ($i=="mean_depth") d=i
    }
    if (!s || !d) exit 2
    print "sample_id","mean_depth","max_dp"
    next
}
{
    mean=$d
    if (mean=="" || mean=="NA" || mean !~ /^[0-9]+([.][0-9]+)?$/) {
        maxdp=100
    } else {
        raw=mean*m
        maxdp=int(raw)
        if (maxdp<raw) maxdp++
        if (maxdp<30) maxdp=30
    }
    print $s,mean,maxdp
}' "${MAPPING_SUMMARY}" > "${MAX_DP_TABLE}"

declare -A MAX_DP
while IFS=$'\t' read -r sample mean maxdp; do
    [[ "${sample}" == "sample_id" ]] && continue
    MAX_DP["${sample}"]="${maxdp}"
done < "${MAX_DP_TABLE}"

# ------------------------------------------------------------
# 3. Family intersection callable masks
# ------------------------------------------------------------

FAMILY_SUMMARY="${OUT}/family_callable_summary.tsv"
printf "family\tcallable_bases_all_6\tcallable_Mb_all_6\tcallable_fraction_of_primary_ACGT\toffspring_site_opportunities\n" > "${FAMILY_SUMMARY}"

TOTAL_OPPORTUNITIES=0

for family in "${FAMILIES[@]}"; do
    samples=(
        "${family}_Q"
        "${family}_1"
        "${family}_2"
        "${family}_3"
        "${family}_4"
        "${family}_5"
    )
    bams=()
    maxima=()

    for sample in "${samples[@]}"; do
        bam="${BAM_DIR}/${sample}.markdup.bam"
        [[ -s "${bam}" ]] || {
            echo "ERROR: missing BAM: ${bam}" >&2
            exit 1
        }
        if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
            samtools index "${bam}"
        fi
        [[ -n "${MAX_DP[${sample}]:-}" ]] || {
            echo "ERROR: no mean-depth threshold for ${sample}" >&2
            exit 1
        }
        bams+=("${bam}")
        maxima+=("${MAX_DP[${sample}]}")
    done

    mask="${OUT}/family_masks/${family}.queen_plus_5_offspring.callable.bed"
    count_tmp="${OUT}/.${family}.count.tmp"

    echo "Calculating ${family}: ${samples[*]}"
    echo "Maximum DP: ${maxima[*]}"

    # -q is minimum base quality; -Q is minimum mapping quality.
    # -s avoids counting overlapping mates twice.
    # samtools depth excludes duplicate, secondary, QC-failed and
    # unmapped reads by default.
    samtools depth -aa -s -q "${MIN_BQ}" -Q "${MIN_MQ}" -b "${ACGT_BED}" "${bams[@]}" 2> "${OUT}/logs/${family}.samtools_depth.log" |
    awk -v min="${MIN_DP}" -v m1="${maxima[0]}" -v m2="${maxima[1]}" -v m3="${maxima[2]}" -v m4="${maxima[3]}" -v m5="${maxima[4]}" -v m6="${maxima[5]}" -v count_file="${count_tmp}" '
BEGIN {OFS="\t"}
function flush() {
    if (active) {
        print rchr,rstart,rend
        total += rend-rstart
        active=0
    }
}
{
    pass=($3>=min && $3<=m1 &&
          $4>=min && $4<=m2 &&
          $5>=min && $5<=m3 &&
          $6>=min && $6<=m4 &&
          $7>=min && $7<=m5 &&
          $8>=min && $8<=m6)
    if (pass) {
        s=$2-1
        e=$2
        if (active && $1==rchr && s==rend) {
            rend=e
        } else {
            flush()
            rchr=$1
            rstart=s
            rend=e
            active=1
        }
    } else {
        flush()
    }
}
END {
    flush()
    printf "%.0f\n",total > count_file
}' > "${mask}"

    callable=$(<"${count_tmp}")
    rm -f "${count_tmp}"

    callable_mb=$(awk -v n="${callable}" 'BEGIN {printf "%.6f",n/1e6}')
    callable_fraction=$(
        awk -v n="${callable}" -v total="${PRIMARY_ACGT_BASES}" 'BEGIN {printf "%.8f",n/total}'
    )
    opportunities=$(
        awk -v n="${callable}" -v k="${N_OFFSPRING}" 'BEGIN {printf "%.0f",n*k}'
    )
    TOTAL_OPPORTUNITIES=$(
        awk -v a="${TOTAL_OPPORTUNITIES}" -v b="${opportunities}" 'BEGIN {printf "%.0f",a+b}'
    )

    printf "%s\t%s\t%s\t%s\t%s\n" "${family}" "${callable}" "${callable_mb}" "${callable_fraction}" "${opportunities}" >> "${FAMILY_SUMMARY}"

    echo "${family}: ${callable} callable bases; ${opportunities} offspring-site opportunities"
done

# ------------------------------------------------------------
# 4. Validate mutation count and calculate rate
# ------------------------------------------------------------

if [[ -s "${CONFIRMED}" ]]; then
    table_count=$(
        gzip -dc "${CONFIRMED}" |
        awk 'NR>1 && NF {n++} END {print n+0}'
    )
    if [[ "${table_count}" -ne "${N_MUTATIONS}" ]]; then
        echo "ERROR: confirmed table has ${table_count} rows, but N_MUTATIONS=${N_MUTATIONS}" >&2
        exit 1
    fi
else
    echo "WARNING: confirmed table not found; using N_MUTATIONS=${N_MUTATIONS}"
fi

RAW_RATE=$(
    awk -v n="${N_MUTATIONS}" -v d="${TOTAL_OPPORTUNITIES}" 'BEGIN {printf "%.12g",n/d}'
)
EXPOSURE=$(
    awk -v d="${TOTAL_OPPORTUNITIES}" -v s="${DETECTION_SENSITIVITY}" 'BEGIN {printf "%.12f",d*s}'
)
CORRECTED_RATE=$(
    awk -v n="${N_MUTATIONS}" -v d="${TOTAL_OPPORTUNITIES}" -v s="${DETECTION_SENSITIVITY}" 'BEGIN {printf "%.12g",n/(d*s)}'
)

CI_LOW="NA"
CI_HIGH="NA"
if command -v Rscript >/dev/null 2>&1; then
    read -r CI_LOW CI_HIGH < <(
        Rscript -e '
        a <- commandArgs(TRUE)
        z <- poisson.test(as.numeric(a[1]), T=as.numeric(a[2]))
        cat(format(z$conf.int[1],scientific=TRUE,digits=12),
            format(z$conf.int[2],scientific=TRUE,digits=12))
        ' "${N_MUTATIONS}" "${EXPOSURE}"
    )
fi

RATE_SUMMARY="${OUT}/mutation_rate_summary.tsv"
printf "n_confirmed_SNPs\ttotal_callable_offspring_site_opportunities\tdetection_sensitivity\traw_rate_per_site_per_generation\tcorrected_rate_per_site_per_generation\tpoisson_exact_95CI_lower\tpoisson_exact_95CI_upper\n" > "${RATE_SUMMARY}"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${N_MUTATIONS}" "${TOTAL_OPPORTUNITIES}" "${DETECTION_SENSITIVITY}" "${RAW_RATE}" "${CORRECTED_RATE}" "${CI_LOW}" "${CI_HIGH}" >> "${RATE_SUMMARY}"

PARAMETERS="${OUT}/callable_site_parameters.tsv"
{
    printf "parameter\tvalue\n"
    printf "families\tC1,C2,C4\n"
    printf "minimum_depth\t%s\n" "${MIN_DP}"
    printf "minimum_base_quality\t%s\n" "${MIN_BQ}"
    printf "minimum_mapping_quality\t%s\n" "${MIN_MQ}"
    printf "maximum_depth_rule\tmax(30,ceil(mean_depth*2.5))\n"
    printf "reference_bases\tprimary_chromosomes_ACGT_only\n"
    printf "family_callable_rule\tqueen_and_all_5_offspring_callable\n"
    printf "denominator\t5_times_sum_of_family_callable_bases\n"
    printf "GQ_AD_AF_SB\tnot_available_at_invariant_sites; assessed by detection simulation\n"
} > "${PARAMETERS}"

echo
echo "Completed."
echo "Family summary: ${FAMILY_SUMMARY}"
echo "Rate summary:   ${RATE_SUMMARY}"
column -t -s $'\t' "${RATE_SUMMARY}" 2>/dev/null || cat "${RATE_SUMMARY}"
