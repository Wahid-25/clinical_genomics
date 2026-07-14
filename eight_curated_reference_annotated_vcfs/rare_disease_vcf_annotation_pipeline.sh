#!/usr/bin/env bash
# rare_disease_vcf_annotation_pipeline.sh
# Purpose: annotate a human WGS VCF for rare/genetic-disorder interpretation.
#
# This is a pipeline TEMPLATE. It wires together common open-source tools and
# local reference databases. You must download/install the tools and databases
# first, then edit annotation_resources.env and run this script.
#
# Major branches:
#   1) SNV/indel annotation: bcftools -> VEP -> SnpEff/SnpSift -> gnomAD/ClinVar/ClinGen -> SpliceAI -> ANNOVAR -> InterVar/ACMG
#   2) CNV/SV interpretation: AnnotSV -> ClassifyCNV -> ISV-CNV -> optional Horizon command placeholder
#
# Not clinical advice. Use for research/training only unless validated in your lab.

set -Eeuo pipefail
IFS=$'\n\t'

usage() {
  cat <<'USAGE'
Usage:
  bash rare_disease_vcf_annotation_pipeline.sh \
    -i sample.small_variants.vcf.gz \
    -o results_dir \
    -c annotation_resources.env \
    [-n sample.cnv.vcf.gz|sample.cnv.bed] \
    [-s SAMPLE_ID] \
    [-a GRCh38|GRCh37] \
    [-t THREADS]

Required:
  -i  SNV/indel VCF, bgzipped VCF, or BCF
  -o  Output directory
  -c  Resource configuration file; see annotation_resources.env.example

Optional:
  -n  CNV/SV input as VCF or BED. For ClassifyCNV/ISV-CNV, BED must be:
      chrom  start  end  DEL|DUP
  -s  Sample name/prefix. Default: derived from input filename
  -a  Assembly: GRCh38 or GRCh37. Default: GRCh38
  -t  Threads. Default: 8

Examples:
  bash rare_disease_vcf_annotation_pipeline.sh \
    -i synthetic_sickle_cell_wgs_unannotated.vcf.gz \
    -o annotated_results \
    -c annotation_resources.env \
    -s SCD_SYNTHETIC \
    -a GRCh38 \
    -t 8

  bash rare_disease_vcf_annotation_pipeline.sh \
    -i proband.snvs_indels.vcf.gz \
    -n proband.cnvs.bed \
    -o proband_annotation \
    -c annotation_resources.env
USAGE
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command in PATH: $1"; }
need_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }
maybe_file() { [[ -n "${1:-}" && -s "$1" ]]; }

INPUT_VCF=""
CNV_INPUT=""
OUTDIR=""
CONFIG=""
SAMPLE=""
ASSEMBLY="GRCh38"
THREADS=8

while getopts ":i:n:o:c:s:a:t:h" opt; do
  case "$opt" in
    i) INPUT_VCF="$OPTARG" ;;
    n) CNV_INPUT="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    c) CONFIG="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    a) ASSEMBLY="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) die "Option -$OPTARG requires an argument" ;;
    \?) die "Unknown option: -$OPTARG" ;;
  esac
done

[[ -n "$INPUT_VCF" ]] || { usage; die "Missing -i input VCF"; }
[[ -n "$OUTDIR" ]] || { usage; die "Missing -o output directory"; }
[[ -n "$CONFIG" ]] || { usage; die "Missing -c config file"; }
need_file "$INPUT_VCF"
need_file "$CONFIG"
source "$CONFIG"

case "$ASSEMBLY" in
  GRCh38|grch38|hg38) ASSEMBLY="GRCh38"; ANNOVAR_BUILDVER="${ANNOVAR_BUILDVER:-hg38}"; CLASSIFY_BUILDVER="${CLASSIFY_BUILDVER:-hg38}"; SPLICEAI_ASSEMBLY="${SPLICEAI_ASSEMBLY:-grch38}" ;;
  GRCh37|grch37|hg19) ASSEMBLY="GRCh37"; ANNOVAR_BUILDVER="${ANNOVAR_BUILDVER:-hg19}"; CLASSIFY_BUILDVER="${CLASSIFY_BUILDVER:-hg19}"; SPLICEAI_ASSEMBLY="${SPLICEAI_ASSEMBLY:-grch37}" ;;
  *) die "Unsupported assembly: $ASSEMBLY. Use GRCh38 or GRCh37." ;;
esac

if [[ -z "$SAMPLE" ]]; then
  base=$(basename "$INPUT_VCF")
  SAMPLE="${base%.vcf.gz}"
  SAMPLE="${SAMPLE%.vcf}"
  SAMPLE="${SAMPLE%.bcf}"
fi

mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/work" "$OUTDIR/snv" "$OUTDIR/cnv" "$OUTDIR/acmg" "$OUTDIR/reports"
LOGFILE="$OUTDIR/logs/${SAMPLE}.pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log "Starting rare disease VCF annotation pipeline"
log "Sample: $SAMPLE"
log "Assembly: $ASSEMBLY"
log "Input VCF: $INPUT_VCF"
log "Output directory: $OUTDIR"

# Required core tools
need_cmd bcftools
need_cmd bgzip
need_cmd tabix
need_cmd awk
need_cmd sed
need_cmd sort
need_cmd java
need_cmd python3
need_file "${REF_FASTA:?Set REF_FASTA in config}"
[[ -s "${REF_FASTA}.fai" ]] || die "Missing FASTA index: ${REF_FASTA}.fai. Run: samtools faidx $REF_FASTA"

# ---------- STEP 1: standardize input VCF ----------
RAW_VCFGZ="$OUTDIR/work/${SAMPLE}.input.vcf.gz"
NORM_VCFGZ="$OUTDIR/work/${SAMPLE}.normalized.split.vcf.gz"

log "Step 1: bgzip/index input and normalize/split multiallelic SNV/indel records"
case "$INPUT_VCF" in
  *.vcf.gz)
    cp -f "$INPUT_VCF" "$RAW_VCFGZ"
    ;;
  *.bcf)
    bcftools view -Oz -o "$RAW_VCFGZ" "$INPUT_VCF"
    ;;
  *.vcf)
    bgzip -c "$INPUT_VCF" > "$RAW_VCFGZ"
    ;;
  *)
    die "Input must be .vcf, .vcf.gz, or .bcf"
    ;;
esac

tabix -f -p vcf "$RAW_VCFGZ"

# This normalizes SNVs/indels. If your VCF contains symbolic SV/CNV records,
# keep those in a separate CNV/SV input via -n because bcftools norm may not be
# appropriate for all symbolic variants.
bcftools norm \
  -f "$REF_FASTA" \
  -m -any \
  -Oz \
  -o "$NORM_VCFGZ" \
  "$RAW_VCFGZ"
tabix -f -p vcf "$NORM_VCFGZ"
log "Normalized VCF: $NORM_VCFGZ"

CURRENT_VCF="$NORM_VCFGZ"

# ---------- STEP 2: VEP annotation + plugins ----------
VEP_VCFGZ="$OUTDIR/snv/${SAMPLE}.vep.vcf.gz"
if [[ "${RUN_VEP:-1}" == "1" ]]; then
  need_cmd vep
  [[ -n "${VEP_CACHE_DIR:-}" ]] || die "Set VEP_CACHE_DIR in config"
  [[ -d "$VEP_CACHE_DIR" ]] || die "VEP_CACHE_DIR does not exist: $VEP_CACHE_DIR"

  log "Step 2: run Ensembl VEP with optional plugins/scores"
  VEP_PLUGIN_ARGS=()
  [[ -n "${VEP_PLUGIN_DIR:-}" ]] && VEP_PLUGIN_ARGS+=(--dir_plugins "$VEP_PLUGIN_DIR")

  if maybe_file "${ALPHAMISSENSE_TSV_GZ:-}"; then
    VEP_PLUGIN_ARGS+=(--plugin "AlphaMissense,file=${ALPHAMISSENSE_TSV_GZ}")
  fi
  if maybe_file "${REVEL_TSV_GZ:-}"; then
    VEP_PLUGIN_ARGS+=(--plugin "REVEL,file=${REVEL_TSV_GZ}")
  fi
  if maybe_file "${CADD_SNV_TSV_GZ:-}" && maybe_file "${CADD_INDEL_TSV_GZ:-}"; then
    VEP_PLUGIN_ARGS+=(--plugin "CADD,snv=${CADD_SNV_TSV_GZ},indels=${CADD_INDEL_TSV_GZ}")
  elif maybe_file "${CADD_SNV_TSV_GZ:-}"; then
    VEP_PLUGIN_ARGS+=(--plugin "CADD,snv=${CADD_SNV_TSV_GZ}")
  fi

  vep \
    --input_file "$CURRENT_VCF" \
    --output_file "$VEP_VCFGZ" \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --species homo_sapiens \
    --assembly "$ASSEMBLY" \
    --cache \
    --offline \
    --dir_cache "$VEP_CACHE_DIR" \
    --fasta "$REF_FASTA" \
    --fork "$THREADS" \
    --everything \
    --symbol \
    --canonical \
    --mane \
    --hgvs \
    --numbers \
    --protein \
    --biotype \
    "${VEP_PLUGIN_ARGS[@]}"
  tabix -f -p vcf "$VEP_VCFGZ"
  CURRENT_VCF="$VEP_VCFGZ"
else
  log "Step 2: skipped VEP because RUN_VEP=0"
fi

# ---------- STEP 3: SnpEff ANN field ----------
SNPEFF_VCFGZ="$OUTDIR/snv/${SAMPLE}.vep.snpeff.vcf.gz"
if [[ "${RUN_SNPEFF:-1}" == "1" ]]; then
  need_file "${SNPEFF_JAR:?Set SNPEFF_JAR in config}"
  [[ -n "${SNPEFF_GENOME:-}" ]] || die "Set SNPEFF_GENOME in config, e.g. GRCh38.99 or GRCh37.75 depending on your SnpEff database"

  log "Step 3: run SnpEff to add ANN consequence field"
  java -Xmx"${JAVA_MEM:-8g}" -jar "$SNPEFF_JAR" ann \
    -v \
    -canon \
    -hgvs \
    "$SNPEFF_GENOME" \
    "$CURRENT_VCF" | bgzip -c > "$SNPEFF_VCFGZ"
  tabix -f -p vcf "$SNPEFF_VCFGZ"
  CURRENT_VCF="$SNPEFF_VCFGZ"
else
  log "Step 3: skipped SnpEff because RUN_SNPEFF=0"
fi

# ---------- STEP 4: ClinVar annotation ----------
CLINVAR_VCFGZ="$OUTDIR/snv/${SAMPLE}.clinvar.vcf.gz"
if [[ "${RUN_CLINVAR:-1}" == "1" && -n "${CLINVAR_VCF_GZ:-}" ]]; then
  need_file "$CLINVAR_VCF_GZ"
  need_file "${SNPSIFT_JAR:?Set SNPSIFT_JAR in config if RUN_CLINVAR=1}"

  log "Step 4: annotate ClinVar fields with SnpSift"
  java -Xmx"${JAVA_MEM:-8g}" -jar "$SNPSIFT_JAR" annotate \
    -id \
    -info "${CLINVAR_INFO_FIELDS:-CLNSIG,CLNREVSTAT,CLNDN,CLNDISDB,CLNHGVS,CLNVC,CLNVCSO,GENEINFO}" \
    "$CLINVAR_VCF_GZ" \
    "$CURRENT_VCF" | bgzip -c > "$CLINVAR_VCFGZ"
  tabix -f -p vcf "$CLINVAR_VCFGZ"
  CURRENT_VCF="$CLINVAR_VCFGZ"
else
  log "Step 4: skipped ClinVar because RUN_CLINVAR=0 or CLINVAR_VCF_GZ not set"
fi

# ---------- STEP 5: gnomAD frequency annotation ----------
GNOMAD_VCFGZ="$OUTDIR/snv/${SAMPLE}.gnomad.vcf.gz"
if [[ "${RUN_GNOMAD:-1}" == "1" && -n "${GNOMAD_VCF_GZ:-}" ]]; then
  need_file "$GNOMAD_VCF_GZ"
  log "Step 5: annotate gnomAD population frequencies with bcftools"

  GNOMAD_HEADER="$OUTDIR/work/gnomad.header.txt"
  cat > "$GNOMAD_HEADER" <<'HDR'
##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency copied from source VCF INFO/AF">
##INFO=<ID=GNOMAD_AC,Number=A,Type=Integer,Description="gnomAD allele count copied from source VCF INFO/AC">
##INFO=<ID=GNOMAD_AN,Number=1,Type=Integer,Description="gnomAD allele number copied from source VCF INFO/AN">
##INFO=<ID=GNOMAD_AF_POPMAX,Number=A,Type=Float,Description="gnomAD popmax allele frequency if present in source VCF">
HDR

  # Edit GNOMAD_ANNOTATE_COLUMNS in config if your gnomAD file uses different tags.
  bcftools annotate \
    -a "$GNOMAD_VCF_GZ" \
    -h "$GNOMAD_HEADER" \
    -c "${GNOMAD_ANNOTATE_COLUMNS:-INFO/GNOMAD_AF:=INFO/AF,INFO/GNOMAD_AC:=INFO/AC,INFO/GNOMAD_AN:=INFO/AN}" \
    -Oz \
    -o "$GNOMAD_VCFGZ" \
    "$CURRENT_VCF"
  tabix -f -p vcf "$GNOMAD_VCFGZ"
  CURRENT_VCF="$GNOMAD_VCFGZ"
else
  log "Step 5: skipped gnomAD because RUN_GNOMAD=0 or GNOMAD_VCF_GZ not set"
fi

# ---------- STEP 6: ClinGen dosage / gene-disease region annotation ----------
CLINGEN_VCFGZ="$OUTDIR/snv/${SAMPLE}.clingen.vcf.gz"
if [[ "${RUN_CLINGEN:-1}" == "1" && -n "${CLINGEN_DOSAGE_BED_GZ:-}" ]]; then
  need_file "$CLINGEN_DOSAGE_BED_GZ"
  log "Step 6: annotate ClinGen dosage sensitivity BED overlaps"

  CLINGEN_HEADER="$OUTDIR/work/clingen.header.txt"
  cat > "$CLINGEN_HEADER" <<'HDR'
##INFO=<ID=CLINGEN_REGION,Number=.,Type=String,Description="ClinGen dosage sensitivity region/gene name from BED overlap">
##INFO=<ID=CLINGEN_HAPLO,Number=.,Type=String,Description="ClinGen haploinsufficiency dosage score from BED overlap">
##INFO=<ID=CLINGEN_TRIPLO,Number=.,Type=String,Description="ClinGen triplosensitivity dosage score from BED overlap">
HDR

  # Expected BED.GZ columns: chrom, start0, end, CLINGEN_REGION, CLINGEN_HAPLO, CLINGEN_TRIPLO
  # Must be bgzip-compressed and tabix-indexed with: tabix -p bed file.bed.gz
  bcftools annotate \
    -a "$CLINGEN_DOSAGE_BED_GZ" \
    -h "$CLINGEN_HEADER" \
    -c CHROM,FROM,TO,INFO/CLINGEN_REGION,INFO/CLINGEN_HAPLO,INFO/CLINGEN_TRIPLO \
    -Oz \
    -o "$CLINGEN_VCFGZ" \
    "$CURRENT_VCF"
  tabix -f -p vcf "$CLINGEN_VCFGZ"
  CURRENT_VCF="$CLINGEN_VCFGZ"
else
  log "Step 6: skipped ClinGen because RUN_CLINGEN=0 or CLINGEN_DOSAGE_BED_GZ not set"
fi

# ---------- STEP 7: SpliceAI ----------
SPLICEAI_VCFGZ="$OUTDIR/snv/${SAMPLE}.spliceai.vcf.gz"
if [[ "${RUN_SPLICEAI:-1}" == "1" ]]; then
  if [[ "${SPLICEAI_MODE:-standalone}" == "standalone" ]]; then
    need_cmd spliceai
    log "Step 7: run standalone SpliceAI"
    spliceai \
      -I "$CURRENT_VCF" \
      -O "$OUTDIR/work/${SAMPLE}.spliceai.raw.vcf" \
      -R "$REF_FASTA" \
      -A "$SPLICEAI_ASSEMBLY"
    bgzip -c "$OUTDIR/work/${SAMPLE}.spliceai.raw.vcf" > "$SPLICEAI_VCFGZ"
    tabix -f -p vcf "$SPLICEAI_VCFGZ"
    CURRENT_VCF="$SPLICEAI_VCFGZ"
  elif [[ "${SPLICEAI_MODE:-}" == "precomputed" && -n "${SPLICEAI_SCORES_VCF_GZ:-}" ]]; then
    need_file "$SPLICEAI_SCORES_VCF_GZ"
    log "Step 7: transfer precomputed SpliceAI INFO/SpliceAI from indexed VCF"
    bcftools annotate \
      -a "$SPLICEAI_SCORES_VCF_GZ" \
      -c INFO/SpliceAI \
      -Oz \
      -o "$SPLICEAI_VCFGZ" \
      "$CURRENT_VCF"
    tabix -f -p vcf "$SPLICEAI_VCFGZ"
    CURRENT_VCF="$SPLICEAI_VCFGZ"
  else
    log "Step 7: skipped SpliceAI because SPLICEAI_MODE is not configured"
  fi
else
  log "Step 7: skipped SpliceAI because RUN_SPLICEAI=0"
fi

# ---------- STEP 8: optional bcftools transfer of MAVERICK score ----------
MAVERICK_VCFGZ="$OUTDIR/snv/${SAMPLE}.maverick.vcf.gz"
if [[ "${RUN_MAVERICK:-1}" == "1" && -n "${MAVERICK_TSV_GZ:-}" ]]; then
  need_file "$MAVERICK_TSV_GZ"
  log "Step 8: annotate MAVERICK score from tabix-indexed TSV"

  MAVERICK_HEADER="$OUTDIR/work/maverick.header.txt"
  cat > "$MAVERICK_HEADER" <<'HDR'
##INFO=<ID=MAVERICK_SCORE,Number=A,Type=Float,Description="MAVERICK pathogenicity score transferred from local TSV">
HDR
  # Expected TSV.GZ columns by default: CHROM POS REF ALT MAVERICK_SCORE
  # Must be tabix-indexed with sequence/position coordinates.
  bcftools annotate \
    -a "$MAVERICK_TSV_GZ" \
    -h "$MAVERICK_HEADER" \
    -c "${MAVERICK_COLUMNS:-CHROM,POS,REF,ALT,INFO/MAVERICK_SCORE}" \
    -Oz \
    -o "$MAVERICK_VCFGZ" \
    "$CURRENT_VCF"
  tabix -f -p vcf "$MAVERICK_VCFGZ"
  CURRENT_VCF="$MAVERICK_VCFGZ"
else
  log "Step 8: skipped MAVERICK because RUN_MAVERICK=0 or MAVERICK_TSV_GZ not set"
fi

FINAL_SMALL_VCFGZ="$OUTDIR/snv/${SAMPLE}.final.small_variants.annotated.vcf.gz"
cp -f "$CURRENT_VCF" "$FINAL_SMALL_VCFGZ"
cp -f "${CURRENT_VCF}.tbi" "${FINAL_SMALL_VCFGZ}.tbi" 2>/dev/null || tabix -f -p vcf "$FINAL_SMALL_VCFGZ"
log "Final SNV/indel annotated VCF: $FINAL_SMALL_VCFGZ"

# ---------- STEP 9: ANNOVAR ----------
if [[ "${RUN_ANNOVAR:-1}" == "1" ]]; then
  need_file "${ANNOVAR_DIR:?Set ANNOVAR_DIR in config}/table_annovar.pl"
  [[ -d "${ANNOVAR_HUMANDB:?Set ANNOVAR_HUMANDB in config}" ]] || die "ANNOVAR_HUMANDB directory does not exist: $ANNOVAR_HUMANDB"
  [[ -n "${ANNOVAR_PROTOCOLS:-}" ]] || die "Set ANNOVAR_PROTOCOLS in config. Example: refGene,clinvar_YYYYMMDD,gnomad41_genome,dbnsfp47a"
  [[ -n "${ANNOVAR_OPERATIONS:-}" ]] || die "Set ANNOVAR_OPERATIONS in config. Example: g,f,f,f"

  log "Step 9: run ANNOVAR table_annovar"
  perl "$ANNOVAR_DIR/table_annovar.pl" \
    "$NORM_VCFGZ" \
    "$ANNOVAR_HUMANDB" \
    -buildver "$ANNOVAR_BUILDVER" \
    -out "$OUTDIR/acmg/${SAMPLE}.annovar" \
    -remove \
    -protocol "$ANNOVAR_PROTOCOLS" \
    -operation "$ANNOVAR_OPERATIONS" \
    -nastring . \
    -vcfinput \
    -polish \
    -thread "$THREADS"
else
  log "Step 9: skipped ANNOVAR because RUN_ANNOVAR=0"
fi

# ---------- STEP 10: InterVar / ACMG-AMP SNV/indel classification ----------
if [[ "${RUN_INTERVAR:-1}" == "1" ]]; then
  [[ -d "${INTERVAR_DIR:?Set INTERVAR_DIR in config}" ]] || die "INTERVAR_DIR does not exist: $INTERVAR_DIR"
  need_file "$INTERVAR_DIR/Intervar.py"
  need_file "$ANNOVAR_DIR/table_annovar.pl"
  need_file "$ANNOVAR_DIR/convert2annovar.pl"
  need_file "$ANNOVAR_DIR/annotate_variation.pl"

  log "Step 10: run InterVar for automated ACMG/AMP SNV/indel classification"
  python3 "$INTERVAR_DIR/Intervar.py" \
    -b "$ANNOVAR_BUILDVER" \
    -i "$NORM_VCFGZ" \
    --input_type=VCF \
    -o "$OUTDIR/acmg/${SAMPLE}.intervar" \
    --table_annovar "$ANNOVAR_DIR/table_annovar.pl" \
    --convert2annovar "$ANNOVAR_DIR/convert2annovar.pl" \
    --annotate_variation "$ANNOVAR_DIR/annotate_variation.pl" \
    -d "$ANNOVAR_HUMANDB" \
    -t "$INTERVAR_DIR/intervardb" || \
      log "InterVar failed or does not support your current build/database setup. Check $OUTDIR/logs."
else
  log "Step 10: skipped InterVar because RUN_INTERVAR=0"
fi

# ---------- CNV/SV branch ----------
# Helper: convert VCF CNVs to ClassifyCNV/ISV-CNV BED format
cnv_to_bed() {
  local input="$1"
  local output="$2"
  case "$input" in
    *.vcf|*.vcf.gz|*.bcf)
      bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\n' "$input" \
        | awk 'BEGIN{OFS="\t"} $3 != "." && ($4=="DEL" || $4=="DUP") {print $1,$2-1,$3,$4}' \
        > "$output"
      ;;
    *.bed|*.bed.gz)
      if [[ "$input" == *.gz ]]; then zcat "$input"; else cat "$input"; fi \
        | awk 'BEGIN{OFS="\t"} NF>=4 && ($4=="DEL" || $4=="DUP") {print $1,$2,$3,$4}' \
        > "$output"
      ;;
    *)
      die "CNV input must be VCF/BCF/BED/BED.GZ"
      ;;
  esac
  [[ -s "$output" ]] || die "No DEL/DUP CNVs were extracted from: $input"
}

if [[ -n "$CNV_INPUT" ]]; then
  need_file "$CNV_INPUT"
  CNV_BED="$OUTDIR/cnv/${SAMPLE}.cnv.classifycnv_isv.bed"
  log "CNV branch: convert CNV input to BED"
  cnv_to_bed "$CNV_INPUT" "$CNV_BED"
  log "CNV BED: $CNV_BED"

  # ---------- STEP 11: AnnotSV ----------
  if [[ "${RUN_ANNOTSV:-1}" == "1" ]]; then
    [[ -d "${ANNOTSV_DIR:-}" ]] || die "Set ANNOTSV_DIR in config"
    need_file "$ANNOTSV_DIR/bin/AnnotSV"
    log "Step 11: run AnnotSV for CNV/SV annotation and ranking"
    "$ANNOTSV_DIR/bin/AnnotSV" \
      -SVinputFile "$CNV_INPUT" \
      -outputFile "$OUTDIR/cnv/${SAMPLE}.AnnotSV.tsv" \
      -genomeBuild "$ASSEMBLY" \
      > "$OUTDIR/logs/${SAMPLE}.AnnotSV.stdout.log" \
      2> "$OUTDIR/logs/${SAMPLE}.AnnotSV.stderr.log" || \
        log "AnnotSV failed. Check logs and verify AnnotSV human annotations are installed."
  else
    log "Step 11: skipped AnnotSV because RUN_ANNOTSV=0"
  fi

  # ---------- STEP 12: ClassifyCNV ----------
  if [[ "${RUN_CLASSIFYCNV:-1}" == "1" ]]; then
    [[ -d "${CLASSIFYCNV_DIR:-}" ]] || die "Set CLASSIFYCNV_DIR in config"
    need_file "$CLASSIFYCNV_DIR/ClassifyCNV.py"
    need_cmd bedtools
    log "Step 12: run ClassifyCNV ACMG-style CNV classification"
    (
      cd "$CLASSIFYCNV_DIR"
      python3 ClassifyCNV.py \
        --infile "$CNV_BED" \
        --GenomeBuild "$CLASSIFY_BUILDVER" \
        --cores "$THREADS" \
        --precise \
        --outdir "${SAMPLE}_ClassifyCNV"
    ) || log "ClassifyCNV failed. Check input BED and ClinGen resource update."
  else
    log "Step 12: skipped ClassifyCNV because RUN_CLASSIFYCNV=0"
  fi

  # ---------- STEP 13: ISV-CNV ----------
  if [[ "${RUN_ISV_CNV:-1}" == "1" ]]; then
    log "Step 13: run ISV-CNV Python package if installed"
    ISV_OUT="$OUTDIR/cnv/${SAMPLE}.ISV_CNV.tsv"
    python3 - "$CNV_BED" "$ISV_OUT" <<'PY'
import sys
import pandas as pd

bed, out = sys.argv[1], sys.argv[2]
cnv = pd.read_csv(bed, sep='\t', header=None, names=['chromosome','start','end','cnv_type'])
try:
    import isv
except Exception as e:
    sys.stderr.write(f"ISV-CNV Python package is not installed or failed to import: {e}\n")
    sys.exit(2)

# ISV expects GRCh38-like columns in many examples: chromosome, start, end, cnv_type.
# Depending on installed package version, wrapper names may differ; this covers the common API.
try:
    if hasattr(isv, 'isv'):
        result = isv.isv(cnv, proba=True, shap=False)
    elif hasattr(isv, 'ISV'):
        model = isv.ISV()
        result = model.predict(cnv)
    else:
        raise AttributeError('Could not find isv.isv() or isv.ISV() API')
    pd.DataFrame(result).to_csv(out, sep='\t', index=False)
except Exception as e:
    sys.stderr.write(f"ISV-CNV run failed: {e}\n")
    sys.exit(3)
PY
    if [[ -s "$ISV_OUT" ]]; then
      log "ISV-CNV output: $ISV_OUT"
    else
      log "ISV-CNV did not produce output. Install/configure the isv package, or set RUN_ISV_CNV=0."
    fi
  else
    log "Step 13: skipped ISV-CNV because RUN_ISV_CNV=0"
  fi

  # ---------- STEP 14: Horizon placeholder ----------
  if [[ "${RUN_HORIZON:-0}" == "1" ]]; then
    if [[ -n "${HORIZON_CMD:-}" ]]; then
      log "Step 14: run local Horizon command provided by HORIZON_CMD"
      # The public literature describes Horizon as an ACMG-aligned CNV model, but
      # there is not a single canonical public CLI at the time this template was written.
      # Set HORIZON_CMD in config to your local command, using {input} and {output} placeholders.
      HORIZON_OUT="$OUTDIR/cnv/${SAMPLE}.Horizon.tsv"
      cmd="${HORIZON_CMD//\{input\}/$CNV_BED}"
      cmd="${cmd//\{output\}/$HORIZON_OUT}"
      bash -lc "$cmd"
      log "Horizon output expected at: $HORIZON_OUT"
    else
      log "RUN_HORIZON=1 but HORIZON_CMD is empty. Add a local command or set RUN_HORIZON=0."
    fi
  else
    log "Step 14: skipped Horizon because RUN_HORIZON=0"
  fi
else
  log "No CNV input supplied with -n; skipped AnnotSV/ClassifyCNV/ISV-CNV/Horizon branch"
fi

# ---------- STEP 15: simple report ----------
REPORT="$OUTDIR/reports/${SAMPLE}.annotation_outputs.txt"
{
  echo "Sample: $SAMPLE"
  echo "Assembly: $ASSEMBLY"
  echo "Final annotated SNV/indel VCF: $FINAL_SMALL_VCFGZ"
  echo "ANNOVAR prefix: $OUTDIR/acmg/${SAMPLE}.annovar"
  echo "InterVar prefix: $OUTDIR/acmg/${SAMPLE}.intervar"
  if [[ -n "$CNV_INPUT" ]]; then
    echo "CNV BED: $OUTDIR/cnv/${SAMPLE}.cnv.classifycnv_isv.bed"
    echo "AnnotSV: $OUTDIR/cnv/${SAMPLE}.AnnotSV.tsv"
    echo "ClassifyCNV: $CLASSIFYCNV_DIR/ClassifyCNV_results/${SAMPLE}_ClassifyCNV"
    echo "ISV-CNV: $OUTDIR/cnv/${SAMPLE}.ISV_CNV.tsv"
    echo "Horizon: $OUTDIR/cnv/${SAMPLE}.Horizon.tsv"
  fi
  echo "Pipeline log: $LOGFILE"
} > "$REPORT"

log "Finished. Output summary: $REPORT"
