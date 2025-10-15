#!/bin/bash
set -euo pipefail
shopt -s nullglob

threads=8

# --- References and resources ---
ref="supporting_files/hg38/hg38.fa"
dbsnp_vcf="supporting_files/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
known_sites="$dbsnp_vcf"
targets40_vcf="supporting_files/targets_40.vcf.gz"

# --- sample names ---
RECEPTOR= [...]
DONOR= [...]

# --- folders ---
raw_fastq="0_raw_fastq"
reads="1_filtered_fastq"
aligned_reads="2_aligned_reads"
dedup="3_dedup_reads"
recal="4_recal"
bqsr="5_bqsr"
size_metrics="6_size_metrics"
results="results"

mkdir -p "$reads" "$aligned_reads" "$dedup" "$recal" "$bqsr" "$size_metrics" "$results"

# --- Intermediate results ---
MIN_DP=10
MIN_ALT_READS=3
OUTDIR_MT2="${results}/mutect2"
OUTDIR_TAB="${results}/tables"
OUTDIR_SIG="${results}/signatures"
OUTDIR_SCR="${results}/scripts"
mkdir -p "$OUTDIR_MT2" "$OUTDIR_TAB" "$OUTDIR_SIG" "$OUTDIR_SCR"

# 0) PRE-CHECK
echo "# 0. PRE-CHECK: referencias e índices"
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"
[[ -f "${ref%.*}.dict" ]] || gatk CreateSequenceDictionary -R "$ref" -O "${ref%.*}.dict"
[[ -f "${dbsnp_vcf}.tbi" ]] || tabix -p vcf "$dbsnp_vcf"
[[ -f "${targets40_vcf}.tbi" ]] || tabix -p vcf "$targets40_vcf"
if [[ ! -f "${ref}.bwt" || ! -f "${ref}.sa" ]]; then
  echo "# Indexando referencia para BWA (primera vez)"; bwa index -a bwtsw "$ref"; fi

# 1) FASTP
echo "# 1. FASTQ QC / TRIMMING WITH FASTP"
R1_list=("${raw_fastq}"/*_R1_001.fastq.gz)
if (( ${#R1_list[@]} == 0 )); then
  echo "[WARN] No hay FASTQ en ${raw_fastq}; saltando FASTP"
else
  for R1 in "${R1_list[@]}"; do
    sample=$(basename "$R1" _R1_001.fastq.gz)
    R2="${raw_fastq}/${sample}_R2_001.fastq.gz"
    echo "Procesando muestra: $sample"
    if [[ ! -f "$R2" ]]; then echo "[ERROR] Falta R2: $R2"; exit 1; fi
    fastp -i "$R1" -I "$R2" \
      -o "${reads}/${sample}_R1_001.filtered.fastq.gz" \
      -O "${reads}/${sample}_R2_001.filtered.fastq.gz" \
      --max_len1 110 --max_len2 110 --disable_quality_filtering \
      --thread "$threads" \
      -h "${reads}/${sample}.fastp_report.html" -j "${reads}/${sample}.fastp_report.json"
  done
fi

# 2) BWA-MEM
echo "# 2. ALIGNMENT WITH BWA-MEM"
F1_list=("${reads}"/*_R1_001.filtered.fastq.gz)
if (( ${#F1_list[@]} == 0 )); then
  echo "[ERROR] No hay FASTQ filtrados en ${reads}. ¿Se ejecutó FASTP?"; exit 1
fi
for R1 in "${F1_list[@]}"; do
  sample=$(basename "$R1" _R1_001.filtered.fastq.gz)
  R2="${reads}/${sample}_R2_001.filtered.fastq.gz"
  [[ -f "$R2" ]] || { echo "[ERROR] Falta $R2"; exit 1; }
  RG="@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}"
  bwa mem -t "$threads" -R "$RG" "$ref" "$R1" "$R2" > "${aligned_reads}/${sample}.sam"
 done

# Comprobación de SAMs creados
sam_files=("${aligned_reads}"/*.sam)
if (( ${#sam_files[@]} == 0 )); then
  echo "[ERROR] No se generaron SAMs en ${aligned_reads}"; exit 1
else
  echo "SAMs creados:"; ls -1 "${aligned_reads}"/*.sam
fi

# 3) MARK DUPLICATES
echo "# 3. MARK DUPLICATES"
for sam in "${sam_files[@]}"; do
  sample=$(basename "$sam" .sam)
  samtools view -@ "$threads" -bS "$sam" | samtools sort -@ "$threads" -o "${dedup}/${sample}_sorted.bam"
  gatk MarkDuplicates -I "${dedup}/${sample}_sorted.bam" -O "${dedup}/${sample}_sorted_dedup.bam" -M "${dedup}/${sample}_dup_metrics.txt" --CREATE_INDEX true
 done

# 4) BQSR
echo "# 4. BASE QUALITY SCORE RECALIBRATION (BQSR)"
for bam in "$dedup"/*_sorted_dedup.bam; do
  sample=$(basename "$bam" _sorted_dedup.bam)
  gatk BaseRecalibrator -R "$ref" -I "$bam" --known-sites "$known_sites" -O "${recal}/${sample}_recal_data.table"
  gatk ApplyBQSR -R "$ref" -I "$bam" --bqsr-recal-file "${recal}/${sample}_recal_data.table" -O "${bqsr}/${sample}_bqsr.bam"
  samtools index "${bqsr}/${sample}_bqsr.bam"
 done

# 5) Métricas de QC profundas (FastQC + samtools stats + Picard + MultiQC)
echo "# 5. Métricas de QC profundas (FastQC + samtools stats + Picard + MultiQC)"

# --- FastQC en RAW y FILTRADOS ---
fq_raw_out="${size_metrics}/fastqc_raw"
fq_filt_out="${size_metrics}/fastqc_filtered"
mkdir -p "$fq_raw_out" "$fq_filt_out"

# RAW
for fq in "${raw_fastq}"/*.fastq.gz; do
  [[ -e "$fq" ]] || continue
  base=$(basename "${fq}" .fastq.gz)
  if [[ ! -s "${fq_raw_out}/${base}_fastqc.zip" ]]; then
    fastqc -t "${threads}" -o "$fq_raw_out" "$fq"
  else
    echo "[SKIP] FastQC RAW ya existe: ${base}"
  fi
done

# FILTRADOS
for fq in "${reads}"/*_R[12]_001.filtered.fastq.gz; do
  [[ -e "$fq" ]] || continue
  base=$(basename "${fq}" .fastq.gz)
  if [[ ! -s "${fq_filt_out}/${base}_fastqc.zip" ]]; then
    fastqc -t "${threads}" -o "$fq_filt_out" "$fq"
  else
    echo "[SKIP] FastQC FILT ya existe: ${base}"
  fi
done

# --- samtools stats por BAM ---
for bam in "${bqsr}"/*_bqsr.bam; do
  [[ -e "$bam" ]] || continue
  sample=$(basename "$bam" _bqsr.bam)
  stats_out="${size_metrics}/${sample}.samtools.stats"
  if [[ ! -s "$stats_out" ]]; then
    samtools stats -@ "${threads}" "$bam" > "$stats_out"
  else
    echo "[SKIP] samtools stats ya existe: ${stats_out}"
  fi
done

# --- Picard: GC bias, base distribution by cycle, quality yield ---
for bam in "${bqsr}"/*_bqsr.bam; do
  [[ -e "$bam" ]] || continue
  sample=$(basename "$bam" _bqsr.bam)

  gc_sum="${size_metrics}/${sample}.gc_bias_summary.txt"
  gc_chart="${size_metrics}/${sample}.gc_bias_chart.pdf"
  gc_metrics="${size_metrics}/${sample}.gc_bias_metrics.txt"

  base_cycle="${size_metrics}/${sample}.base_dist_by_cycle.txt"
  base_cycle_pdf="${size_metrics}/${sample}.base_dist_by_cycle.pdf"

  qy_metrics="${size_metrics}/${sample}.quality_yield_metrics.txt"

  if [[ ! -s "$gc_metrics" || ! -s "$gc_chart" || ! -s "$gc_sum" ]]; then
    gatk CollectGcBiasMetrics -I "$bam" -R "$ref" -O "$gc_metrics" -S "$gc_sum" --CHART "$gc_chart"
  else
    echo "[SKIP] Picard GC Bias ya existe para $sample"
  fi

  if [[ ! -s "$base_cycle" || ! -s "$base_cycle_pdf" ]]; then
    gatk CollectBaseDistributionByCycle -I "$bam" -O "$base_cycle" --CHART "$base_cycle_pdf"
  else
    echo "[SKIP] BaseDistributionByCycle ya existe para $sample"
  fi

  if [[ ! -s "$qy_metrics" ]]; then
    gatk CollectQualityYieldMetrics -I "$bam" -O "$qy_metrics"
  else
    echo "[SKIP] QualityYieldMetrics ya existe para $sample"
  fi
done

# --- MultiQC: agregamos todo ---
mqc_out="${size_metrics}"
mqc_report="${mqc_out}/multiqc_report.html"
if [[ ! -s "$mqc_report" ]]; then
  multiqc \
    "$fq_raw_out" \
    "$fq_filt_out" \
    "$reads" \
    "$size_metrics" \
    -o "$mqc_out"
else
  echo "[SKIP] MultiQC ya existe: $mqc_report"
fi

# 6) Firma R & D
echo "# 6. FIRMA R&D"

gatk HaplotypeCaller -R "$ref" -I "${bqsr}/${RECEPTOR}_bqsr.bam" -ERC GVCF -L "$targets40_vcf" -O "${OUTDIR_SIG}/${RECEPTOR}.g.vcf.gz" --dbsnp "$dbsnp_vcf" --native-pair-hmm-threads "$threads"
gatk HaplotypeCaller -R "$ref" -I "${bqsr}/${DONOR}_bqsr.bam"     -ERC GVCF -L "$targets40_vcf" -O "${OUTDIR_SIG}/${DONOR}.g.vcf.gz"    --dbsnp "$dbsnp_vcf" --native-pair-hmm-threads "$threads"

gatk CombineGVCFs -R "$ref" -V "${OUTDIR_SIG}/${RECEPTOR}.g.vcf.gz" -V "${OUTDIR_SIG}/${DONOR}.g.vcf.gz" -L "$targets40_vcf" -O "${OUTDIR_SIG}/RD.cohort.g.vcf.gz"

gatk GenotypeGVCFs -R "$ref" -V "${OUTDIR_SIG}/RD.cohort.g.vcf.gz" -L "$targets40_vcf" -O "${OUTDIR_SIG}/RD.signature.vcf.gz"

tabix -p vcf "${OUTDIR_SIG}/RD.signature.vcf.gz" || true

gatk VariantsToTable -V "${OUTDIR_SIG}/RD.signature.vcf.gz" -O "${OUTDIR_SIG}/RD.signature.tsv" -F CHROM -F POS -F ID -F REF -F ALT -GF GT --show-filtered true

# 7) Donor/Receptor map extendido
echo "# 7. DONOR/RECEPTOR MAP"
cat > "${OUTDIR_SCR}/make_donor_map.py" << 'PY'
import sys, gzip

vcf_path, receptor, donor, out_tsv = sys.argv[1:5]

open_vcf = (lambda p: gzip.open(p, 'rt')) if vcf_path.endswith('.gz') else open

with open_vcf(vcf_path) as f, open(out_tsv, 'w') as out:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            samples = header[9:]
            assert receptor in samples and donor in samples, 'Receptor o Donor no en VCF'
            rec_i = samples.index(receptor) + 9
            don_i = samples.index(donor) + 9
            out.write('\t'.join([
                'CHROM','POS','ID','REF','ALTS','GT_R','GT_D',
                'donor_alt_indices','receptor_alt_indices','informative_type'
            ]) + '\n')
            continue

        fields = line.strip().split('\t')
        chrom, pos, vid, ref, alts, _, _, _, fmt = fields[:9]
        alts_list = alts.split(',') if alts != '.' else []
        fmt_keys = fmt.split(':')

        def parse_gt(s):
            vals = s.split(':')
            d = dict(zip(fmt_keys, vals))
            gt = d.get('GT', './.')
            alle = [a for a in gt.replace('|','/').split('/') if a != '.']
            return gt, set(int(a) for a in alle if a.isdigit())

        gt_r, set_r = parse_gt(fields[rec_i])
        gt_d, set_d = parse_gt(fields[don_i])

        donor_only = set()
        receptor_only = set()

        if 0 in set_d and 0 not in set_r:
            donor_only.add(0)
        if 0 in set_r and 0 not in set_d:
            receptor_only.add(0)
        for a in set_d:
            if a != 0 and a not in set_r:
                donor_only.add(a)
        for a in set_r:
            if a != 0 and a not in set_d:
                receptor_only.add(a)

        donor_idx_sorted = sorted(donor_only)
        receptor_idx_sorted = sorted(receptor_only)

        if donor_idx_sorted and receptor_idx_sorted:
            itype = 'BOTH'
        elif donor_idx_sorted:
            itype = 'DONOR_ONLY'
        elif receptor_idx_sorted:
            itype = 'RECEPTOR_ONLY'
        else:
            itype = 'NONE'

        out.write('\t'.join([
            chrom, pos, vid, ref, alts, gt_r, gt_d,
            ','.join(map(str, donor_idx_sorted)),
            ','.join(map(str, receptor_idx_sorted)),
            itype
        ]) + '\n')
PY


python3 "${OUTDIR_SCR}/make_donor_map.py" "${OUTDIR_SIG}/RD.signature.vcf.gz" "$RECEPTOR" "$DONOR" "${OUTDIR_SIG}/donor_map.tsv"
echo "[INFO] donor_map.tsv -> ${OUTDIR_SIG}/donor_map.tsv"

# 8) Detectar mezclas
echo "# 8. DETECTAR MEZCLAS"
MIX_BAMS=( $(ls "${bqsr}"/*_bqsr.bam | grep -v -E "/(${RECEPTOR}_bqsr|${DONOR}_bqsr)\.bam$" || true) )
if (( ${#MIX_BAMS[@]} == 0 )); then echo "[ERROR] No hay BAMs de mezcla en ${bqsr}"; exit 1; fi

echo "[INFO] Mezclas detectadas:"; for x in "${MIX_BAMS[@]}"; do echo "  - $(basename "$x")"; done

# 9) Mutect2 por mezcla
echo "# 9. MUTECT2 POR MEZCLAS"
for mix_bam in "${MIX_BAMS[@]}"; do
  mix_base=$(basename "$mix_bam" _bqsr.bam)
  gatk Mutect2 -R "$ref" -I "$mix_bam" -tumor "$mix_base" -I "${bqsr}/${RECEPTOR}_bqsr.bam" -normal "$RECEPTOR" \
    -L "$targets40_vcf" --alleles "${OUTDIR_SIG}/RD.signature.vcf.gz" --genotype-germline-sites true --max-reads-per-alignment-start 0 \
    -O "${OUTDIR_MT2}/${mix_base}.unfiltered.vcf.gz"
  tabix -p vcf "${OUTDIR_MT2}/${mix_base}.unfiltered.vcf.gz" || true
  bcftools query -s "$mix_base" -f '%CHROM	%POS	%REF	%ALT	[%AD]	[%DP]
' "${OUTDIR_MT2}/${mix_base}.unfiltered.vcf.gz" > "${OUTDIR_MT2}/${mix_base}.AD.tsv"
 done

# 10) SITIOS ANCLA (hetero + dosificación)
echo "# 10. SITIOS ANCLA (desde donor_map con heterocigosidad y dosificación)"
cat > "${OUTDIR_SCR}/create_anchor_sites_from_map.py" << 'PY'
import sys, csv, argparse

def parse_idx_list(s):
    return [int(x) for x in s.split(',')] if s else []

def parse_gt(gt):
    gt = gt.replace('|','/')
    return [int(a) for a in gt.split('/') if a.isdigit()]

def dosage_for(indices, gt_list):
    return sum(1 for a in gt_list if a in set(indices))

p=argparse.ArgumentParser()
p.add_argument('--donormap', required=True)
p.add_argument('--out', required=True)
args=p.parse_args()

rows=[]
with open(args.donormap) as f:
    r=csv.DictReader(f, delimiter='\t')
    for row in r:
        didx=parse_idx_list(row.get('donor_alt_indices',''))
        ridx=parse_idx_list(row.get('receptor_alt_indices',''))
        gtR=parse_gt(row.get('GT_R',''))
        gtD=parse_gt(row.get('GT_D',''))

        if didx:
            mode='DONOR'; idx_used=didx; dosage=dosage_for(idx_used, gtD)
        elif ridx:
            mode='RECEPTOR'; idx_used=ridx; dosage=dosage_for(idx_used, gtR)
        else:
            continue
        if dosage==0: continue

        rows.append((row['CHROM'], row['POS'], ','.join(map(str,idx_used)), mode, str(dosage), row.get('GT_R',''), row.get('GT_D',''), '', '', '50.0000'))

with open(args.out,'w') as out:
    w=csv.writer(out, delimiter='\t')
    w.writerow(['CHROM','POS','idx_used','mode','dosage_used','GT_R','GT_D','anchor_reads_sel','anchor_DP','anchor_pct_corr'])
    for rr in rows: w.writerow(rr)
PY

python3 "${OUTDIR_SCR}/create_anchor_sites_from_map.py" \
  --donormap "${OUTDIR_SIG}/donor_map.tsv" \
  --out "${OUTDIR_SIG}/anchor_sites.tsv"

N_ANCHOR=$(($(wc -l < "${OUTDIR_SIG}/anchor_sites.tsv")-1))
echo "[INFO] Sitios ancla (heterocigosidad/dosis): ${OUTDIR_SIG}/anchor_sites.tsv (${N_ANCHOR} sitios)"

# 11) Resumen % Donante (anclado + ponderado por DP)
echo "# 11. RESUMEN % DONANTE (ancla y 0% si no detecta)"
cat > "${OUTDIR_SCR}/summarize_mixes_anchored.py" << 'PY'
import sys, csv, os, argparse

def load_anchor(path):
    keys = []; D = {}
    with open(path, newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            key = (row['CHROM'], row['POS'])
            idx = [int(x) for x in row['idx_used'].split(',')] if row['idx_used'] else []
            D[key] = {
                'idx': idx,
                'mode': row['mode'],
                'dosage': int(row['dosage_used']),
                'GT_R': row.get('GT_R', ''),
                'GT_D': row.get('GT_D', '')
            }
            keys.append(key)
    return keys, D

p = argparse.ArgumentParser()
p.add_argument('--anchor_sites', required=True)
p.add_argument('--ad_files', nargs='+', required=True)
p.add_argument('--min_dp', type=int, default=10)
p.add_argument('--min_alt_reads', type=int, default=3)
p.add_argument('--outdir_tables', required=True)
p.add_argument('--summary_out', required=True)
a = p.parse_args()

ANCHOR_KEYS, ANCHOR = load_anchor(a.anchor_sites)
os.makedirs(a.outdir_tables, exist_ok=True)
ALL = []

for path in a.ad_files:
    mix = os.path.basename(path).replace('.AD.tsv','')
    AD = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            chrom, pos, ref, alts, ad_str, dp_str = line.split('\t')
            key = (chrom, pos)
            AD_list = [int(x) for x in ad_str.split(',')] if ad_str not in ('.','') else []
            DP = int(dp_str) if dp_str not in ('.','') else sum(AD_list)
            AD[key] = {'AD': AD_list, 'DP': DP, 'REF': ref, 'ALTS': alts}

    per_site = os.path.join(a.outdir_tables, f'{mix}.per_site.tsv')
    with open(per_site, 'w', newline='') as fo:
        w = csv.writer(fo, delimiter='\t')
        w.writerow([
            'MIX','CHROM','POS','mode','idx_used','dosage_used',
            'GT_R','GT_D',
            'AD_list','DP','reads_sel','pct_raw','pct_corr','status'
        ])
        sum_raw = 0.0
        sum_corr = 0.0
        n_sites = 0
        W = 0
        SUMW = 0.0

        for key in ANCHOR_KEYS:
            info = ANCHOR[key]
            rec = AD.get(key, {'AD': [], 'DP': 0})
            AD_list, DP = rec['AD'], rec['DP']
            reads_sel = sum(AD_list[i] for i in info['idx'] if i < len(AD_list)) if AD_list else 0

            if DP < a.min_dp or reads_sel < a.min_alt_reads or info['dosage'] == 0:
                f_raw = 0.0
                f_corr = 0.0
                status = 'ZERO_BY_ANCHOR'
            else:
                f_raw = reads_sel / DP
                f_corr = f_raw * (2.0 / info['dosage'])
                if info['mode'] == 'RECEPTOR':
                    f_raw = 1.0 - f_raw
                    f_corr = 1.0 - f_corr
                # Limitar por seguridad a [0,1]
                f_corr = max(0.0, min(1.0, f_corr))
                status = 'DETECTED'

            w.writerow([
                mix, key[0], key[1], info['mode'], ','.join(map(str, info['idx'])),
                info['dosage'], info['GT_R'], info['GT_D'],
                ','.join(map(str, AD_list)) if AD_list else '', DP, reads_sel,
                f"{100*f_raw:.4f}", f"{100*f_corr:.4f}", status
            ])

            sum_raw += f_raw
            sum_corr += f_corr
            n_sites += 1
            W += DP
            SUMW += f_corr * DP

        anchored_mean_raw  = (100.0 * sum_raw / n_sites) if n_sites > 0 else 0.0
        anchored_mean_corr = (100.0 * sum_corr / n_sites) if n_sites > 0 else 0.0
        dp_weighted_corr   = (100.0 * SUMW / W) if W > 0 else 0.0

    ALL.append([mix, n_sites, anchored_mean_raw, anchored_mean_corr, dp_weighted_corr])

with open(a.summary_out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['MIX','n_anchor_sites','anchored_mean_raw','anchored_mean_corr','dp_weighted_corr'])
    for row in ALL:
        w.writerow(row)
PY

python3 "${OUTDIR_SCR}/summarize_mixes_anchored.py" \
  --anchor_sites "${OUTDIR_SIG}/anchor_sites.tsv" \
  --ad_files ${OUTDIR_MT2}/*.AD.tsv \
  --min_dp "$MIN_DP" \
  --min_alt_reads "$MIN_ALT_READS" \
  --outdir_tables "$OUTDIR_TAB" \
  --summary_out "${OUTDIR_TAB}/final_summary.tsv"

# FIN
echo -e "
===============================================
[DONE] Resultados clave (v8):
  - Firma R&D (GTs):        ${OUTDIR_SIG}/RD.signature.tsv
  - Mapa donante/receptor:  ${OUTDIR_SIG}/donor_map.tsv
  - Sitios ancla (hom-op):  ${OUTDIR_SIG}/anchor_sites.tsv
  - VCF Mutect2 por mezcla: ${OUTDIR_MT2}/*unfiltered.vcf.gz
  - AD por mezcla:          ${OUTDIR_MT2}/*.AD.tsv
  - Por sitio por mezcla:   ${OUTDIR_TAB}/*.per_site.tsv
  - Resumen anclado:        ${OUTDIR_TAB}/final_summary.tsv
  - MultiQC (profundo):     ${size_metrics}/multiqc_report.html
==============================================="
