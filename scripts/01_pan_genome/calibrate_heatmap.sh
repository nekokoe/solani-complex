#!/usr/bin/env bash
# Calibrate pheatmap cell sizes so the heatmap BODY renders SQUARE and at a
# target size (>= 2x the original). Iterates: render -> measure body bbox (via
# Python/PIL on a 300-dpi render) -> correct cellheight (aspect) and scale both
# cell sizes (size) -> re-render until converged.
set -euo pipefail

BASE="/share/pasteur/gaosh/Fusarium/phase2/solaniPathogenicityPaper"
TARGET_IN=10.5          # target body width in inches (>=2x the ~5in original)
RES=150

MEASURE=$(cat <<'PY'
from PIL import Image
import numpy as np
from scipy import ndimage
import sys
im = np.array(Image.open(sys.argv[1]).convert('RGB')).astype(int)
purple = np.array([0x44,0x01,0x54]); yellow = np.array([0xFD,0xE7,0x25])
m = (np.abs(im-purple).sum(axis=2) < 200) | (np.abs(im-yellow).sum(axis=2) < 200)
m = ndimage.binary_dilation(m, iterations=3)
ys, xs = np.where(m)
gw = xs.max()-xs.min()+1; gh = ys.max()-ys.min()+1
print(f"{gw} {gh}")
PY
)

calibrate() {
  local sp="$1"
  local ncol nrow
  ncol=$(wc -l < "$BASE/09.rf_clinical_genes/top_significant_hogs.$sp.list")
  nrow=$(wc -l < "$BASE/05.pan-genome/classify_med_env_shorter/all.$sp.sample.list")
  local cw=0.9 ch
  ch=$(python3 -c "print(f'{$cw*$ncol/$nrow:.4f}')")
  for i in 1 2 3 4 5 6; do
    Rscript "$BASE/scripts/plot_rf_heatmap.R" \
      --feature-matrix "$BASE/05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv" \
      --top-features "$BASE/09.rf_clinical_genes/top_significant_hogs.$sp.list" \
      --sample-list "$BASE/05.pan-genome/classify_med_env_shorter/all.$sp.sample.list" \
      --group-file "$BASE/bmu_ref_groups.211.txt" \
      --effector-hog "$BASE/08.pan_020_genes/effector_master.hog_level.tsv" \
      --species "$sp" --fig-w 24 --fig-h 24 \
      --cellwidth "$cw" --cellheight "$ch" \
      --output-prefix "/tmp/hmcal_$sp" >/dev/null 2>&1
    pdftoppm -png -r "$RES" -f 1 -l 1 "/tmp/hmcal_$sp.pdf" "/tmp/hmc_$sp" >/dev/null 2>&1
    read -r gw gh <<< "$(python3 -c "$MEASURE" /tmp/hmc_$sp-1.png)"
    echo "  $sp iter$i: cw=$cw ch=$ch -> body ${gw}x${gh}px ($(python3 -c "print(f'{$gw/$RES:.1f}')")x$(python3 -c "print(f'{$gh/$RES:.1f}')")in, aspect $(python3 -c "print(f'{$gw/$gh:.3f}')"))"
    # 1) fix aspect via cellheight
    ch=$(python3 -c "print(f'{$ch*$gw/$gh:.4f}')")
    # 2) fix size via both cells
    scale=$(python3 -c "print(f'{$TARGET_IN*$RES/$gw:.4f}')")
    cw=$(python3 -c "print(f'{$cw*$scale:.4f}')")
    ch=$(python3 -c "print(f'{$ch*$scale:.4f}')")
    if python3 -c "exit(0 if 0.97 <= $gw/$gh <= 1.03 else 1)"; then
      # verify size once more
      Rscript "$BASE/scripts/plot_rf_heatmap.R" \
        --feature-matrix "$BASE/05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv" \
        --top-features "$BASE/09.rf_clinical_genes/top_significant_hogs.$sp.list" \
        --sample-list "$BASE/05.pan-genome/classify_med_env_shorter/all.$sp.sample.list" \
        --group-file "$BASE/bmu_ref_groups.211.txt" \
        --effector-hog "$BASE/08.pan_020_genes/effector_master.hog_level.tsv" \
        --species "$sp" --fig-w 24 --fig-h 24 \
        --cellwidth "$cw" --cellheight "$ch" \
        --output-prefix "$BASE/09.rf_clinical_genes/heatmap_rf.$sp" >/dev/null 2>&1
      pdftoppm -png -r "$RES" -f 1 -l 1 "$BASE/09.rf_clinical_genes/heatmap_rf.$sp.pdf" "/tmp/hmf_$sp" >/dev/null 2>&1
      read -r gw gh <<< "$(python3 -c "$MEASURE" /tmp/hmf_$sp-1.png)"
      echo "  FINAL $sp: cw=$cw ch=$ch -> body ${gw}x${gh}px ($(python3 -c "print(f'{$gw/$RES:.1f}')")x$(python3 -c "print(f'{$gh/$RES:.1f}')")in, aspect $(python3 -c "print(f'{$gw/$gh:.3f}')"))"
      return 0
    fi
  done
  echo "  WARNING: $sp did not converge"
}

for sp in solani falciformis; do
  echo "== $sp =="
  calibrate "$sp"
done
