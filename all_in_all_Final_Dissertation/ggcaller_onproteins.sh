ggcaller \
  -i ggcaller_input/genome1.faa \
     ggcaller_input/genome2.faa \
     ggcaller_input/genome3.faa \
     ggcaller_input/genome4.faa \
  --identity-cutoff 0.95 \
  --len-diff-cutoff 0.9 \
  --family-threshold 0.7 \
  --alignment pan \
  --aligner def \
  --out ggc_mtbc \
  --threads 8