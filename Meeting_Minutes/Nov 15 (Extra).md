Confirmed that the feature classifier generation code is correct! Reported to Chris that I ran QIIME2 downstreat analysis up until the human-readable file generation with the fixed dataset and rep_features and stats.qza given by the TAs, but cannot proceed that step due tomissing table.qza output

Chris provided the code necessary to rerun denoise:

qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux_seqs.qza \
--p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG \
--p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-discard-untrimmed \
--o-trimmed-sequences 1a_demux_seqs-trimmed.qza \
--verbose > cutadapt-log-2.txt

screen -S denoise
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 1a_demux_seqs-trimmed.qza \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

For next week, re-run denoise and generate the output necessary for R analysis
