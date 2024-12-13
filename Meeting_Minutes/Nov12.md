# November 12, 2024

## Meeting Agenda
- QIIME troubleshooting

## Meeting Notes
- Confirmed that the feature classifier generation code is correct
- QIIME2 downstream analysis was ran up until the human-readable file generation with the fixed dataset and rep_features and stats.qza given by the TAs, but cannot proceed that step due to missing table.qza output

- Chris provided the code necessary to rerun denoise:

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

## Action Items
- Re-run denoise and generate the output necessary for R analysis
