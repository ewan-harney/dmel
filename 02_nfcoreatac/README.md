# 02_nfcoreatac
Within this 02_nfcoreatac folder is 1 subdirectory, which contains further subdirectories

- macs2
  - broad_peak subdirectoryy
    - consensus (subdirectory)
    - qc (subdirectory)
    - various peaks files for all samples

Abrreviations appearing in individual file names: 
Manz = Manzanares
Heat = Heat shocked (F3 only)
Ctrl = Control (F3 or F6)
HSan = Heat shocked in the ancestral generation (F6 only)

These output were generated by the running the nf-core/atac pipeline (see InterChromaTE_nf-core-atacseq.sh script). Peaks files provide a useful summary of atac-seq results for each sample. However subsequent intepareto analyses used individual sample .bam files created by the pipeline. Beacuse these files are so large (uncompressed size > 24 Gb) they are not provided in this data repository.

Within the broad_peak subdirectory are four files  for ever sample
".mLb.clN_peaks" summary of locations and characteristics of accessible chromatin broad peaks
".mLb.clN_peaks.broadPeak" broadPeaks are more informative for general signals of chromatin accessibility
".mLb.clN_peaks.gappedPeak" gappedPeak is a representation of narrow peaks as blocks over a broad peak
".mLb.clN_peaks.annotatePeaks" annotations of broad peaks

Within the broad_peak subdirectory is a consensus subdirectory that contains various output files concerning consensus peaks across all samples, including:
- consensus_peaks.mLb.clN.bed : This bed file provides an overview of all possible open chromatin peaks in the data.

See https://nf-co.re/atacseq/2.1.2/ for more information on nf-core/atacseq
