This is from: https://www.biostars.org/p/112049/

The BSMAP program comes with a methylation extractor that (A) uses features in samtools that are no longer supported and (B) will output nothing if --pair is used (denoting to discard discordant and singleton alignments) with a single-end dataset. This are both easy fixes and this (only slightly tested) version of methratio.py should fix these issues.
