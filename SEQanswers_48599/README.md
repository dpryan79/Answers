Inspired by http://seqanswers.com/forums/showthread.php?t=48599

The task is to take a sorted/index BAM file or RNAseq data and produce two WIG files with per-base coverage that handles strand correctly (assuming strand-specific data). This solution uses the htslib C API as a possible example of doing this.
