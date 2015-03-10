The original question is [here](https://www.biostars.org/p/133825). In short, someone aligned fastq files with phred+64 quality encodings without telling the aligner that. In an ideal world, one would simply redo the alignments, though in reality they're not likely to change much if at all and it'd be convenient to simply modify the SAM/BAM/CRAM file to convert the quality scores to phred+33 (otherwise, GATK breaks).

The accompanying script is meant to perform this conversion. It will accept a SAM/BAM/CRAM file and produce either a BAM or CRAM file, converting the quality scores as it goes. This script will assume that quality is incorrectly encoded as phred+64. If a SAM/BAM/CRAM file with proper quality score encoding is input, then this program will likely abort with an error ("likely" only because if there are no low quality scores then it's impossible to determine the encoding.

You can install this script by:

    git clone https://github.com/dpryan79/Answers.git
    cd Answers
    git submodule init
    git submodule update
    cd biostars_133825
    make
    mv ConvertPhredQuals /somewhere/in/your/path

Obviously, "/somewhere/in/your/path" should be changed.
