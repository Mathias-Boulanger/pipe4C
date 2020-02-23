VERSION <- '1.50'

get_script_path <- function(path=NULL) {
	if(is.null(path)){
		cmdArgs = commandArgs(trailingOnly = FALSE)
		needle = "--file="
		match = grep(needle, cmdArgs)
		if (length(match) > 0) {
			# Rscript
			return(normalizePath(sub(needle, "", cmdArgs[match])))
		} else {
			ls_vars = ls(sys.frames()[[1]])
			if ("fileName" %in% ls_vars) {
				# Source'd via RStudio
				return(normalizePath(sys.frames()[[1]]$fileName))
			} else {
				# Source'd via R console
				return(normalizePath(sys.frames()[[1]]$ofile))
			}
		}
	} else {
		return(path)
	}
}

#################################################################################################################
### PARSING THE INPUT ###########################################################################################
#################################################################################################################
if(!suppressMessages(require("optparse", character.only = TRUE))) stop("Package not found: optparse")

option_list = list(
	make_option(c("-v", "--vpFile"), type="character", default=NULL, 
							help="path to viewpoint file [required]", metavar="/path/to/vp_file"),
	make_option(c("-f", "--fqFolder"), type="character", default=NULL, 
							help="path to folder containing the FASTQ files [required]", metavar="/path/to/FASTQ_folder/"),
	make_option(c("-o", "--outFolder"), type="character", default=NULL, 
							help="path to the output folder [required]", metavar="/path/to/output_folder/"),
	make_option(c("-c", "--confFile"), type="character", default="conf.yml", 
							help="path to configuration file [default %default]", metavar="/path/to/conf.yml/"),

	
	make_option(c("-d", "--mismatchMax"), type="integer", default=NULL,
							help="Maximum number of mismatches allowed with primer sequence during demultiplexing.",metavar="number"),  
	make_option(c("-q", "--qualityCutoff"), type="integer", default=NULL,
							help="Q-score. Trim 3-end of all sequences using a sliding window as soon as 2 of 5 nucleotides has quality encoding less than the Q-score.",metavar="number"),
	make_option(c("-l", "--trimLength"), type="integer", default=NULL,
							help="Trim reads to defined capture length from 3-end.",metavar="number"),
	make_option(c("-m", "--minAmountReads"), type="integer", default=NULL,
							help="Minimum amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed.",metavar="number"),
	make_option(c("-r", "--readsQuality"), type="integer", default=NULL,
							help="Bowtie2 minimum quality mapped reads.",metavar="number"),
	make_option(c("-z", "--cores"), type="integer", default=NULL,
							help="Number of cores for parallelization.",metavar="number"),
	make_option(c("-s", "--wSize"), type="integer", default=NULL,
							help="The running mean window size.",metavar="number"),
	make_option(c("-n", "--nTop"), type="integer", default=NULL,
							help="Top fragments discarded for normalization.",metavar="number"),

	make_option(c("-vp", "--vpRegion"), type="integer", default=NULL,
							help="Windows of observation for the peakC analysis.",metavar="number"),
	make_option(c("-a", "--alphaFDR"), type="integer", default=NULL,
							help="alphaFDR for the peakC analysis.",metavar="number"),
	make_option(c("-qwd", "--qWd"), type="integer", default=NULL,
							help="qWd for the peakC analysis.",metavar="number"),
	make_option(c("-qwr", "--qWr"), type="integer", default=NULL,
							help="qWr for the peakC analysis.",metavar="number"),
	make_option(c("-md", "--minDist"), type="integer", default=NULL,
							help="minDist for the peakC analysis.",metavar="number"),
	make_option(c("-mg", "--min.gapwidth"), type="integer", default=NULL,
							help="min distance between significative peakC peaks to create 2 different regions of interaction.",metavar="number"),
	make_option(c("-fs", "--fixedStepWigBin"), type="integer", default=NULL,
							help="bin to use to write the fixedStep Wig file.",metavar="number"),


	make_option(c("-u", "--mapUnique"), action="store_true", default=FALSE,
							help="Extract uniquely mapped reads, based on the lack of XS tag."),
	make_option(c("-nb", "--nonBlind"), action="store_true", default=FALSE,
							help="Only keep non-blind fragments"),
	make_option(c("-bdg", "--bdg"), action="store_true", default=FALSE,
							help="create bedGraph files for all samples"),	
	make_option(c("-w", "--wig"), action="store_true", default=FALSE,
							help="create wig files for all samples"),
	make_option(c("-wf", "--fixedStepWig"), action="store_true", default=FALSE,
							help="create fixedStepWig files for all samples"),
	
	make_option(c("-p", "--plot"), action="store_true", default=FALSE,
							help="Create viewpoint coverage plot for all samples."),
	make_option(c("-g", "--genomePlot"), action="store_true", default=FALSE,
							help="Create genomeplot for all samples (only possible if analysis is all in vpFile)."),
	make_option(c("-t", "--tsv"), action="store_true", default=FALSE,
							help="Create tab separated values file for all samples"),
	make_option(c("-b", "--bins"), action="store_true", default=FALSE,
							help="Co-runt reads for binned regions."),
	make_option(c("-chr", "--chr_random"), action="store_true", default=FALSE,
							help="chr_random is in bowtie2 index genome."),
	make_option(c("-chrUn", "--chrUn"), action="store_true", default=FALSE,
							help="chrUn is in bowtie2 index genome."),
	make_option(c("-chrM", "--chrM"), action="store_true", default=FALSE,
							help="chrM is in bowtie2 index genome."),

	make_option(c("-pe", "--peakC"), action="store_true", default=FALSE, 
							help="Perform peakC analysis."),
	make_option(c("-de", "--DESeq2"), action="store_true", default=FALSE, 
							help="Perform DESeq2 analysis on peakC results."),
	make_option(c("-rep", "--replicates"), action="store_true", default=FALSE, 
							help="take in count replicates for peakC and output files.")

)

helptext<-"Values stored in the configuration file (conf.yml) are used by default."

opt_parser = OptionParser(option_list=option_list, description=helptext)
argsL <- parse_args(opt_parser)

if (is.null(argsL$vpFile)){
	print_help(opt_parser)
	stop("vpFile is required.\n", call.=FALSE)
}
if (is.null(argsL$fqFolder)){
	print_help(opt_parser)
	stop("fqFolder is required.\n", call.=FALSE)
}
if (is.null(argsL$outFolder)){
	print_help(opt_parser)
	stop("outFolder is required.\n", call.=FALSE)
}


#################################################################################################################
### Load packages ###############################################################################################
#################################################################################################################
message('\n------ Loading functions and configuration file')

#.libPaths(c(.libPaths(), "/pipeline4C/BSgenomes/installed/"))

if(!suppressMessages(require("ShortRead", character.only=TRUE))) stop("Package not found: ShortRead")
if(!suppressMessages(require("GenomicRanges", character.only=TRUE))) stop("Package not found: GenomicRanges")
if(!suppressMessages(require("GenomicAlignments", character.only=TRUE))) stop("Package not found: GenomicAlignments")
if(!suppressMessages(require("caTools", character.only=TRUE))) stop("Package not found: caTools")
if(!suppressMessages(require("config", character.only=TRUE))) stop("Package not found: config")
if(!suppressMessages(require("rtracklayer", character.only=TRUE))) stop("Package not found: rtracklayer")
if(!suppressMessages(require("dplyr", character.only=TRUE))) stop("Package not found: dyplr")

source(sub(pattern='pipe4C\\.R', replacement='functions\\.R', x=get_script_path()))

if (argsL$confFile=='conf.yml'){
	argsL$confFile<-sub(pattern='pipe4C\\.R', replacement='conf.yml', x=get_script_path()) 
}

configOpt <- createConfig(confFile=argsL$confFile)
configOpt$pipeline.version <- VERSION


if (!is.null(argsL$qualityCutoff)){
	configOpt$qualityCutoff<-argsL$qualityCutoff
}
if (!is.null(argsL$trimLength)){
	configOpt$trimLength<-argsL$trimLength
}
if (!is.null(argsL$minAmountReads)){
	configOpt$minAmountReads<-argsL$minAmountReads
}
if (!is.null(argsL$readsQuality)){
	configOpt$readsQuality<-argsL$readsQuality
}
if (!is.null(argsL$cores)){
	configOpt$cores<-argsL$cores
}
if (!is.null(argsL$wSize)){
	configOpt$wSize<-argsL$wSize
}
if (!is.null(argsL$nTop)){
	configOpt$nTop<-argsL$nTop
}
if (!is.null(argsL$mismatchMax)){
	configOpt$mmMax<-argsL$mismatchMax
}

if (!is.null(argsL$vpRegion)){
	configOpt$vpRegion<-argsL$vpRegion
}
if (!is.null(argsL$alphaFDR)){
	configOpt$alphaFDR<-argsL$alphaFDR
}
if (!is.null(argsL$qWd)){
	configOpt$qWd<-argsL$qWd
}
if (!is.null(argsL$qWr)){
	configOpt$qWr<-argsL$qWr
}
if (!is.null(argsL$minDist)){
	configOpt$minDist<-argsL$minDist
}
if (!is.null(argsL$min.gapwidth)){
	configOpt$min.gapwidth<-argsL$min.gapwidth
}
if (!is.null(argsL$fixedStepWigBin)){
	configOpt$fixedStepWigBin<-argsL$fixedStepWigBin
}


if (argsL$mapUnique){
	configOpt$mapUnique<-argsL$mapUnique
}
if (argsL$nonBlind){
	configOpt$nonBlind<-argsL$nonBlind
}
if (argsL$bdg){
	configOpt$bdg<-argsL$bdg
}
if (argsL$wig){
	configOpt$wig<-argsL$wig
}
if (argsL$fixedStepWig){
	configOpt$fixedStepWig<-argsL$fixedStepWig
}
if (argsL$plot){
	configOpt$cisplot<-argsL$plot
}
if (argsL$genomePlot){
	configOpt$genomePlot<-argsL$genomePlot
}
if (argsL$tsv){
	configOpt$tsv<-argsL$tsv
}
if (argsL$bins){
	configOpt$bins<-argsL$bins
}
if (argsL$chr_random){
	configOpt$chr_random<-argsL$chr_random
}
if (argsL$chrUn){
	configOpt$chrUn<-argsL$chrUn
}
if (argsL$chrM){
	configOpt$chrM<-argsL$chrM
}

if (argsL$peakC){
	configOpt$peakC<-argsL$peakC
}
if (argsL$DESeq2){
	configOpt$DESeq2<-argsL$DESeq2
}
if (argsL$replicates){
	configOpt$replicates<-argsL$replicates
}

#Check whether Bowtie2 and samtools are installed...
message("\n------ Checking the version of bowtie2:\n")
b <- suppressWarnings(system("command -v bowtie2 >/dev/null"))
# if bowtie2 is not installed b = 127
if (b == 0){
	system("bowtie2 --version")
} else {
	stop("Bowtie2 is not installed!\n")	
}

message("\n------ Checking the version of samtools:\n")
s <- suppressWarnings(system("command -v samtools >/dev/null"))
# if Samtools is not installed s = 127
if (s == 0){
	system("samtools --version")
} else {
	stop("samtools is not installed!\n")	
}


#################################################################################################################
### PIPELINE ####################################################################################################
#################################################################################################################

Run.4Cpipeline(VPinfo.file = argsL$vpFile, FASTQ.F = argsL$fqFolder, OUTPUT.F = argsL$outFolder, configuration=configOpt)

#################################################################################################################
### END #########################################################################################################
#################################################################################################################
