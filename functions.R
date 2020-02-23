# FUNCTIONS 4C PIPELINE

createConfig <- function(confFile=argsL$confFile){
	configF <- config::get(file=confFile)
	baseFolder <- configF$fragFolder
	normFactor <- configF$normalizeFactor
	plotView <- configF$plot$plotView
	maxY <- configF$plot$maxY
	xaxisUnit <- configF$plot$xaxisUnit
	plotType <- configF$plot$plotType
	binSize <- configF$plot$binSize

	qualityCutoff <- configF$qualityCutoff
	trimLength <- configF$trimLength
	minAmountReads <- configF$minAmountReads
	readsQuality <- configF$readsQuality
	cores <- configF$cores
	wSize <- configF$wSize
	nTop <- configF$nTop
	mapUnique <- configF$mapUnique
	nonBlind <- configF$nonBlind
	bdg <- configF$bdg
	wig <- configF$wig
	fixedStepWig <- configF$fixedStepWig
	fixedStepWigBin <- configF$fixedStepWigBin
	cisplot <- configF$cisplot
	genomePlot <- configF$genomePlot
	tsv <- configF$tsv
	bins <- configF$bins
	mmMax <- configF$mismatchMax

	peakC <- configF$peakC
	replicates <- configF$replicates

	chr_random <- configF$chr_random
	chrUn <- configF$chrUn
	chrM <- configF$chrM

	vpRegion = configF$PeakC$vpRegion
	alphaFDR = configF$PeakC$alphaFDR
	qWd = configF$PeakC$qWd
	qWr = configF$PeakC$qWr
	minDist = configF$PeakC$minDist
	min.gapwidth = configF$PeakC$min.gapwidth
	DESeq2 = configF$PeakC$DESeq2
	ctrl = configF$PeakC$ctrl
	
	enzymes <- data.frame(
		name=as.character(sapply(configF$enzymes, function(x) strsplit(x, split=' ')[[1]][1])),
		RE.seq=as.character(sapply(configF$enzymes, function(x) strsplit(x, split=' ')[[1]][2])),
		stringsAsFactors=FALSE, row.names=1)

	genomes <- data.frame(
		genome=as.character(sapply(configF$genomes, function(x) strsplit(x, split=' ')[[1]][1])),
		BSgenome=as.character(sapply(configF$genomes, function(x) strsplit(x, split=' ')[[1]][2])),
		stringsAsFactors=FALSE, row.names=1)

	bt2Genomes <- data.frame(
		genome=as.character(sapply(configF$bowtie2, function(x) strsplit(x, split=' ')[[1]][1])),
		path=as.character(sapply(configF$bowtie2, function(x) strsplit(x, split=' ')[[1]][2])),
		stringsAsFactors=FALSE, row.names=1)

	return(list(baseFolder=baseFolder, normFactor=normFactor, enzymes=enzymes, genomes=genomes, bt2Genomes=bt2Genomes, plotView=plotView,
		maxY=maxY, xaxisUnit=xaxisUnit,	plotType=plotType, binSize=binSize,	qualityCutoff=qualityCutoff, trimLength=trimLength, minAmountReads=minAmountReads,
		readsQuality=readsQuality, cores=cores, wSize=wSize, nTop=nTop, mapUnique=mapUnique, nonBlind=nonBlind, bdg=bdg, wig=wig, fixedStepWig=fixedStepWig,
		fixedStepWigBin=fixedStepWigBin, cisplot=cisplot, genomePlot=genomePlot, tsv=tsv, bins=bins, mmMax=mmMax, chr_random=chr_random, chrUn=chrUn, chrM=chrM,
		peakC=peakC, replicates=replicates, vpRegion=vpRegion, alphaFDR=alphaFDR, qWd=qWd, qWr=qWr, minDist=minDist, min.gapwidth=min.gapwidth, DESeq2=DESeq2, ctrl=ctrl))
}

Read.VPinfo <- function(VPinfo.file, replicates){
	#Indentify Fastq files, experiment names and primer sequence from VPinfo file
	VPinfo <- read.table(VPinfo.file, sep="\t", stringsAsFactors=FALSE, header=TRUE)

	if (! "spacer" %in% colnames(VPinfo)){
		VPinfo$spacer<-0
	}
	
	if (ncol(VPinfo)==10){
		defaultNames<-c('expname', 'spacer', 'primer', 'firstenzyme', 'secondenzyme', 'genome', 'vpchr', 'vppos', 'analysis', 'fastq')
		headerCount <- sum(colnames(VPinfo) %in% defaultNames)
		if (replicates == TRUE){
			message("Replicates option activated! Condition and replicate column are required.")
			return()
		}
	}

	if (ncol(VPinfo)==12){
		defaultNames<-c('expname','spacer', 'primer', 'firstenzyme', 'secondenzyme', 'genome', 'vpchr', 'vppos', 'analysis', 'fastq', 'condition', 'replicate')
		headerCount <- sum(colnames(VPinfo) %in% defaultNames) 
	}

	
	if (headerCount < ncol(VPinfo)){
		message("Headers not correct in VPinfo file.")
		return()
	}
	
	# Modification of duplicated exp.name depending of condition and replicate
	if (replicates == TRUE){
		while (any(duplicated(VPinfo$expname)) == TRUE){
			dupName <- VPinfo$expname[which(duplicated(VPinfo$expname) == TRUE)[1]]
			dupVPinfo <- VPinfo[VPinfo$expname == dupName,]
			identicalNameRows <- rownames(dupVPinfo)
			if (length(unique(dupVPinfo$vpchr)) != 1 || length(unique(dupVPinfo$vppos)) != 1){
				message("Some duplicated exp.names do not present the same viewpoint (vpchr or vppos).")
				return()
			}

			for (i in identicalNameRows){
				newName = paste0(dupVPinfo$expname[rownames(dupVPinfo) == i], "_", dupVPinfo$condition[rownames(dupVPinfo) == i], "_Rep", dupVPinfo$replicate[rownames(dupVPinfo) == i])
				dupVPinfo$expname[rownames(dupVPinfo) == i] <- newName
				if (newName %in% unique(dupVPinfo$expname[duplicated(dupVPinfo$expname)])){ 
					message("Some duplicated exp.names present the same condition and replicate.")
					return()
				} else {
					VPinfo$expname[rownames(VPinfo) == i] <- newName
				}
			}
		}
	}

	expname <- as.character(gsub("[^A-Za-z0-9_]", "", VPinfo$expname))
	primer <-   toupper(as.character(gsub("[^A-Za-z]", "", VPinfo$primer)))
	firstenzyme <-  as.character(gsub("[^A-Za-z0-9]", "", VPinfo$firstenzyme))
	secondenzyme <- as.character(gsub("[^A-Za-z0-9]", "", VPinfo$secondenzyme))
	genome <- as.character(gsub("[^A-Za-z0-9_]", "", VPinfo$genome))
	 
	#vpchr <- as.character(gsub("[^0-9XYMxym]", "",VPinfo$vpchr))
	#To be sure that any kind of chromosome name for different species can be used I only remove spaces.
	vpchr <- as.character(gsub(" ", "",VPinfo$vpchr)) 
	
	if(is.numeric(VPinfo$vppos)){
		vppos <- as.numeric(gsub("\\D+", "", VPinfo$vppos))
	} else {
		message("One of the vppos is not numeric.")
		return()
	}
	
	analysis <- as.character(gsub("[^a-z]", "", VPinfo$analysis))
	fastq <- as.character(VPinfo$fastq)
	spacer <- as.numeric(gsub("\\D+", "", VPinfo$spacer))

	if (ncol(VPinfo)==10){
		VPinfo$condition<-NA
		VPinfo$replicate<-1
	}
	condition <- as.character(gsub("[^A-Za-z0-9_]", "", VPinfo$condition))

	if(is.numeric(VPinfo$replicate)){
		replicate <- as.numeric(gsub("\\D+", "", VPinfo$replicate))
	} else {
		message("One of the replicate value is not numeric.")
		return()
	}
	

	return(data.frame(expname=expname, spacer=spacer, primer=primer, firstenzyme=firstenzyme, secondenzyme=secondenzyme, genome=genome, 
		vpchr=vpchr, vppos=vppos, analysis=analysis, fastq=fastq, condition=condition, replicate=replicate,	stringsAsFactors=FALSE))         
}

demux.FASTQ <- function(VPinfo, FASTQ.F, FASTQ.demux.F, demux.log.path, overwrite = FALSE, mmMax){
	# reading functions
	if(!suppressMessages(require("ShortRead", character.only=TRUE))) stop("Package not found: ShortRead")
		
	# Demultiplex FastQ files: 
	# 1.  This will remove weird characters from VPinfo name and primer info.
	# 2.  Exclude experiments that have non-unique names
	# 3.  Demultiplex each Fastq file and save it as a new Fastq file which can be uploaded to GEO.
	# If a spacer sequence is used with random nucleotides before the primer, the spacer option should be set.

	#Demultiplex per FASTQ file

	file.fastq <- sort(unique(VPinfo$fastq))
	for (i in 1:length(file.fastq)){
		VPinfo.fastq <- VPinfo[VPinfo$fastq == file.fastq[i],]

		# Remove weird characters from VPinfo file
		exp.name.all <- as.character(VPinfo$expname)
		exp.name <- as.character(VPinfo.fastq$expname)
		primer <- toupper(as.character(VPinfo.fastq$primer))
		spacer <- as.numeric(VPinfo.fastq$spacer)

		# Remove experiments with duplicated names
		exp.name.unique <- NULL
		primer.unique <- NULL
		spacer.unique <- NULL
		for (a in 1:length(exp.name)){
			if (exp.name[a] %in% exp.name.all[duplicated(exp.name.all)]){
				error.msg <- paste("      ### ERROR: Experiment name not unique for ", exp.name[a])
				message(error.msg)
				write(error.msg, demux.log.path, append = TRUE)
			} else {
				exp.name.unique <- c(exp.name.unique, exp.name[a])
				primer.unique <- c(primer.unique, primer[a])
				spacer.unique <- c(spacer.unique, spacer[a])
			}
		}
		fq.df <- data.frame()
		for (b in 1:length(exp.name.unique)){
			outfile <- paste0(FASTQ.demux.F, exp.name.unique[b], ".fastq.gz")
			if (file.exists(outfile)){
				if (overwrite == TRUE){
					unlink(outfile)
					error.msg <- paste("      ### WARNING: File", outfile, "exists and will be overwritten.")
					write(error.msg, demux.log.path, append = TRUE)
					message(error.msg)
					newrow <- data.frame(exp.name = exp.name.unique[b], primer = primer.unique[b], spacer = spacer.unique[b])
					fq.df <- rbind(fq.df, newrow)
				} else {
					message(paste("      ### WARNING: File", outfile, "exists. continuing with existing file."))
				}
			} else {
				newrow <- data.frame(exp.name = exp.name.unique[b], primer = primer.unique[b], spacer = spacer.unique[b])
				fq.df <- rbind(fq.df, newrow)
			}
		}
		
		#Check whether primers has enough distance with maximum allowed mismatches.
		#
		#https://www.rdocumentation.org/packages/DescTools/versions/0.99.19/topics/StrDist
		#The function computes the Hamming and the Levenshtein (edit) distance of two given strings (sequences).
		#The Hamming distance between two vectors is the number mismatches between corresponding entries.
		#In case of the Hamming distance the two strings must have the same length.
		#

		if (nrow(fq.df)>1){
			message("Check whether primer can be seperated")
			primers.mm<-DNAStringSet(fq.df$primer)
			names(primers.mm)<-fq.df$exp.name

			for (c in seq_along(primers.mm)){
				primers.mm2<-primers.mm[-c]
				#dist <- srdistance(primers.mm2,primers.mm[c])[[1]]
				dist <- srdistance(primers.mm2, primers.mm[c], method = "Hamming")[[1]]

				#srdistanceFilter()
				#http://manuals.bioinformatics.ucr.edu/home/ht-seq

				if (min(dist)<= mmMax){
					error.msg <- paste("      ### WARNING: primer sequence not unique for", fq.df$exp.name[c], "with", mmMax, "mismatches allowed.")
					message(error.msg)
					write(error.msg, demux.log.path, append = TRUE)
					error.msg2 <- paste("      ### WARNING: primer overlaps with primer ", names(primers.mm2[dist<=mmMax]),"\n")
					message(error.msg2)
					write(error.msg2, demux.log.path, append = TRUE)
				}
			}
		}

		# Read FastQ
		if (nrow(fq.df) > 0){
			message(paste("      >>> Reading Fastq: ", file.fastq[i], " <<<"))
			if (file.exists(paste0(FASTQ.F, "/", file.fastq[i]))){

				# then we stream the fastq files at 1,000,000 reads each time (default)
				stream <- FastqStreamer((paste0(FASTQ.F, "/", file.fastq[i])))
				message(paste0("         ### ", fq.df$exp.name, " check for primer sequence ", fq.df$primer, 
					" spacer:", fq.df$spacer, " max mismatch:", mmMax, "\n"))
				while (length(fq <- yield(stream))){
					for (i in 1:nrow(fq.df)){
						primer.seq <- as.character(fq.df$primer[i])
						spacer <- as.numeric(fq.df$spacer[i])
						if (mmMax == 0){
							demultiplex.primer = srFilter(function(x){
								substr(sread(x), spacer + 1, spacer + nchar(primer.seq)) == primer.seq
								}, name = "demultiplex.primer")
							demux.fq <- fq[demultiplex.primer(fq)]
						} else {
							shortreads <- narrow(sread(fq), start = spacer + 1, end = spacer + nchar(primer.seq))
							dist <- srdistance(shortreads, primer.seq, method = "Hamming")[[1]]
							demux.fq <- fq[dist <= mmMax]
						}
						writeFastq(demux.fq, paste0(FASTQ.demux.F, fq.df$exp.name[i], ".fastq.gz"), mode = "a")
					}
				}
				close(stream)
			}else{
				error.msg <- paste("      ### WARNING: File", file.fastq[i], " does not exist. Experiment skipped.")
				write(error.msg, demux.log.path, append = TRUE)
				message(error.msg)
			}
		}
	}
	message("\n")
}

trim.FASTQ <- function(exp.name, firstcutter, secondcutter, file.fastq, trim.F, cutoff, trim.length=0, log.path, min.amount.reads=1000){
	txt.tmp <- paste0(trim.F, exp.name, ".txt")
	info.file <- paste0(trim.F, exp.name, ".info.rds")
	if (file.exists(txt.tmp) & file.exists(info.file)){
		error.msg <- paste0("         ### WARNING: trimmed file ", exp.name, " already exists, continuing with existing file.")
		write(error.msg, log.path, append=TRUE)
		message(error.msg)
		return(readRDS(info.file)) 
	} else {
		message(paste("         ### Reading Fastq: ", file.fastq))
		demux.fq <- readFastq(file.fastq)
		# Qualtity trimming::remove? -- lets test it
		if (cutoff > 0){
			# Trim ends with qualtiy score < cutoff
			# trimTailw remove low-quality reads from the right end using a sliding window
			# trim as soon as 2 of 5 nucleotides has quality encoding less than phred score.
			Phred.cutoff <- rawToChar(as.raw(cutoff + 33))
			demux.fq <- trimTailw(demux.fq, 2, Phred.cutoff, 2)
			# Alternatively use trimTails
			# remove a tally of (successive) nucleotides falling at or below a quality threshold (trimTails).
			# fq.trimTails.10 <- trimTails(fq, k=2, a=Phred.cutoff, successive=FALSE)
			# successive: indicating whether failures can occur anywhere in the sequence, or must be successive
			#Mean average base quality score>cutoff
			#Does this make sense after end trimming? Should the cutoff be higher?
			demux.fq <- demux.fq[(alphabetScore(demux.fq) / width(demux.fq)) > cutoff]
			#Drop reads that are less than 30nt after quality filtering. Discuss lenght with Geert/Valerio
			demux.fq <- demux.fq[width(demux.fq) >= 30]
			#Maybe do this filtering afterwards? read len=Cap length?? (or a few nt less?)
		} #close cutoff

		sequences <- sread(demux.fq)
		nReads <- length(sequences)

		if (nReads < min.amount.reads){
			error.msg <- paste0("         ### ERROR: ",exp.name," - Less reads in FASTQ than set as minimum. Reads: ", nReads)
			write(error.msg, log.path, append=TRUE)
			message(error.msg)
			message("         ### To continue alter minAmountRead argument")
			return()
		} else {
			message(paste0("         ### Total Reads: ", nReads))
		}

		if (trim.length > 0){
			sequences <- substr( sequences, 1, ( trim.length-1+motif.1st.pos ) )
			read.length.table <- sort(table(width(sequences)), decreasing=TRUE)[1]
		} else {
			read.length.table <- sort(table(width(demux.fq)), decreasing=TRUE)[1]
		}

		read.length<-as.numeric(names(read.length.table))
		read.length.perc<-round(as.numeric(read.length.table/nReads*100),2)

		if (read.length.perc < 60){
			error.msg <- paste0("         ### WARNING:", exp.name," - Max read length in only ", read.length.perc, "% of the reads.")
			write(error.msg, log.path, append=TRUE)
			message(error.msg)
		}

		#Find the most occuring position of the firstcutter
		motifPos<-sort(table(regexpr(firstcutter, sequences)), decreasing=TRUE)[1]
		motif.1st.pos<-as.numeric(names(motifPos))
		motifPos.perc<-round(as.numeric(motifPos/nReads*100),2)

		motif.1st.pos.2nd <- FALSE
		if (motif.1st.pos > 0){
			#Check whether firstcutter motif is within the first 4 nts (as part of a barcode). If so take 2nd firstcutter pos.
			if (motif.1st.pos < 5){
				motif.1st.pos.2nd <- TRUE
				#motif.1st.pos <- as.numeric(names(sort(table(gregexpr(firstcutter, sequences)[[1]][2]), decreasing=TRUE))[1])
				
				motifPos<-sort(table(gregexpr(firstcutter, sequences)[[1]][2]), decreasing=TRUE)[1]
				motif.1st.pos<-as.numeric(names(motifPos))
				motifPos.perc<-round(as.numeric(motifPos/nReads*100),2)
			}

			if (motifPos.perc < 90){
				error.msg <- paste0("         ### WARNING: ", exp.name, " - Most occuring position RE1 found in only", motifPos.perc, "% of the reads.")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
			}

			# Trim non-blind fragend sequences based on 2nd cutter 
			motif.2ndRE.pos <- as.numeric(names(sort(table(regexpr(secondcutter, sequences)), decreasing=TRUE))[1]) #Determine whether part of barcode
			if (motif.2ndRE.pos > 5 | motif.2ndRE.pos == -1){
				sequences.no.2nd.enzyme <- sub(paste0(secondcutter, ".*$"), secondcutter, sequences) # Changes only the 1st pattern match per string
			} else { #part of bardcode. use the 2nd motif
				sequences.no.2nd.enzyme <- sub(paste0("(", secondcutter, ".*?)", secondcutter, ".*$"), paste0("\\1", secondcutter), sequences)
			}
			# The . matches any character. .* matches zero or more instances of any character. 
			# But the matching is greedy; it would match as much as possible. I want matching that is not greedy (only up until the second .) 
			# so I added ? to suppress the greedy match and .*? matches any group of characters up until we hit the next thing in the regex.
			# Because the first part was enclosed in parentheses (\\..*?) it is stored as \1,
			# so the substitution pattern \\1 restores everything before the second . and the second . is replaced with the secondcutter.
			# Trim 2nd first cutter motif (Trim blind fragments)
			if (motif.1st.pos.2nd == FALSE){
				sequences.no.2nd.enzyme <-sub(paste0("(", firstcutter, ".*?)", firstcutter, ".*$"), paste0("\\1", firstcutter), sequences.no.2nd.enzyme)  
			}else{
			# Trim 3rd first cutter motif (Trim blind fragments) if present in barcode
			sequences.no.2nd.enzyme <-sub(paste0("(", firstcutter, ".*?", firstcutter, ".*?)", firstcutter, ".*$"), paste0("\\1", firstcutter), sequences.no.2nd.enzyme)
			}
			# Trim sequences based on 1st cutter
			keep <- substr(sequences.no.2nd.enzyme, motif.1st.pos, nchar(sequences.no.2nd.enzyme))
		} # close if (motif.1st.pos>0)

		if(motif.1st.pos == -1){
			rm(sequences, demux.fq)
			error.msg <- paste("         ### ERROR:", exp.name, "1st RE motif not found in reads")
			write(error.msg, log.path, append=TRUE)
			message(error.msg)
			return()
		}
		write(keep, txt.tmp)
		rm(sequences, demux.fq)
		captureLen <- read.length - motif.1st.pos + 1
		saveRDS(list(captureLen=captureLen, nReads=nReads, motifPosperc=motifPos.perc, readlenperc=read.length.perc ), file=info.file)
		return(list(captureLen=captureLen, nReads=nReads, motifPosperc=motifPos.perc, readlenperc=read.length.perc))
	}
}

makeBAM <- function(exp.name, BAM.F, NCORES, Bowtie2IndexFile, txt.tmp, log.path, bowtie.log.path, bamFile, map.unique, readsQual=1){
	TEMPfile <- tempfile(pattern = "aln.", tmpdir = tempdir(), fileext = "")
	check.trunc <- 1
	while (check.trunc < 10){
		message(paste0("         ### attempt ", check.trunc))
		if (map.unique == TRUE){
			bamFile.tmp <- paste0(BAM.F, exp.name, "_tmp.bam")
			CMD <- paste0("(bowtie2 -p ", NCORES, " -r -x ", Bowtie2IndexFile, " -U ", txt.tmp, 
				" | samtools view -hbSu - | samtools sort -T ", TEMPfile, " -o ", bamFile.tmp, ") 2>&1")
			
			bowtie.output <- system(command=CMD, intern=TRUE)
			for(i in bowtie.output){ 
				message(paste0('          ', i)) 
			}

			write(paste("Bowtie2 output:",exp.name), bowtie.log.path, append=TRUE) 
			write(bowtie.output, bowtie.log.path, append=TRUE)   

			#Check whether BAM in truncated (Due to pipe problems?). If truncated try 5 times.
			if (length(grep("truncated file", bowtie.output)) == 1){
				error.msg <- paste("         ### ERROR: ",exp.name," - Truncated BAM file. Repeating Bowtie mapping...")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
				check.trunc <- check.trunc + 1
				unlink(bamFile.tmp)
			} else {
				check.trunc <- 10
				system(paste0("samtools view -H ", bamFile.tmp, " > header.sam")) #Output the header only
				system(paste0("samtools view -F 4 ",bamFile.tmp, " | grep -v \"XS:\" | cat header.sam - | samtools view -b - > ", bamFile ))
				system( paste0( "samtools index ",  bamFile ))
				unlink(bamFile.tmp)
			}
		} else{
			# map with quality, default 1
			CMD <- paste0("(bowtie2 -p ", NCORES, " -r -x ", Bowtie2IndexFile, " -U ", txt.tmp, 
				" | samtools view -q ", readsQual, " -bSu - | samtools sort -T ", TEMPfile, " -o ", bamFile, ") 2>&1")
			bowtie.output <- system(command=CMD, intern=TRUE)
			system( paste0("samtools index ", bamFile))

			for(i in bowtie.output){ 
				message(paste0('          ', i)) 
			}
			write(paste("Bowtie2 output:",exp.name), bowtie.log.path, append=TRUE)
			write(bowtie.output, bowtie.log.path, append=TRUE)
			#Check whether BAM in truncated. If truncated try 10 times.
			if (length(grep("truncated file", bowtie.output)) == 1){
				error.msg <-paste0("         ### ERROR: ",exp.name," - Truncated BAM file. Repeating Bowtie mapping...")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
				check.trunc <- check.trunc + 1
				unlink(bamFile)
			} else{
				check.trunc <- 10
			}
		}
	}
	unlink(TEMPfile)
}

norm4C <- function(readsGR, nReads=1e6, nTop=2, wSize=21){
	readsGR$normReads <- 0
	sumTop <- sum(-sort(-readsGR$reads)[1:nTop])
	wNorm <- nReads/(sum(readsGR$reads)-sumTop)
	readsGR$normReads <- wNorm*readsGR$reads
	readsGR$norm4C <- runmean(x=readsGR$normReads, k=wSize, endrule="mean")
	tmp <- readsGR$norm4C 
	tmp[tmp < 0] <- 0
	readsGR$norm4C <- tmp
	return(readsGR)
}

do.wigToBigWig <- function(wig_file, config_genomes, assemblyName){
	if(!suppressMessages(require(as.character(config_genomes[assemblyName,]), character.only=TRUE))) stop(paste0("Package not found: ", as.character(config_genomes[assemblyName,])))
	if(!suppressMessages(require("BSgenome", character.only=TRUE))) stop("Package not found: BSgenome")

	do.call(require, args=list(config_genomes[assemblyName,]))
	assign('genome', base::get(config_genomes[assemblyName,]))
	info <- seqinfo(genome)
	wigToBigWig(wig_file, info)
}

exportWig <- function(gR, expName, filename, vpPos, vpChr, plotView, config_genomes, assemblyName, fixed=FALSE, bin=50){
	gzname<- paste0(filename, ".gz")
	gz1 <- gzfile(gzname, "w")
	
	browserPos <- paste0(vpChr,":", vpPos-plotView, "-", vpPos+plotView)
	positionLine <- paste0("browser position ", browserPos, "\n")
	cat(positionLine, file=gz1)
	trackLine <- paste0("track type=wiggle_0 name=", expName, " visibility=full autoScale=off  viewLimits=0.0:2500.0 maxHeightPixels=50:50:8 color=0,0,200  priority=10\n")
	cat(trackLine, file=gz1, append=TRUE)
	if (fixed){
		data <- GRanges(seqnames=seqnames(gR), ranges=ranges(gR), strand = strand(gR), score=gR$norm4C)
		#Remove the last tile of each chromosome if its length is not equal to the bin size 
		data <- data[-which(end(data) - start(data) != bin-1)]
	} else {
		data <- GRanges(seqnames=seqnames(gR), ranges=gR$pos, strand = strand(gR), score=round(gR$norm4C, digits=5))
	}
	export.wig(data, gz1, append = TRUE)
	close(gz1)
	
	do.wigToBigWig(wig_file=gzname, config_genomes=config_genomes, assemblyName=assemblyName)
}

exportBDG <- function(gR, expName, filename, vpPos, vpChr, plotView, interval=FALSE){
	gzname<- paste0(filename, ".gz")
	gz1 <- gzfile(gzname, "w")
	positionLine <- paste0("browser position ", vpChr,":", vpPos-plotView, "-", vpPos+plotView, "\n")
	cat(positionLine, file=gz1)
	trackLine <- paste0("track type=bedGraph name=", expName, " visibility=full autoScale=off  viewLimits=0.0:2500.0 maxHeightPixels=50:50:8 color=0,0,200  priority=10\n")
	cat(trackLine, file=gz1, append=TRUE)
	
	if (interval == TRUE){
		data_bdg <- data.frame(seqnames=seqnames(gR), start=start(gR), end=end(gR), score4C=round(gR$norm4C, digits=5))
	} else {
		data_bdg <- data.frame(seqnames=seqnames(gR), start=gR$pos, end=gR$pos, score4C=round(gR$norm4C, digits=5))
	}	
	write.table(data_bdg, file=gz1, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	close(gz1)
}

exportTSV <- function(gR, filename, merge=FALSE){
	gzname<- paste0(filename, ".gz")
	gz1 <- gzfile(gzname, "w")
	if (merge) {
		tsvdat <- data.frame(chr=seqnames(gR), start=start(gR), end=end(gR), pos=gR$pos, type=gR$type, fe_strand=gR$fe_strand,
			fe_id=gR$fe_id, normReads=gR$normReads, norm4C=gR$norm4C)
	} else {
		tsvdat <- data.frame(chr=seqnames(gR), start=start(gR), end=end(gR), pos=gR$pos, type=gR$type, fe_strand=gR$fe_strand,
			fe_id=gR$fe_id, reads=gR$reads, normReads=gR$normReads, norm4C=gR$norm4C)
	}
	write.table(tsvdat, file=gz1, append=FALSE, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	close(gz1)
}

olRanges <- function(query, subject){
	## Find overlapping ranges
	olindex <- as.matrix(findOverlaps(query, subject))
	query <- query[olindex[,1]]
	subject <- subject[olindex[,2]]
	olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))

	## Pre-queries for overlaps
	startup <- olma[,"Sstart"] < olma[,"Qstart"]
	enddown <- olma[,"Send"] > olma[,"Qend"]
	startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
	endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]

	## Overlap types
	olup <- startup & endin
	oldown <- startin & enddown
	inside <- startin & endin
	contained <- startup & enddown
	
	## Overlap types in one vector
	OLtype <- rep("", length(olma[,"Qstart"]))
	OLtype[olup] <- "olup"
	OLtype[oldown] <- "oldown"
	OLtype[inside] <- "inside"
	OLtype[contained] <- "contained"
	
	## Overlap positions
	OLstart <- rep(0, length(olma[,"Qstart"]))
	OLend <- rep(0, length(olma[,"Qstart"]))
	OLstart[olup] <- olma[,"Qstart"][olup]
	OLend[olup] <- olma[,"Send"][olup]
	OLstart[oldown] <- olma[,"Sstart"][oldown]
	OLend[oldown] <- olma[,"Qend"][oldown]
	OLstart[inside] <- olma[,"Sstart"][inside]
	OLend[inside] <- olma[,"Send"][inside]
	OLstart[contained] <- olma[,"Qstart"][contained]
	OLend[contained] <- olma[,"Qend"][contained]

	## Absolute and relative length of overlaps
	OLlength <- (OLend - OLstart) + 1
	OLpercQ <- OLlength/width(query)*100
	OLpercS <- OLlength/width(subject)*100

	## Output type
	oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
	elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf[,-c(3,4)])
	return(query)
}

alignToFragends <- function(gAlign, fragments, firstcut){
	strand(fragments) <- ifelse(fragments$fe_strand==3, "-", "+")
	ovl <- olRanges(query=fragments, subject=as(gAlign, "GRanges"))
	#remove sequences that overlap only with RE motif
	ovl <- ovl[ovl$OLlength >= nchar(as.character(firstcut))+1]
	#Sequence need to start within unique fragment
	ovl <- ovl[strand(ovl)=="+" & ovl$Sstart >= start(ovl) | strand(ovl)=="-" & ovl$Send <= end(ovl)]

	#To do: Make it more strict. Reads should start directly at a restriction enzyme cutting site. 
	#ovl <- ovl[strand(ovl)=="+" & ovl$Sstart >= start(ovl) & ovl$Sstart <= start(ovl)+nchar(as.character(firstcut))| strand(ovl)=="-" & ovl$Send <= end(ovl) & & ovl$Send >= start(ovl)-nchar(as.character(firstcut))]
	
	nReads <- tapply(ovl$Sindex, ovl$Qindex, length)
	fragments$reads <- 0
	fragments$reads[as.numeric(names(nReads))] <- nReads
	return(fragments)
}

Digest <- function(assemblyName, firstcutter_Digest, secondcutter_Digest, baseFolder_Digest, config_genomes, chr_random, chrUn, chrM){
	if(!suppressMessages(require(as.character(config_genomes[assemblyName,]), character.only=TRUE))) stop(paste0("Package not found: ", as.character(config_genomes[assemblyName,])))
	if(!suppressMessages(require("BSgenome", character.only=TRUE))) stop("Package not found: BSgenome")

	firstcutter <- as.character(firstcutter_Digest)
	secondcutter <- as.character(secondcutter_Digest)

	rdsFile <- paste0(baseFolder_Digest, assemblyName, "_", firstcutter, "_", secondcutter, ".rds")

	if (!file.exists(rdsFile)){
		do.call(require, args=list(config_genomes[assemblyName,]))
		assign('frag.genome', base::get(config_genomes[assemblyName,]))

		chr <- unique(seqnames(frag.genome))

		#The "hap" ones are alternate assemblies for certain regions.DO NOT USE THE *hap* files !!!!
		chr <- chr[grep(pattern="_hap", x=chr, invert=TRUE)]
		chr <- chr[grep(pattern="_alt", x=chr, invert=TRUE)]

		#Make sure Bowtie2 index does not contain these chr
		if (chr_random==FALSE){
			chr <- chr[grep(pattern="_random", x=chr, invert=TRUE)]
		}
		if (chrUn==FALSE){
			chr <- chr[grep(pattern="chrUn_", x=chr, invert=TRUE)]
		}
		if (chrM==FALSE){
			chr <- chr[grep(pattern="chrM", x=chr, invert=TRUE)]
		}

		outFrags <- GRanges()
		for (i in seq_along(chr)){
			chrom <- chr[i]
			print(paste("Processing ", chrom, sep=""))

			RE1 <- GRanges(seqnames=chrom, ranges(matchPattern(pattern=firstcutter, subject=frag.genome[[chrom]])))
			RE2 <- GRanges(seqnames=chrom, ranges(matchPattern(pattern=secondcutter, subject=frag.genome[[chrom]])))

			frag.RE1 <- RE1[1:(length(RE1)-1)]
			end(frag.RE1)[1:(length(RE1)-1)] <- end(RE1)[2:length(RE1)]

			# Only keep RE2 sites that completely overlap with RE1 fragment
			RE2.ol <- olRanges(frag.RE1, RE2)
			RE2.ol <- RE2[RE2.ol[RE2.ol$OLpercS==100]$Sindex]

			# blind = all RE1 frg which do not contain RE2 site
			blinds <- frag.RE1[!(frag.RE1 %over% RE2.ol)]
			blinds$type <- "blind"
			blinds <- rep(blinds[blinds$type=="blind"], each = 2)
			blinds$fe_strand <- rep(c(5,3), length.out = length(blinds))

			# nonBlind = all RE1 frg which contain RE2 site
			nonBlinds <- frag.RE1[unique(findOverlaps(frag.RE1,RE2.ol)@from)]
			nonBlinds$type <- "non_blind"

			# 
			nonBlinds.fe5 <- nonBlinds
			end(nonBlinds.fe5) <- end(RE2.ol[findOverlaps(nonBlinds, RE2.ol, select="first")])
			nonBlinds.fe5$fe_strand <- 5

			# 
			nonBlinds.fe3 <- nonBlinds
			start(nonBlinds.fe3) <- start(RE2.ol[findOverlaps(nonBlinds, RE2.ol, select="last")])
			nonBlinds.fe3$fe_strand <- 3

			nonBlinds.fe5.start <- GRanges()
			nonBlinds.fe3.end <- GRanges()

			#Add first and last fragends
			if (start(head(RE1, 1)) > start(head(RE2, 1))){
				frag <- GRanges(seqnames=chrom, IRanges(1, end(RE1[1])))
				RE2.ol.start <- olRanges(frag, RE2)
				RE2.ol.start <- RE2[RE2.ol.start[RE2.ol.start$OLpercS==100]$Sindex]
				nonBlinds.fe5.start <- RE1[1]
				start(nonBlinds.fe5.start) <- start(RE2.ol.start[findOverlaps(frag, RE2.ol.start,select="last")])
				nonBlinds.fe5.start$type <- "non_blind"
				nonBlinds.fe5.start$fe_strand <- 5
			}

			if (end(tail(RE1, 1)) < end(tail(RE2, 1))){
				frag <- GRanges(seqnames=chrom, IRanges(start(RE1[length(RE1)]), end(RE2[length(RE2)])))
				RE2.ol.start <- olRanges(frag, RE2)
				RE2.ol.start <- RE2[RE2.ol.start[RE2.ol.start$OLpercS==100]$Sindex]
				nonBlinds.fe3.end <- RE1[length(RE1)]
				end(nonBlinds.fe3.end) <- end(RE2.ol.start[findOverlaps(frag, RE2.ol.start,select="first")])
				nonBlinds.fe3.end$type <- "non_blind"
				nonBlinds.fe3.end$fe_strand <- 3
			}

			outFrags <- sort(c(outFrags,blinds, nonBlinds.fe5,nonBlinds.fe3,nonBlinds.fe5.start,nonBlinds.fe3.end))
		}

		outFrags$fe_id <- 1:length(outFrags)
		outFrags$pos <- ifelse(outFrags$fe_strand == 5, start(outFrags)+nchar(firstcutter), end(outFrags)-nchar(firstcutter))
		#outFrags$posEnd <- ifelse(outFrags$fe_strand == 3, start(outFrags)+nchar(firstcutter), end(outFrags)-nchar(firstcutter))

		#Save new RDS
		saveRDS(outFrags, file=rdsFile)
		return(outFrags)
	} else {
		message('Frag genome file already exists.')
		return(0)
	}
}

getUniqueFragends <- function(fragsGR_Unique, firstcutter_Unique, secondcutter_Unique, captureLen_Unique=50, nThreads_Unique=4, baseFolder_Unique, 
	Bowtie2IndexFile, assemblyName, config_genomes){
	if(!suppressMessages(require(as.character(config_genomes[assemblyName,]), character.only=TRUE))) stop(paste0("Package not found: ", as.character(config_genomes[assemblyName,])))
	if(!suppressMessages(require("BSgenome", character.only=TRUE))) stop("Package not found: BSgenome")

	##If captureLen is not in colnames, generate this column
	if (paste0("len", captureLen_Unique) %in% colnames(elementMetadata(fragsGR_Unique)) == FALSE){
		message(paste('          ### Create FASTA file for capture length', captureLen_Unique))
		fragends <- fragsGR_Unique
		fragFile <- paste0(baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, ".rds")
		fastaFile <- paste0(baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, "_", captureLen_Unique, ".fa")
		bamFile <- paste0(baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, "_", captureLen_Unique, ".bam")
		uniquesFile <- paste0(baseFolder_Unique, assemblyName, "_", firstcutter_Unique, "_", secondcutter_Unique, "_", captureLen_Unique, "_unique.txt")

		## Create FASTA files for mapping back to REFGENOMEs

		#Extract sequence 3'ends
		start(fragends)[fragends$fe_strand == 3 & end(fragends) > (start(fragends)+captureLen_Unique - 1)]<-
		end(fragends)[fragends$fe_strand == 3 & end(fragends) > (start(fragends)+captureLen_Unique - 1)] - captureLen_Unique + 1

		#Extract sequence 5'ends
		end(fragends)[fragends$fe_strand == 5 & start(fragends) < (end(fragends)-captureLen_Unique + 1)]<-
		start(fragends)[fragends$fe_strand == 5 & start(fragends) < (end(fragends)-captureLen_Unique + 1)] + captureLen_Unique - 1

		do.call(require, args=list(config_genomes[assemblyName,]))
		seqs <- getSeq(x=base::get(config_genomes[assemblyName,]), names=fragends)
		names(seqs) <- fragends$fe_id
		writeXStringSet(seqs,fastaFile)

		message(paste('          ### Mapping FASTA file back to',assemblyName))

		#cmd <- paste0("bowtie2 -p ",nThreads," -x ",Bowtie2IndexFile," -f ",fastaFile," | samtools view -q 1 -hbSu - | samtools sort -T /tmp - > ",bamFile)
		#-q 1: Skip alignments with MAPQ smaller than 1
		#-h Include the header in the output.
		#-b Output in BAM format --> why do we need to output as a BAM file? This is just a temp file.
		#-S Input is in SAM format

		#-u Output uncompressed BAM. 
		#This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command. 
		#No sorting is quicker

		cmd <- paste0("(bowtie2 -p ", nThreads_Unique, " -x ", Bowtie2IndexFile, " -f ", fastaFile, " | samtools view -q 1 -hbSu - > ", bamFile, ") 2>&1")

		check.trunc <- 0
		while (check.trunc < 1){
			message(paste("          ### Bowtie2:", captureLen_Unique))
			bowtie.output <- system(cmd, intern=TRUE)
			if (length(grep("[main_samview] truncated file", bowtie.output)) == 1){
				message(paste("          ### ERROR: Truncated BAM file: repeating Bowtie mapping"))
			} else { 
				check.trunc <- 1
			}
		}

		system(paste0("samtools view ", bamFile, " | grep -v \"XS:\" | cut -f 1 > ", uniquesFile))
		ids <- as.numeric(readLines(uniquesFile))
		elementMetadata(fragsGR_Unique)[[paste0("len", captureLen_Unique)]] <- fragsGR_Unique$fe_id%in%ids #Use Fragend coordinates. not coordinates for unique mapping. 
		fragsGR_Unique <- fragsGR_Unique[, order(colnames(elementMetadata(fragsGR_Unique)))]
		saveRDS(object=fragsGR_Unique, file=fragFile)
		unlink(c(fastaFile, bamFile, uniquesFile))
	}
	# if the unique capture len is already in fragGR, we just use that column and rename it to 'unique'
	fragsGR2 <- fragsGR_Unique[, c("pos", "type", "fe_strand", "fe_id", paste0("len", captureLen_Unique))]
	colnames(elementMetadata(fragsGR2))[5] <- "unique"
	return(fragsGR2)
}

getFragMap <- function(vpChr_FragMap=NULL, firstcutter_FragMap, secondcutter_FragMap, genome_FragMap, captureLen_FragMap=60, nThreads_FragMap=10, 
	baseFolder_FragMap, Bowtie2IndexFile, config_genomes, chr_random=chr_random, chrUn=chrUn, chrM = chrM){
	# here we want to check if we have the genome available for this analysis, in case yes, loading it, otherwise, advice on the 
	# missing genome and switch to the next guy
	# message(paste("      >>> Create frag map for genome :",genome[i],"with re1 :",firstcutter," and re2 :",secondcutter,"and captureLen", captureLen, "<<<"))
	fragFile <- paste0(baseFolder_FragMap, genome_FragMap,"_", firstcutter_FragMap, "_", secondcutter_FragMap, ".rds")

	if (file.exists(fragFile)){
		# File with unique fragends exists, use it to check for unique fragends for given capture length 
		message('         ### Loading: ', fragFile)
		fragsGR <- readRDS(fragFile) 
	} else {
		# Make the digest first
		message(paste0("         ### Creating new digest for genome: ", genome_FragMap, " with re1: ", firstcutter_FragMap, " and re2: ", secondcutter_FragMap))
		fragsGR <- Digest(assemblyName=genome_FragMap, firstcutter_Digest=firstcutter_FragMap, secondcutter_Digest=secondcutter_FragMap, 
			baseFolder_Digest=baseFolder_FragMap, config_genomes=config_genomes, chr_random=chr_random, chrUn=chrUn, chrM = chrM) 
	}

	message('         ### Retrieve and store the unique fragends')
	fragsGR <- getUniqueFragends(fragsGR_Unique=fragsGR, firstcutter_Unique=firstcutter_FragMap, secondcutter_Unique=secondcutter_FragMap,
		captureLen_Unique=captureLen_FragMap, nThreads_Unique=nThreads_FragMap, baseFolder_Unique=baseFolder_FragMap, Bowtie2IndexFile=Bowtie2IndexFile,
		assemblyName=genome_FragMap, config_genomes=config_genomes)

	if (!is.null(vpChr_FragMap)){
		message('         ### Selecting fragments from ', vpChr_FragMap)
		fragsGR <- subset(fragsGR, as.character(seqnames(fragsGR)) == vpChr_FragMap & unique == TRUE)
	} else {
		message('         ### Selecting fragments from whole genome', vpChr_FragMap)
		fragsGR <- fragsGR[fragsGR$unique]
	}
	return(fragsGR)
}

plot.chroms <- function(exp, reads, cutoff=0.999, yMax=2500){
	chroms <- as.character(unique(seqnames(reads)))
	total.reads <- sum(reads$norm4C)
	abs.cut <- quantile(reads$norm4C, cutoff, na.rm=T)
	reads$norm4C[reads$norm4C > abs.cut] <- abs.cut
	reads$norm4C <- reads$norm4C/abs.cut
	reads$chr <- sub("chr","",as.character(seqnames(reads)))
	reads[seqnames(reads)=="chrX"]$chr<-length(chroms)-1
	reads[seqnames(reads)=="chrY"]$chr<-length(chroms)
	reads$chr <-as.numeric(reads$chr)
	layout(cbind(c(rep(1,3),2)))
	par(mar=c(5, 4, 4, 2))
	plot(c(0,max(reads$pos)), c(0,length(chroms)), type='n', axes=F, xlab="Chromosome position (Mb)", ylab="", main = exp)
	segments(reads$pos, reads$chr, reads$pos, reads$chr + reads$norm4C, lwd=0.5)
	lab <- seq(0,ceiling(max(reads$pos)/10e6)*10, by=10)
	axis(1, at = lab*1e6, label=lab, las=2)
	axis(2, at= 1:length(chroms)+0.5, lab=chroms, lwd=0, las=2)
}

make.reads.and.bins <- function(reads, assemblyName, res=25e3, config_genomes, spacer=" "){

	if(!suppressMessages(require(as.character(config_genomes[assemblyName,]), character.only=TRUE))) stop(paste0("Package not found: ", as.character(config_genomes[assemblyName,])))
	if(!suppressMessages(require("BSgenome", character.only=TRUE))) stop("Package not found: BSgenome")
	
	#Extract fragends with reads
	Captures <- reads[reads$norm4C>0]
	Captures.rds<-Captures[,1]
	Captures.rds$normReads<-round(Captures$normReads,3)
	Captures.rds$norm4C<-round(Captures$norm4C,3)

	#make bins
	#Get genome info
	do.call(require, args=list(config_genomes[assemblyName,]))
	assign('genome', base::get(config_genomes[assemblyName,]))

	#Chr size
	chr <- as.character(unique(seqnames(reads)))
	x<-seqinfo(genome)
	seqlevels(x) <- chr

	spacer<- paste0(spacer, "  ")
	message(paste0(spacer, "making bins"))
	bin.GR <- tileGenome(x, tilewidth=res, cut.last.tile.in.chrom=TRUE)
	bin.GR$pos <- start(resize(bin.GR,width=1,fix="center"))

	#overlap
	message(paste0(spacer, "Overlapping reads with bins"))
	reads <- Captures.rds
	hits <- findOverlaps(reads,bin.GR)
	reads$bin <- 0
	reads[queryHits(hits),]$bin <- subjectHits(hits)
	
	sum.reads <- aggregate(reads$normReads, by = list(reads$bin), FUN=sum)
	bin.GR$normReads<-0
	bin.GR$normReads[sum.reads$Group.1] <- sum.reads$x

	sum.reads2 <- aggregate(reads$norm4C, by = list(reads$bin), FUN=sum)
	bin.GR$norm4C<-0
	bin.GR$norm4C[sum.reads2$Group.1] <- sum.reads2$x

	return(list(reads=Captures.rds,bins=bin.GR))
}

doBins <- function(file, expname, bin, reads, assemblyName, config_genomes, report, vpInfo, log.path){
	if (nrow(vpInfo) != 1){
		spacer = "                 "
		msg <- paste0(spacer, "--- ")
	} else {
		spacer = "           "
		msg <- paste0(spacer, "### ")
	}

	if (file.exists(file)){
		if (bin < 1e3){
			error.msg <- paste0(msg, "WARNING: rds bin file ", expname, " already exist with a bin of ", bin, ", continuing with exisiting file.")
		} else {
			error.msg <- paste0(msg, "WARNING: rds bin file ", expname, " already exist with a bin of ", (bin/1e3), "kb, continuing with exisiting file.")
		}
	
		write(error.msg, log.path, append=TRUE)
		message(error.msg)
		rds <- readRDS(file)
		return(list(reads=rds$reads,bins=rds$bins))
	} else {
		message(paste0(msg, "counting normalized reads on bins."))
		bin.GR <- make.reads.and.bins(reads=reads, assemblyName=assemblyName, res=bin, config_genomes=config_genomes, spacer=spacer)
		saveRDS(object=list(reads=bin.GR$reads, bins=bin.GR$bins, report=report, vpInfo=vpInfo), file=file, compress=TRUE)
		return(bin.GR)
	}
}

export.report <- function(RDS.F, OUTPUT.F, logname){ 
	files <- list.files(path=RDS.F, pattern=".rds", recursive=FALSE)
	report.df <- data.frame()
	for(i in 1:length(files)){    
		vpname <- gsub("[.].*$", "", files[i])
		rds <- readRDS(paste0(RDS.F, files[i]))
		report <- rds$report
		newrow <- data.frame(vpname=vpname,report)
		report.df <- rbind(report.df, newrow)
	}
	write.table(report.df, file=paste0(OUTPUT.F,"/report_", logname,".txt"), row.names=FALSE, quote=FALSE)
}

createReport <- function(allReads, mapReads, demuxReads, chromosome, vpPos, normFactor=1e6, wSize, nTop, motifPosperc, readlenperc){
	nMapped <- length(mapReads)
	uniqueReads <- sum(allReads$reads)
	uniqueCaptures <- sum(allReads$reads>0)

	reads <- allReads[as.vector(seqnames(allReads)) == chromosome]

	cisReads <- sum(reads$reads)
	percCis<-round(100*cisReads/uniqueReads, 2) #this does not exclude the Top reads..
	cisCaptures <- sum(reads$reads>0)

	nMappedCis <- length(mapReads[seqnames(mapReads) == chromosome])
	nMappedCisperc <- round(100*nMappedCis/nMapped, 2)

	topIdx <- 1:length(reads) %in% order(-reads$reads)[1:nTop]
	topReads <- sum(reads$reads[topIdx])
	topPct <- round(100*topReads/cisReads, digits=2)
	readsGR <- reads[!topIdx]

	allReads.nTop <- subsetByOverlaps(allReads, reads[topIdx], invert=T)
	uniqueReads.nTop <- sum(allReads.nTop$reads)
	cisReads.nTop <- sum(readsGR$reads)
	percCis.nTop <- round(100*cisReads.nTop/uniqueReads.nTop, 2)

	vpWithin100Kb <- readsGR[unique(queryHits(findOverlaps(ranges(readsGR), resize(IRanges(vpPos, vpPos), width=2e5, fix="center"))))]
	capt100Kb <- round(100*mean(vpWithin100Kb$reads>0), digits=2)
	cov100Kb <- round(100*sum(vpWithin100Kb$reads)/sum(readsGR$reads), digits=2)

	vpWithin1Mb <- readsGR[unique(queryHits(findOverlaps(ranges(readsGR), resize(IRanges(vpPos, vpPos), width=2e6, fix="center"))))]
	capt1Mb <- round(100*mean(vpWithin1Mb$reads>0), digits=2)
	cov1Mb <- round(100*sum(vpWithin1Mb$reads)/sum(readsGR$reads), digits=2)

	report <- data.frame(
		nReads=demuxReads, # total demux reads
		motifPosperc=motifPosperc,
		readlenperc=readlenperc,
		nMapped=nMapped, # Bowtie mapped read
		nMappedCis=nMappedCis, # Bowtie mapped reads in Cis
		nMappedCisperc=nMappedCisperc, # Bowtie mapped % reads in Cis
		fragMapped=uniqueReads, # Toal unique reads
		fragMappedCis=cisReads,
		fragMappedCisPerc=percCis,
		fragMappedCisCorr=cisReads.nTop,
		fragMappedCisPercCorr=percCis.nTop,
		nCaptures=uniqueCaptures,
		nCisCaptures=cisCaptures,
		topPct=topPct,
		capt100Kb=capt100Kb,
		cov100Kb=cov100Kb,
		capt1Mb=capt1Mb,
		cov1Mb=cov1Mb,
		stringsAsFactors=FALSE
	)

	return(report)
}

createPlot <- function(plotTitle, vpPos, chromosome, fragGR, plotLegend=NULL, plotView=1e6, maxY=2500, minY=0, xaxisUnit=c('Mb', 'Kb', 'bp'), 
	plotRegion='cis', foldOut='./', plotType=c('PDF', 'PNG')){
	xaxisUnit <- match.arg(xaxisUnit, c('Mb', 'Kb', 'bp'))
	if(plotRegion == 'cis'){
		if(xaxisUnit == 'Mb'){
			scaleValue <- 1e6
		} else if (xaxisUnit == 'Kb'){
				scaleValue <- 1e3
		} else {
			scaleValue <- 1
		}
		zoom <- GRanges(seqnames=chromosome, IRanges(start=vpPos-plotView, end=vpPos+plotView))
		if (start(zoom) < 1){
			start(zoom) <- 1
		}
		reads.zoom <- fragGR[unique(findOverlaps(fragGR, zoom)@from)]
		if(plotType == 'PDF'){
			pdf(file=paste0(foldOut, plotTitle, ".pdf"))
		} else {
			png(file=paste0(foldOut, plotTitle, ".png"), width=3200, height=3200, res=300)
		}
		plot(reads.zoom$pos/scaleValue, reads.zoom$norm4C,
			type = "h",
			xlim = c(start(zoom)/scaleValue, end(zoom)/scaleValue),
			ylim = c(minY, maxY),
			xlab = paste0(chromosome," - Chromosome position (", xaxisUnit, ")"),
			ylab = paste0("4C Coverage / ", scaleValue, " mapped reads"),
			frame.plot = FALSE,
			main = plotTitle)
		if (!is.null(plotLegend)){
			legend(x="topleft", bty="n", legend=paste0(names(plotLegend), " : ", plotLegend))
		}
		dev.off()
	} else if(plotRegion == "all"){
		if(plotType == 'PDF'){
			pdf(file=paste0(foldOut, plotTitle, ".pdf"))
		} else {
			png(file=paste0(foldOut, plotTitle, ".png"), width = 1200, height = 600, res = 300)
		}
		plot.chroms(exp=plotTitle, reads=fragGR$bins, cutoff=0.999, yMax=maxY)
		dev.off()
	}
}

getVPReads <- function(rds, vpRegion=2e6){
	reads <- rds$reads
	vppos <- rds$vpInfo$pos
	vpChr <- rds$vpInfo$chr

	zoom <- GRanges(seqnames=vpChr, resize(IRanges(vppos, vppos), width=vpRegion, fix="center"))
	vpGR <- reads[unique(queryHits(findOverlaps(reads,zoom)))]
	peakCDat <- data.frame(pos=vpGR$pos, reads=vpGR$normReads)

	return(peakCDat)
}

getPeakCPeaks <- function(resPeakC, min.gapwidth=0){

	vpChr <- resPeakC$vpChr

	peakRanges <- reduce(IRanges(resPeakC$peak, resPeakC$peak), min.gapwidth=min.gapwidth)
	peakGR <- GRanges(seqnames=vpChr, ranges=peakRanges)

	return(peakGR)
}

exportPeakCPeaks <- function(resPeakC, bedFile, name=NULL, desc=NULL, includeVP=TRUE, min.gapwidth=0){

	vpChr <- resPeakC$vpChr
	vpPos <- resPeakC$vpPos

	if(is.null(name)){
		name <- "peakC_track"
	}

	if(is.null(desc)){
		desc <- paste0("peakC peaks on ", vpChr)
	}	

	browserPos <- paste0(vpChr, ":", vpPos-1e6,"-", vpPos+1e6)
	browserPosLine <- paste("browser position", browserPos, "\n")
	cat(browserPosLine, file=bedFile)

	trackLine <- paste0('track name=\"', name, '\" description=\"', desc, '\" visibility=2 itemRgb=\"On\"', "\n")
	cat(trackLine, file=bedFile, append=TRUE)

	if(min.gapwidth>0){
		resPeakC$exportPeakGR <- getPeakCPeaks(resPeakC, min.gapwidth=min.gapwidth)
	}

	if(includeVP){
		bedDF <- as.data.frame(resPeakC$exportPeakGR)[,1:3]
		colnames(bedDF)[1] <- "chr"
		write.table(bedDF, file=bedFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
	} else {
		vpGR <- resize(GRanges(seqnames=vpChr, IRanges(vpPos,vpPos)), width=1e3, fix="center")
		bedDF <- as.data.frame(resPeakC$exportPeakGR)[,1:3]
		colnames(bedDF)[1] <- "chr"
		write.table(bedDF, file=bedFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
	}

	return(bedDF)
}

doPeakC <- function(rdsFiles, vpRegion=2e6, wSize=21, alphaFDR=0.05, qWd=1.5, qWr=1, minDist=15e3){
	if(!suppressMessages(require("GenomicRanges", character.only=TRUE))) stop("Package not found: GenomicRanges")
	if(!suppressMessages(require("peakC", character.only=TRUE))) stop("Package not found: peakC")
	if(!suppressMessages(require("isotone", character.only=TRUE))) stop("Package not found: isotone")

	# Analyze num.exp replicated 4C-Seq experiments with peakC 
	num.exp <- length(rdsFiles)
	message(paste0("         ### Performing peakC on ",num.exp, " experiments:"))

	if(num.exp>1){
		peakCDat <- list()
		vppos <- vector()
		for(i in 1:num.exp){
			if(file.exists(rdsFiles[i])){
				message("          Loading data for experiment: ", basename(rdsFiles[i]))
				rds <- readRDS(rdsFiles[i])
				vppos[i] <- rds$vpInfo$pos
				vpChr <- as.vector(rds$vpInfo$chr)
				peakCDat[[i]] <- getVPReads(rds=rds, vpRegion=vpRegion)
			} else {
				stop(paste0("File not found: ",rdsFiles[i]))
			}
		}
		vppos <- unique(vppos)
		vpChr <- unique(vpChr)

		message("          viewpoint position: ", vpChr, ":", vppos)

			if(length(vppos)==1){
				resPeakC <- suppressWarnings(combined.analysis(data=peakCDat, num.exp=num.exp, vp.pos=vppos, wSize=wSize, 
					alphaFDR=alphaFDR, qWr=qWr, minDist=minDist))
			} else {
				stop("          Multiple viewpoints have been found in your rds files !")
			}
		} else if(num.exp==1){
			message("          Loading data for experiment: ", basename(rdsFiles[1]))
			rds <- readRDS(rdsFiles[1])
			vppos <- rds$vpInfo$pos
			vpChr <- as.vector(rds$vpInfo$chr)
			message("          viewpoint position: ", vpChr, ":", vppos)
			peakCDat <- getVPReads(rds=rds,vpRegion=vpRegion)
			resPeakC <- suppressWarnings(single.analysis(data=peakCDat,vp.pos=vppos,wSize=wSize,qWd=qWd,qWr=qWr,minDist=minDist))
		} else {
			stop("No rds file!")
		}

		resPeakC$vpPos <- vppos
		resPeakC$vpChr <- vpChr
		if(length(resPeakC$peak)>0){
			resPeakC$exportPeakGR <- getPeakCPeaks(resPeakC, min.gapwidth=0)
		} else {
			resPeakC$exportPeakGR <- NULL
		}
		return(resPeakC)
}

readConditions <- function(vpInfo, log.path){
	cond <- unique(vpInfo$condition)
	to_process <- list()
	for (k in cond){
		condVPinfo <- vpInfo[vpInfo$condition == k,]
		condVPinfo$VP <- paste0(condVPinfo$vpchr, "-", condVPinfo$vppos)
		name_list <- list()
		for (z in as.character(unique(condVPinfo$VP))) {
			tmp <- condVPinfo[condVPinfo$VP == z,]
			names <- as.vector(tmp$expname)
			if (length(names) == 1){
				next
			} else if(length(names) == length(unique(tmp$replicate))){
				name_list[[z]] <- names
			} else {
				msg <-""
				for (t in names) {
					msg <- paste0(msg, " ", t)
				}
				error.msg <- paste0("         ### WARNING: these experiments:", msg, "present identical viewpoint, condition and replicate.")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
			}				
		}
		to_process[[k]] <- name_list
	}
	return(to_process)
}

testRDS <- function(rdsFile, parameters, expname, analysis){
	if(!suppressMessages(require("dplyr", character.only=TRUE))) stop("Package not found: dyplr")
	rds <- readRDS(rdsFile)
	if (is.null(rds$run.par)){
		return(FALSE)
	}
	options(scipen = 999)

	tmp <- parameters[grep(pattern='^[0-9]', parameters$value),]
	numParameters <- data.frame(param=as.character(tmp$param), value=as.character(tmp$value), stringsAsFactors=FALSE)
	numParameters <- numParameters[-which(numParameters[,1] == "pipeline.version" | numParameters[,1] == "nThreads"),]
		
	tmp <- rds$run.par[grep(pattern='^[0-9]', rds$run.par$value),]
	numRdsParameters <- data.frame(param=as.character(tmp$param), valueRDS=as.character(tmp$value), stringsAsFactors=FALSE)
	numRdsParameters <- numRdsParameters[-which(numRdsParameters[,1] == "pipeline.version" | numRdsParameters[,1] == "nThreads"),]

	param <- left_join(numParameters, numRdsParameters, by="param")
	param[which(is.na(param$valueRDS) == TRUE), which(colnames(param) =="valueRDS")] <-  "NotFound"
	param$test <- ifelse(param$value == param$valueRDS, TRUE, FALSE)

	if(all(param$test) == TRUE){
		if (analysis == rds$vpInfo$analysis){
			return(TRUE)
		} else {
			return(FALSE)
		}		
	} else {
		diffParam <- param[which(param$test == FALSE),1]
		if ("min.amount.reads" %in% diffParam | "trim.length" %in% diffParam | "cutoff" %in% diffParam ){
			txt.tmp <- paste0(trim.F, expname, ".txt")
			info.file <- paste0(trim.F, expname, ".info.rds")
			if (file.exists(txt.tmp)){
				file.remove(txt.tmp)
			}
			if (file.exists(info.file)){
				file.remove(info.file)
			}
		}
		if ("reads.quality" %in% diffParam ){
			bamFile <- paste0(BAM.F, expname, ".bam")
			if (file.exists(bamFile)){
				file.remove(bamFile)
			}
		}
		return(FALSE)
	}	
}

make.DESeq2matrix <- function(subVPinfo, gR_positivePeaks, pathToRDS, region=FALSE, merge=FALSE){
	cond <- unique(subVPinfo$condition)
	gR_matrix <- GRanges(seqnames=seqnames(gR_positivePeaks), ranges=ranges(gR_positivePeaks), strand=strand(gR_positivePeaks))
	num_column = 1
	for (k in cond) {
		if (merge == TRUE){
			VP <- paste0("chr", unique(subVPinfo[which(subVPinfo$condition == k), colnames(VPinfo) == "vpchr"]), 
				"-", unique(subVPinfo[which(subVPinfo$condition == k), colnames(VPinfo) == "vppos"]))
			expname <- paste0("merge_condition_", k, "_viewpoint_", VP)
		} else {
			expname <- subVPinfo[which(subVPinfo$condition == k), colnames(VPinfo) == "expname"]
			num_replicates <- subVPinfo[which(subVPinfo$condition == k), colnames(VPinfo) == "replicate"]
		}

		for (y in 1:length(expname)) {
			rdsFile <- paste0(pathToRDS, expname[y], ".rds")					
			tmpRDS <- readRDS(rdsFile)
			reads <- tmpRDS$reads
			ranges(reads) <- reads$pos
			rm(tmpRDS)

			if (merge == TRUE){
				msg <- paste0("merge_", k)
			} else {
				rep <- num_replicates[y]
				msg <- paste0(k, "_", "rep", rep)
			}

			if(region == TRUE){
				score <- c()
				for (l in 1:length(gR_matrix)) {
					tmp <- gR_matrix[l]
					peaks <- as.data.frame(reads[unique(queryHits(findOverlaps(reads, tmp)))])
					score <- c(score, round(median(peaks$norm4C), digits=5))
				}
				gR_matrix$tmp <- score
				colnames(elementMetadata(gR_matrix))[num_column] <- paste0("norm4C_", msg)
			} else {
				positive <- reads[unique(queryHits(findOverlaps(reads, gR_matrix)))]				
				gR_matrix$tmp <- round(positive$norm4C, digits=5)
				colnames(elementMetadata(gR_matrix))[num_column] <- paste0("norm4C_", msg)
			}
			num_column = num_column + 1
		}
	}

	tmp <- as.data.frame(gR_matrix)
	matrix <- tmp[,-(1:5)]
	if (region == TRUE){
		rownames(matrix) <- paste0(tmp$seqnames, ":", tmp$start, "-", tmp$end)
	} else {
		rownames(matrix) <- paste0(tmp$seqnames, ":", tmp$start)
	}
	return(matrix)
}

do.deseq2 <- function(matrix, metadata, output, expname, ctrlCondition){
	for (i in 1:ncol(matrix)){
		matrix[,i] <- as.integer(matrix[,i])
	}
	metadata <- subVPinfo[order(metadata$condition), ]
	metadata$condition <- as.factor(metadata$condition)
	cond <- unique(metadata$condition)

	dds <- DESeqDataSetFromMatrix(countData = matrix, colData = metadata, design = ~ condition)
	dds$condition <- relevel(dds$condition, ref=ctrlCondition)
	
	# already normalized data
	sizeFactors(dds) <- 1

	dds <- estimateDispersions(dds, fitType="local", quiet=TRUE)
	pdf(paste0(output, expname, "_dispersion_estimation_local.pdf"))
	plotDispEsts(dds)
	invisible(dev.off())

	dds <- nbinomWaldTest(dds)

	# PCA analysis
	rld_dds <- rlog(dds, blind=FALSE)
	pdf(paste0(output, expname, "_PCA_analysis.pdf"))
	plotPCA(rld_dds, intgroup=c("condition")) + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
	invisible(dev.off())

	dataRes <- data.frame(interaction_site=rownames(matrix))
	dataRes$interaction_site <- as.character(dataRes$interaction_site)
	for (c in 1:length(cond)) {
		tmpCond <- as.character(cond[c])
		dataRes$tmp <- rowMeans(counts(dds, normalized=TRUE)[, dds$condition == tmpCond])
		colnames(dataRes)[(c+1)] <- paste0("baseMean_", tmpCond)
	}
	dataRes$tmp <- rowMeans(counts(dds, normalized=TRUE))
	colnames(dataRes)[(length(cond) + 2)] <- "baseMean"

	resNames <- resultsNames(dds)
	resNames <- resNames[-1]
	for (res in resNames) {
		tempRes <- as.data.frame(results(dds, alpha=0.1 , pAdjustMethod="BH", name=res, format="DataFrame"))
		tempRes <- tempRes[-1]
		CN <- gsub("condition_", "", res)
		colnames(tempRes) <- c(paste0(CN, ".log2FoldChange"), paste0(CN, ".lfcSE"), paste0(CN, ".stat"), paste0(CN, ".pvalue"), paste0(CN, ".padj"))
		tempRes$interaction_site <- as.character(rownames(tempRes))
		dataRes <- left_join(dataRes, tempRes, by="interaction_site")
	}
	for (i in 1:ncol(dataRes)) {
		if(class(dataRes[, i]) == "numeric"){
			dataRes[, i] <- round(dataRes[, i], digits=3)
		}
	}
	tableName <- paste0(output, expname, "_DESeq2_results.txt")
	write.table(dataRes, file=tableName, eol="\n", sep="\t", row.names=FALSE, quote=FALSE, col.names = TRUE)
}

Run.4Cpipeline <- function(VPinfo.file, FASTQ.F, OUTPUT.F, configuration){
	#4C-pipeline
	#1. demux.FASTQ
	#2. trim.FASTQ
	#3. mapped trimmer.files using bowtie2
	#4. extract unique reads
	#5. make RDS file
	#6. make plot
	#7. make WIG
	#8. make report


	cutoff = configuration$qualityCutoff
	trim.length = configuration$trimLength
	min.amount.reads=configuration$minAmountReads
	reads.quality=configuration$readsQuality
	map.unique = configuration$mapUnique
	nThreads = configuration$cores
	wSize = configuration$wSize
	nTop = configuration$nTop
	nonBlind = configuration$nonBlind
	make.bdg = configuration$bdg
	make.wig = configuration$wig
	make.fixedStep = configuration$fixedStepWig
	fixedStepWigBin = configuration$fixedStepWigBin
	make.cisplot = configuration$cisplot
	make.gwplot = configuration$genomePlot
	tsv = configuration$tsv
	bins = configuration$bins
	mmMax = configuration$mmMax
	normFactor = configuration$normFactor

	peakC = configuration$peakC
	replicates = configuration$replicates

	chr_random = configuration$chr_random
	chrUn = configuration$chrUn
	chrM = configuration$chrM

	vpRegion = configuration$vpRegion
	alphaFDR = configuration$alphaFDR
	qWd = configuration$qWd
	qWr = configuration$qWr
	minDist = configuration$minDist
	minGapwidth = configuration$min.gapwidth
	DESeq2 = configuration$DESeq2
	ctrl = configuration$ctrl


	# create folders
	# define folder names
	LOG.F <- gsub(x=paste0(OUTPUT.F, "/LOG/"), pattern='//', replacement='/')
	FASTQ.demux.F <- gsub(x=paste0(OUTPUT.F, "/FASTQ/"), pattern='//', replacement='/')
	TRIM.F <-gsub(x=paste0(OUTPUT.F, "/TRIMMED/"), pattern='//', replacement='/')
	BAM.F <-gsub(x=paste0(OUTPUT.F, "/BAM/"), pattern='//', replacement='/')
	RDS.F <-gsub(x=paste0(OUTPUT.F, "/RDS/"), pattern='//', replacement='/')
	RDS.MERGE.F <-gsub(x=paste0(OUTPUT.F, "/RDS/MergeRDS/"), pattern='//', replacement='/')
	RDS.BIN.F <- gsub(x=paste0(OUTPUT.F, "/RDS-BIN/"), pattern='//', replacement='/')
	RDS.BIN.MERGE.F <- gsub(x=paste0(OUTPUT.F, "/RDS-BIN/MergeRDS-BIN/"), pattern='//', replacement='/')
	PLOT.F <- gsub(x=paste0(OUTPUT.F, "/PLOT/"), pattern='//', replacement='/')
	PLOT.MERGE.F <- gsub(x=paste0(OUTPUT.F, "/PLOT/MergePLOT/"), pattern='//', replacement='/')
	BDG.F <- gsub(x=paste0(OUTPUT.F, "/BDG/"), pattern='//', replacement='/')
	BDG.MERGE.F <- gsub(x=paste0(OUTPUT.F, "/BDG/MergeBDG/"), pattern='//', replacement='/')
	WIG.F <- gsub(x=paste0(OUTPUT.F, "/WIG/"), pattern='//', replacement='/')
	WIG.MERGE.F <- gsub(x=paste0(OUTPUT.F, "/WIG/MergeWIG/"), pattern='//', replacement='/')
	GENOMEPLOT.F <- gsub(x=paste0(OUTPUT.F, "/GENOMEPLOT/"), pattern='//', replacement='/')
	GENOMEPLOT.MERGE.F <- gsub(x=paste0(OUTPUT.F, "/GENOMEPLOT/MergeGENOMEPLOT/"), pattern='//', replacement='/')
	TSV.F <- gsub(x=paste0(OUTPUT.F, "/TSV/"), pattern='//', replacement='/')
	TSV.MERGE.F <- gsub(x=paste0(OUTPUT.F, "/TSV/MergeTSV/"), pattern='//', replacement='/')

	PEAKC.F <- gsub(x=paste0(OUTPUT.F, "/PEAKC/"), pattern='//', replacement='/')
	PEAKC.PLOTS.F <- gsub(x=paste0(OUTPUT.F, "/PEAKC/Plots/"), pattern='//', replacement='/')
	PEAKC.BED_BDG.F <- gsub(x=paste0(OUTPUT.F, "/PEAKC/BED_and_BDG/"), pattern='//', replacement='/')
	PEAKC.RDS.F <- gsub(x=paste0(OUTPUT.F, "/PEAKC/RDS/"), pattern='//', replacement='/')
	PEAKC.DESeq2.F <- gsub(x=paste0(OUTPUT.F, "/PEAKC/DESeq2/"), pattern='//', replacement='/')

	logDirs <- list()

	##To do: replace ifelse for if statement or give FALSE statement...
	logDirs$outFolder <- ifelse(!dir.exists(OUTPUT.F), dir.create(OUTPUT.F), FALSE)
	logDirs$logFolder <- ifelse(!dir.exists(LOG.F), dir.create(LOG.F), FALSE)
	logDirs$fastqFolder <- ifelse(!dir.exists(FASTQ.demux.F), dir.create(FASTQ.demux.F), FALSE)
	logDirs$trimFolder <- ifelse(!dir.exists(TRIM.F), dir.create(TRIM.F), FALSE)
	logDirs$bamFolder <- ifelse(!dir.exists(BAM.F), dir.create(BAM.F), FALSE)
	logDirs$rdsFolder <- ifelse(!dir.exists(RDS.F), dir.create(RDS.F), FALSE)
	logDirs$rdsbinFolder <- ifelse(!dir.exists(RDS.BIN.F), dir.create(RDS.BIN.F), FALSE)

	# Extract experiments from VPinfo, remove weird characters and save in output folder. 
	message(paste0("\n------ Reading the VP info file: ", VPinfo.file))
	VPinfo <- Read.VPinfo(VPinfo.file=VPinfo.file, replicates=replicates)
	if (is.null(VPinfo)){
		stop("viewpoint info file (vpFile) not correct.")
	}
	write.table(VPinfo, paste0(OUTPUT.F, "/VPinfo.txt"), sep="\t", row.names=FALSE, quote=F)

	if (replicates == TRUE){
		logDirs$rdsMergeFolder <- ifelse(!dir.exists(RDS.MERGE.F), dir.create(RDS.MERGE.F), FALSE)
		logDirs$rdsbinMergeFolder <- ifelse(!dir.exists(RDS.BIN.MERGE.F), dir.create(RDS.BIN.MERGE.F), FALSE)
	}

	if (make.cisplot == TRUE){
		logDirs$plotFolder <- ifelse(!dir.exists(PLOT.F), dir.create(PLOT.F), FALSE)
		if(replicates == TRUE){
			logDirs$plotMergeFolder <- ifelse(!dir.exists(PLOT.MERGE.F), dir.create(PLOT.MERGE.F), FALSE)
		}
	}
	if (make.bdg == TRUE){
		logDirs$bdgFolder <- ifelse(!dir.exists(BDG.F), dir.create(BDG.F), FALSE)
		if(replicates == TRUE){
			logDirs$bdgMergeFolder <- ifelse(!dir.exists(BDG.MERGE.F), dir.create(BDG.MERGE.F), FALSE)
		}
	}
	if (make.wig == TRUE){
		logDirs$wigFolder <- ifelse(!dir.exists(WIG.F), dir.create(WIG.F), FALSE)
		if(replicates == TRUE){
			logDirs$wigMergeFolder <- ifelse(!dir.exists(WIG.MERGE.F), dir.create(WIG.MERGE.F), FALSE)
		}
	}
	if(any(VPinfo$analysis == 'all') & make.gwplot){
		logDirs$genomeplotFolder <- ifelse(!dir.exists(GENOMEPLOT.F), dir.create(GENOMEPLOT.F), FALSE)
		if(replicates == TRUE){
			logDirs$genomeplotMergeFolder <- ifelse(!dir.exists(GENOMEPLOT.MERGE.F), dir.create(GENOMEPLOT.MERGE.F), FALSE)
		}
	}
	if (tsv == TRUE){
		logDirs$tsvFolder <- ifelse(!dir.exists(TSV.F), dir.create(TSV.F), FALSE)
		if(replicates == TRUE){
			logDirs$tsvMergeFolder <- ifelse(!dir.exists(TSV.MERGE.F), dir.create(TSV.MERGE.F), FALSE)
		}
	}
	if (peakC == TRUE){
		if(!suppressMessages(require("peakC", character.only=TRUE))) stop("Package not found: peakC")
		if(!suppressMessages(require("isotone", character.only=TRUE))) stop("Package not found: isotone")
		logDirs$peakcFolder <- ifelse(!dir.exists(PEAKC.F), dir.create(PEAKC.F), FALSE)
		logDirs$peakcFolderPlots <- ifelse(!dir.exists(PEAKC.PLOTS.F), dir.create(PEAKC.PLOTS.F), FALSE)
		logDirs$peakcFolderBED_BDG <- ifelse(!dir.exists(PEAKC.BED_BDG.F), dir.create(PEAKC.BED_BDG.F), FALSE)
		logDirs$peakcFolderRDS <- ifelse(!dir.exists(PEAKC.RDS.F), dir.create(PEAKC.RDS.F), FALSE)
		if (DESeq2 == TRUE){
			if(!suppressMessages(require("DESeq2", character.only=TRUE))) stop("Package not found: DESeq2")
			if(!suppressMessages(require("ggplot2", character.only=TRUE))) stop("Package not found: ggplot2")
			if(!suppressMessages(require("ggrepel", character.only=TRUE))) stop("Package not found: ggrepel")
			logDirs$peakcFolderDESeq2 <- ifelse(!dir.exists(PEAKC.DESeq2.F), dir.create(PEAKC.DESeq2.F), FALSE)
		}
	}

	# Make Log files
	LOGNAME <- format(Sys.time(), '%Y_%m_%d__%H_%M_%S')
	log.path <- paste0(LOG.F, "Log_pipe4C_", LOGNAME, ".txt")
	demux.log.path <- paste0(LOG.F, "Log_demux_", LOGNAME, ".txt")
	trim.log.path <- paste0(LOG.F, "Log_4Ctrim_", LOGNAME, ".txt")
	bowtie.log.path <- paste0(LOG.F, "Log_Bowtie2", LOGNAME, ".txt")

	# Save Run parameters to logfile
	run.par <- data.frame(
		param=c("pipeline.version", "baseFolder", "VPinfo.file", "FASTQ.F", "OUTPUT.F", "normFactor", "cutoff", "trim.length",
			"minAmountReads", "reads.quality", "map.unique", "wSize", "nTop", "make.bdg", "make.wig", "make.fixedStep", "fixedStepWigBin",
			"make.cisplot", "make.gwplot", "nThreads", "nonBlind", "tsv","bins", "mmMax", "chr_random", "chrUn", "chrM", "peakC",
			"replicates", "vpRegion", "alphaFDR", "qWd", "qWr", "minDist", "min.gapwidth", "DESeq2"),
		value=c(configuration$pipeline.version, configuration$baseFolder, VPinfo.file, FASTQ.F, OUTPUT.F, normFactor, cutoff, trim.length, 
			min.amount.reads, reads.quality, map.unique, wSize, nTop, make.bdg, make.wig, make.fixedStep, fixedStepWigBin, make.cisplot, make.gwplot, 
			nThreads, nonBlind, tsv,bins, mmMax, chr_random, chrUn, chrM, peakC, replicates, vpRegion, alphaFDR, qWd, qWr, minDist, minGapwidth, DESeq2))

	write.table(run.par, log.path, quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

	#4C-seq analysis
	message("\n------ 4Cseq analysis\n")

	exp.name <- as.character(VPinfo$expname)
	primer <- as.character(VPinfo$primer)
	genome <- as.character(VPinfo$genome)
	firstenzyme <-as.character(VPinfo$firstenzyme)
	secondenzyme <- as.character(VPinfo$secondenzyme)
	vpChr <- as.character(VPinfo$vpchr)
	vppos <- as.numeric(VPinfo$vppos)
	analysis <- as.character(VPinfo$analysis)

	#Check for existing analysis
	message(" +++ Checking existing analysis +++")
	needToLoad <- c()
	for (i in 1:length(exp.name)){
		rdsFile<-paste0(RDS.F, exp.name[i], ".rds")
		
		if (file.exists(rdsFile)){
			tmp <- testRDS(rdsFile=rdsFile, parameters=run.par, expname=exp.name[i], analysis=analysis[i])
			needToLoad <- c(needToLoad, tmp)
		} else {
			needToLoad <- c(needToLoad, FALSE)
		}
	}

	#1. demux.FASTQ
	#Demultiplex all fastq files and write in FASTQ folder in OUTPUT.F
	message(" +++ Demultiplexing Fastq files based on VPinfo file +++")
	demux.FASTQ(VPinfo=VPinfo, FASTQ.F=FASTQ.F, FASTQ.demux.F=FASTQ.demux.F, demux.log.path=demux.log.path, mmMax = mmMax)

	if (peakC == TRUE & DESeq2 == TRUE){
		list_peakCPeaksPos_Per_VP <- list()
		list_peakCRegionsPos_Per_VP <- list()
		if (replicates == TRUE){
			list_peakCMergePeaksPos_Per_VP <- list()
			list_peakCMergeRegionsPos_Per_VP <- list()
		}
	}

	for (i in 1:length(exp.name)){
		message(paste0(" +++ experiment: ", exp.name[i], " +++"))

		rdsFile<-paste0(RDS.F, exp.name[i], ".rds")
		CHR <- paste0("chr", vpChr[i])

		if (file.exists(rdsFile) & needToLoad[i]){
			error.msg <- paste0("         ### WARNING: rds file ", exp.name[i], " already exists with identical parameters, continuing with existing file.")
			write(error.msg, log.path, append=TRUE)
			message(error.msg)
			rds <- readRDS(rdsFile)
			reportAnalysis <- rds$report
			vpInfo <- rds$vpInfo

			if (vpInfo$analysis == "cis"){
				reads.cis <- rds$reads				
			} else if (vpInfo$analysis == "all"){
				reads.all <- rds$reads
				reads.cis <- norm4C(readsGR=reads.all[as.vector(seqnames(reads.all)) == CHR,], nReads=normFactor, wSize=wSize, nTop=nTop)
			}
			rm(rds)
		} else {
			if (exp.name[i] %in% exp.name[duplicated(exp.name)]){
				error.msg <- paste0("      ### ERROR: Experiment name not unique for ", exp.name[i])
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}
			primer.sequence <- primer[i]

			if (!(firstenzyme[i] %in% rownames(configuration$enzymes))){
				error.msg <- paste0("      ### ERROR: firstcutter not found for ", exp.name[i])
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}
			firstcutter <- configuration$enzyme[firstenzyme[i],]

			if (!(secondenzyme[i] %in% rownames(configuration$enzymes))){
				error.msg <- paste0("      ### ERROR: secondcutter not found for ", exp.name[i])
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}
			secondcutter <- configuration$enzyme[secondenzyme[i],]

			if (!(genome[i] %in% rownames(configuration$genomes))){
				error.msg <- paste0("      ### ERROR: genome not found for ", exp.name[i])
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}

			if (!(analysis[i] %in% c("cis", "all"))){
				error.msg <- paste0("      ### ERROR: analysis of cis/all not found for ", exp.name[i])
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}

			file.fastq <- paste0(FASTQ.demux.F, exp.name[i], ".fastq.gz")

			if (!file.exists(file.fastq)){
				error.msg <- paste0("      ### ERROR: ", exp.name[i], " FASTQ file does not exist")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}


			#2. trim.FASTQ
			message("      >>> Trim of the fastq <<<")
			trimFASTQ <- trim.FASTQ(exp.name=exp.name[i], firstcutter=firstcutter, secondcutter=secondcutter, file.fastq=file.fastq, trim.F=TRIM.F, 
				cutoff=cutoff, trim.length=trim.length, log.path=trim.log.path, min.amount.reads=min.amount.reads)

			txt.tmp <- paste0(TRIM.F, exp.name[i], ".txt")
			if (!file.exists(txt.tmp)){
				error.msg <- paste0("         ### ERROR: ", exp.name[i], " trimmed file not found")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}

			captureLen <- trimFASTQ$captureLen

			if (captureLen<=nchar(firstcutter)){
				error.msg <- paste0("         ### ERROR: ", exp.name[i], " capture length <= length firstcutter motif.")
				message("         ### Skipping experiment.")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}

			nReads <- trimFASTQ$nReads
			motifPosperc <-trimFASTQ$motifPosperc
			readlenperc <- trimFASTQ$readlenperc


			# 3. make BAM files
			message("      >>> Alignment of reads to reference genome <<<")
			bamFile <- paste0(BAM.F, exp.name[i], ".bam")

			if (!file.exists(bamFile)){
				makeBAM(exp.name=exp.name[i], BAM.F=BAM.F, NCORES=nThreads, Bowtie2IndexFile=configuration$bt2Genomes[genome[i],], txt.tmp=txt.tmp, 
					log.path=log.path, bowtie.log.path=bowtie.log.path, bamFile=bamFile, map.unique=map.unique, readsQual=reads.quality)
			} else {
				error.msg <- paste("         ### WARNING: ", exp.name[i], " BAM file already exist, continuing with exisiting file.")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
			}

			mappedReads <- readGAlignments(file=bamFile, index=bamFile)

			# To do: lock frag map
			#exec 3>hg19_Dpn2_Csp6I # open a file handle; this part will always succeed
			#flock -x 3      # lock the file handle; this part will block
			#To release the lock:
			#exec 3>&-       # close the file handle

			message(paste0("      >>> Create frag map for genome: ", genome[i], " with RE1:", firstcutter, " and RE2:", secondcutter, " and capture length:", captureLen, " <<<"))
			frags <- getFragMap(vpChr_FragMap=NULL,	firstcutter_FragMap=firstcutter, secondcutter_FragMap=secondcutter, genome_FragMap=genome[i],	captureLen_FragMap=captureLen,
				nThreads_FragMap=nThreads, baseFolder_FragMap=configuration$baseFolder,	Bowtie2IndexFile=configuration$bt2Genomes[genome[i],],	config_genomes=configuration$genomes,
				chr_random=chr_random, chrUn=chrUn, chrM = chrM)

			message("      >>> Align reads to fragments <<<")
			readsAln <- alignToFragends(gAlign=mappedReads, fragments=frags, firstcut=firstcutter)
			rm(frags)

			message("      >>> Compute statistics and create report <<<")
			reportAnalysis <- createReport(allReads=readsAln, mapReads=mappedReads, demuxReads=nReads, chromosome=CHR, vpPos=vppos[i], normFactor=normFactor, 
				wSize=wSize, nTop=nTop, motifPosperc=motifPosperc, readlenperc=readlenperc)

			vpInfo <- data.frame(name=exp.name[i], Fastq=file.fastq, primer=primer.sequence, genome=genome[i], firstcutter=firstcutter, secondcutter=secondcutter, 
				chr=CHR, pos=vppos[i], captureLen=captureLen, analysis=analysis[i], date=Sys.time(), Pipeline=configuration$pipeline.version)


			#Extract Cis reads.
			message("      >>> Compute normalized 4C score per fragment <<<")

			#Data normalization
			if (analysis[i] == "cis" | make.cisplot == TRUE){
				reads.cis <- norm4C(readsGR=readsAln[as.vector(seqnames(readsAln)) == CHR], nReads=normFactor, wSize=wSize, nTop=nTop)
				if (nonBlind){
					reads.cis<-reads.cis[reads.cis$type=="non_blind"]
				}
			}		

			if (analysis[i] == "all" | bins == TRUE | make.gwplot == TRUE){
				reads.all <- norm4C(readsGR=readsAln, nReads=normFactor, wSize=wSize, nTop=nTop)
				if (nonBlind){
					reads.all<-reads.all[reads.all$type=="non_blind"]
				}
			}
			rm(readsAln)

			#5. make RDS file
			if (analysis[i] == "cis"){
				saveRDS(list(reads=reads.cis, report=reportAnalysis, vpInfo=vpInfo, run.par=run.par), file=rdsFile)
			} else if (analysis[i] == "all"){
				saveRDS(list(reads=reads.all, report=reportAnalysis, vpInfo=vpInfo, run.par=run.par), file=rdsFile)
			}
		}
		
		if (bins == TRUE | make.gwplot == TRUE){
			if (analysis[i] == "all"){
				message("      >>> Creating bins <<<")	
				if (configuration$binSize < 1e3){
					file.out <- paste0(RDS.BIN.F, exp.name[i], "_bin_res_", configuration$binSize, "bp.rds")
				} else {
					file.out <- paste0(RDS.BIN.F, exp.name[i], "_bin_res_", (configuration$binSize/1e3), "kb.rds")
				}				
				bin.GR <- doBins(file=file.out,  expname=exp.name[i], bin=configuration$binSize , reads=reads.all, assemblyName=genome[i], 
					config_genomes=configuration$genomes, report=reportAnalysis, vpInfo=vpInfo, log.path=log.path)
			}
		}
		
		#6. make plot
		if (make.cisplot == TRUE){
			message("      >>> Creating local 4C Plot <<<")
			createPlot(plotTitle=exp.name[i], vpPos=vppos[i], chromosome=CHR, fragGR=reads.cis, plotLegend=reportAnalysis, plotView=configuration$plotView,
				maxY=configuration$maxY, minY=0, xaxisUnit=configuration$xaxisUnit, plotRegion='cis', foldOut=PLOT.F, plotType=configuration$plotType)
		}

		if (analysis[i]=="all" & make.gwplot == TRUE){
			message("      >>> Creating genome-wide 4C coverage Plot <<<")
			createPlot(plotTitle=exp.name[i], vpPos=vppos[i], chromosome=CHR, fragGR=bin.GR, plotLegend=reportAnalysis, plotView=configuration$plotView,
				maxY=configuration$maxY, minY=0, xaxisUnit=configuration$xaxisUnit, plotRegion='all', foldOut=GENOMEPLOT.F, plotType=configuration$plotType)
		}

		#Export data
		if (make.bdg == TRUE){
			message("      >>> Creating bedGraph file <<<")

			if (nonBlind){
				bdgFile <- paste0(BDG.F, exp.name[i], "_nonblind_WIN", wSize,".bedGraph")
			} else {
				bdgFile <- paste0(BDG.F, exp.name[i], "_WIN", wSize ,".bedGraph")
			}
			gzBDGName <- paste0(bdgFile, ".gz")

			if (file.exists(bdgFile) | file.exists(gzBDGName)){
				error.msg <- paste0("         ### WARNING: ", exp.name[i], " bedGraph file already exists.")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
			} else {
				if (analysis[i] == "cis"){
					reads.bdg <- reads.cis[order(reads.cis$pos)]
				} else {
					reads.bdg <- reads.all[order(seqnames(reads.all), reads.all$pos)]
				}
				message("         ### Generation of bedGraph file.")
				exportBDG(gR=reads.bdg, expName=exp.name[i], filename=bdgFile, vpPos=vppos[i], vpChr=CHR, 
					plotView=configuration$plotView, interval=FALSE)
			}
		}

		if (make.wig == TRUE){
			message("      >>> Creating WIG file <<<")

			if (nonBlind){
				wigFile <- paste0(WIG.F, exp.name[i], "_nonblind_WIN", wSize, ".wig")
			} else {
				wigFile <- paste0(WIG.F, exp.name[i], "_WIN", wSize , ".wig")
			}
			gzWigName <- paste0(wigFile, ".gz")
			bwName <- gsub(".wig", ".bw", wigFile)

			if (file.exists(wigFile) | file.exists(gzWigName) & file.exists(bwName)){
				error.msg <- paste0("         ### WARNING: ", exp.name[i], " WIG and BIGWIG files already exist.")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
			} else if (file.exists(wigFile) | file.exists(gzWigName) & !file.exists(bwName)){
				error.msg <- paste0("         ### WARNING: ", exp.name[i], " WIG file already exists.")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
				if (file.exists(wigFile)){
					do.wigToBigWig(wig_file=wigFile, config_genomes=configuration$genomes, assemblyName=genome[i])
				} else {
					do.wigToBigWig(wig_file=gzWigName, config_genomes=configuration$genomes, assemblyName=genome[i])
				}				
			} else {
				if (analysis[i] == "cis"){
					reads.wig <- reads.cis[order(reads.cis$pos)]
				} else {
					reads.wig <- reads.all[order(seqnames(reads.all),reads.all$pos)]
				}
				message("         ### Generation of variableStep WIG file.")
				exportWig(gR=reads.wig, expName=exp.name[i], filename=wigFile, vpPos=vppos[i], vpChr=CHR, plotView=configuration$plotView, 
					config_genomes=configuration$genomes, assemblyName=genome[i], fixed=FALSE)
			}

			if (make.fixedStep){
				if (nonBlind){
					wigBinFile <- paste0(WIG.F, exp.name[i], "_nonblind_BIN_", fixedStepWigBin, ".wig")
				} else {
					wigBinFile <- paste0(WIG.F, exp.name[i], "_BIN", fixedStepWigBin, ".wig")
				}

				gzWigBinName <- paste0(wigBinFile, ".gz")
				bwBinName <- gsub(".wig", ".bw", wigBinFile)

				if (file.exists(wigBinFile) | file.exists(gzWigBinName) & file.exists(bwBinName)){
					error.msg <- paste0("         ### WARNING: ", exp.name[i], " fixedStep WIG and BIGWIG file already exist.")
					message(error.msg)
					write(error.msg, log.path, append=TRUE)
				} else if (file.exists(wigBinFile) | file.exists(gzWigBinName) & !file.exists(bwBinName)){
					error.msg <- paste0("         ### WARNING: ", exp.name[i], " fixedStep WIG file already exists.")
					message(error.msg)
					write(error.msg, log.path, append=TRUE)
					if (file.exists(wigBinFile)){
						do.wigToBigWig(wig_file=wigBinFile, config_genomes=configuration$genomes, assemblyName=genome[i])
					} else {
						do.wigToBigWig(wig_file=gzWigBinName, config_genomes=configuration$genomes, assemblyName=genome[i])
					}
				} else {
					if (fixedStepWigBin < 1e3){
						file.out2 <- paste0(RDS.BIN.F, exp.name[i], "_bin_res_", fixedStepWigBin, "bp.rds")
					} else {
						file.out2 <- paste0(RDS.BIN.F, exp.name[i], "_bin_res_", (fixedStepWigBin/1e3), "kb.rds")
					}

					bin.GR.fixedWig <- doBins(file=file.out2,  expname=exp.name[i], bin=fixedStepWigBin, reads=reads.all, assemblyName=genome[i], 
					config_genomes=configuration$genomes, report=reportAnalysis, vpInfo=vpInfo, log.path=log.path)

					message("         ### Generation of fixedStep WIG file.")
					exportWig(gR=bin.GR.fixedWig$bins, expName=exp.name[i], filename=wigBinFile, vpPos=vppos[i], vpChr=CHR, plotView=configuration$plotView, 
						config_genomes=configuration$genomes, assemblyName=genome[i], fixed=TRUE, bin=fixedStepWigBin)
				}
			}
		}

		if (tsv == TRUE){
			message("      >>> Creating TSV file <<<")

			if (nonBlind){
				tsvFile <- paste0(TSV.F, exp.name[i], "_nonblind_WIN", wSize, ".tsv")
			} else {
				tsvFile <- paste0(TSV.F, exp.name[i], "_WIN", wSize , ".tsv")
			}
			gzTSVName <- paste0(tsvFile, ".gz")

			if (file.exists(tsvFile) | file.exists(gzTSVName)){
				error.msg <- paste0("         ### WARNING: ", exp.name[i], " TSV file already exists.")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
			} else {
				if (analysis[i] == "cis"){
					reads.tsv <- reads.cis[order(reads.cis$pos)]
				} else {
					reads.tsv <- reads.all[order(seqnames(reads.all),reads.all$pos)]
				}
				exportTSV(gR=reads.tsv, filename=tsvFile, merge=FALSE)
			}
		}

		# PeakC analysis
		if (peakC == TRUE){
			message("      >>> PeakC analysis <<<")
			bedFile = paste0(PEAKC.BED_BDG.F, exp.name[i], "_peakC_peaks.bed")
			bedFileRegion = gsub("peakC_peaks", "peakC_regions", bedFile)
			BedGraphPeakCFile = paste0(PEAKC.BED_BDG.F, exp.name[i],  "_peakC_peaks.bedGraph")
			BedGraphFileRegion = gsub("peakC_peaks", "peakC_regions", BedGraphPeakCFile)
			PeakCRdsFile = paste0(PEAKC.RDS.F, exp.name[i], "_peakC_results.rds")

			if (analysis[i] == "cis"){
				reads <- reads.cis
			} else {
				reads <- reads.all
			}

			if (file.exists(PeakCRdsFile)){
				error.msg <- paste0("         ### WARNING: ", exp.name[i], " peakC analysis already exists, continuing with the existing RDS file.")
				message(error.msg)
				write(error.msg, log.path, append=TRUE)
				resPeakC <- readRDS(PeakCRdsFile)
				VPpos <- resPeakC$vpPos
				VPChr <- resPeakC$vpChr						
			} else {
				VPpos <- vppos[i]
				VPChr <- as.vector(CHR)

				zoom <- GRanges(seqnames=VPChr, resize(IRanges(VPpos, VPpos), width=vpRegion, fix="center"))
				vpGR <- reads[unique(queryHits(findOverlaps(reads,zoom)))]
				peakCDat <- data.frame(pos=vpGR$pos,reads=vpGR$normReads)
				
				resPeakC <- suppressWarnings(single.analysis(data=peakCDat,vp.pos=VPpos,wSize=wSize,qWd=qWd,qWr=qWr,minDist=minDist))

				resPeakC$vpPos <- VPpos
				resPeakC$vpChr <- VPChr

				if(length(resPeakC$peak)>0){
					resPeakC$exportPeakGR <- getPeakCPeaks(resPeakC=resPeakC, min.gapwidth=0)
				} else {
					resPeakC$exportPeakGR <- NULL
				}

				#Save peakC results
				saveRDS(resPeakC, file=PeakCRdsFile)
			}

			#Plot peakC results
			plotFile = paste0(PEAKC.PLOTS.F,"peakC_", exp.name[i],".pdf")
			pdf(file=plotFile)
			plot_C(data=resPeakC, num.exp = 1, y.min = 0, y.max = 750)
			dev.off()

			#Export bed file
			bedPeaks <- exportPeakCPeaks(resPeakC=resPeakC, bedFile=bedFile,
				name=paste0(exp.name[i],"peakC_peaks"), desc=NULL, includeVP=TRUE, min.gapwidth=0)
			bedregion <- exportPeakCPeaks(resPeakC=resPeakC, bedFile=bedFileRegion,
				name=paste0(exp.name[i],"peakC_region"), desc=NULL, includeVP=TRUE, min.gapwidth=minGapwidth)

			#export bedGraph
			reads2 <- reads
			ranges(reads2) <- reads2$pos
			peakCReads <-  GRanges(seqnames=bedPeaks$chr, ranges=bedPeaks$start)
			readsPeakC.bdg <- reads2[unique(queryHits(findOverlaps(reads2, peakCReads)))]
			
			peakCRegions <-  GRanges(seqnames=bedregion$chr, ranges=paste0(bedregion$start, "-", bedregion$end))
			score <- c()
			for (l in 1:length(peakCRegions)) {
				tmp <- peakCRegions[l]
				peaks <- as.data.frame(reads2[unique(queryHits(findOverlaps(reads2, tmp)))])
				score <- c(score, median(peaks$norm4C))
			}
			regionsPeakC.bdg <- peakCRegions
			regionsPeakC.bdg$norm4C <- score
			rm(tmp, peaks, score)

			exportBDG(gR=readsPeakC.bdg, expName=paste0(exp.name[i], "_peakC_peaks"), filename=BedGraphPeakCFile, 
				vpPos=VPpos, vpChr=VPChr, plotView=configuration$plotView, interval=FALSE)
			exportBDG(gR=regionsPeakC.bdg, expName=paste0(exp.name[i], "_peakC_regions"), filename=BedGraphFileRegion,
				vpPos=VPpos, vpChr=VPChr, plotView=configuration$plotView, interval=TRUE)

			#Preparing data to DESeq2 (recovery of all positive position per viewpoint)
			if (DESeq2 == TRUE){
				VP <- paste0(CHR, ":", vppos[i])
				if (is.null(list_peakCPeaksPos_Per_VP[[VP]]) == TRUE){
					dataReads <- as.data.frame(readsPeakC.bdg)
					list_peakCPeaksPos_Per_VP[[VP]] <- paste0(dataReads$seqnames, ":", dataReads$start)
				} else {
					dataReads <- as.data.frame(readsPeakC.bdg)
					dataReadsVP <- paste0(dataReads$seqnames, ":", dataReads$start)
					list_peakCPeaksPos_Per_VP[[VP]] <- unique(c(list_peakCPeaksPos_Per_VP[[VP]], dataReadsVP))
				}
				if (is.null(list_peakCRegionsPos_Per_VP[[VP]]) == TRUE){
					dataRegionsReads <- as.data.frame(regionsPeakC.bdg)
					list_peakCRegionsPos_Per_VP[[VP]] <- paste0(dataRegionsReads$seqnames, ":", dataRegionsReads$start, "-", dataRegionsReads$end)
				} else {
					dataRegionsReads <- as.data.frame(regionsPeakC.bdg)
					dataRegionsReadsVP <- paste0(dataRegionsReads$seqnames, ":", dataRegionsReads$start, "-", dataRegionsReads$end)
					list_peakCRegionsPos_Per_VP[[VP]] <- unique(c(list_peakCRegionsPos_Per_VP[[VP]], dataRegionsReadsVP))
				}
			}
			rm(peakCReads, readsPeakC.bdg, peakCRegions, regionsPeakC.bdg)		
		}
		message('\n')
	}

	if(replicates == TRUE){
		message("\n------ Replicate analysis")
		message("      >>> Recovery replicates informations <<<")
		
		to_process <- readConditions(vpInfo=VPinfo, log.path=log.path)
		cond <- as.character(unique(VPinfo$condition))

		for (k in cond){
			message(paste0("      >>> Processing the condition: ", k, " <<<"))
			for (z in 1:length(to_process[[k]])){
				disabled = 0
				expname <- as.vector(to_process[[k]][[z]])
				ViewP <- names(to_process[[k]][z])
				mergeExpname <- paste0(k, "_viewpoint_chr", ViewP)
				num_replicates <- VPinfo[which(VPinfo$expname == expname), colnames(VPinfo) == "replicate"]
				msg <-":"
				for (t in expname){
					msg <- paste0(msg, " | ", t)
				}
				message(paste0("         ### Experiment list to merge", msg ))
				mergeRdsFile <- paste0(RDS.MERGE.F, "merge_condition_", mergeExpname, ".rds")
				if (file.exists(mergeRdsFile)){
					error.msg <- paste0("            *** WARNING: rds file ", mergeExpname, " already exists, continuing with existing file.")
					write(error.msg, log.path, append=TRUE)
					message(error.msg)
					rds <- readRDS(mergeRdsFile)					
					merge <- rds$reads
					reportAnalysis <- rds$report					
					vpInfo <- rds$vpInfo
					expnameRDS <- vpInfo$name
					rm(rds)
					list_rdsFiles <- c()
					for (y in 1:length(expnameRDS)) {
						rdsFile <- paste0(RDS.F, expnameRDS[y], ".rds")
						list_rdsFiles <- c(list_rdsFiles, rdsFile)
					}
				} else {
					list_rdsFiles <- c()
					n=0
					for (y in 1:length(expname)) {
						rdsFile <- paste0(RDS.F, expname[y], ".rds")
						rep <- num_replicates[y]
						tmpRDS <- readRDS(rdsFile)
						replicateName <- paste0("rep", rep)
						if (y == 1) {
							merge <- tmpRDS$reads[,c("pos", "type", "fe_strand", "fe_id", "unique", "normReads", "norm4C")]
							colnames(elementMetadata(merge))[6] <- paste0("normReads_", replicateName)
							colnames(elementMetadata(merge))[7] <- paste0("norm4C_", replicateName)
							tmp <- as.data.frame(merge[,6:7])
							normRead <- as.data.frame(tmp[, 6])
							readsNorm4C <- as.data.frame(tmp[, 7])
							vpInfo <- tmpRDS$vpInfo
							reportAnalysis <- tmpRDS$report
						} else {
							numCol <- 5 + y + n
							merge$tmp <- tmpRDS$reads$normReads
							merge$tmp2 <- tmpRDS$reads$norm4C
							colnames(elementMetadata(merge))[numCol] <- paste0("normReads_", replicateName)
							colnames(elementMetadata(merge))[numCol+1] <- paste0("norm4C_", replicateName)
							tmp <- as.data.frame(merge[, numCol:(numCol+1)])
							normRead <- cbind(normRead, as.data.frame(tmp[, 6]))
							readsNorm4C <- cbind(readsNorm4C, as.data.frame(tmp[, 7]))
							vpInfo <- rbind(vpInfo, tmpRDS$vpInfo)
							reportAnalysis <- rbind(reportAnalysis, tmpRDS$report)
						}
						n=n+1
						list_rdsFiles <- c(list_rdsFiles, rdsFile)
						rm(tmpRDS)
					}

					if (length(unique(vpInfo$analysis)) != 1){
						error.msg <- paste0("            *** ERROR: ", msg, " files present different type of analysis (cis and all). Please uniform the analysis type between these replicates.")
						message("            *** Skipping replicates")
						write(error.msg, log.path, append=TRUE)
						message(error.msg)
						next
					}

					message("            ***  Merging normalized reads per fragment")
					#### Median as peakC package do for the plot_c function...
					merge$normReads <- apply(normRead, 1, median)
					merge$norm4C <- apply(readsNorm4C, 1, median)

					saveRDS(list(reads=merge, report=reportAnalysis, vpInfo=vpInfo, run.par=run.par), file=mergeRdsFile)
				}
				
				if (unique(vpInfo$analysis) == "cis") {
					disabled <- 1
				}

				if (bins == TRUE | make.gwplot == TRUE & disabled == 0){
					message("            *** Creating bins")	
					if (configuration$binSize < 1e3){
						file.out <- paste0(RDS.BIN.F, "merge_condition_", mergeExpname, "_bin_res_", configuration$binSize, "bp.rds")
					} else {
						file.out <- paste0(RDS.BIN.F, "merge_condition_", mergeExpname, "_bin_res_", (configuration$binSize/1e3), "kb.rds")
					}				
					bin.GR <- doBins(file=file.out,  expname=mergeExpname, bin=configuration$binSize , reads=merge, assemblyName=unique(vpInfo$genome), 
						config_genomes=configuration$genomes, report=reportAnalysis, vpInfo=vpInfo, log.path=log.path)
				}

				# Merge Plot
				if (make.cisplot == TRUE){
					message("            *** Creating local 4C Plot")
					createPlot(plotTitle=mergeExpname, vpPos=unique(vpInfo$pos), chromosome=unique(vpInfo$chr), fragGR=merge[as.vector(seqnames(merge)) == unique(vpInfo$chr)], 
						plotLegend=reportAnalysis, plotView=configuration$plotView,	maxY=configuration$maxY, minY=0, xaxisUnit=configuration$xaxisUnit, plotRegion='cis', 
						foldOut=PLOT.MERGE.F, plotType=configuration$plotType)
				}

				if (disabled == 0 & make.gwplot == TRUE){
					message("            *** Creating genome-wide 4C coverage Plot")
					createPlot(plotTitle=mergeExpname, vpPos=unique(vpInfo$pos), chromosome=unique(vpInfo$chr), fragGR=bin.GR, plotLegend=reportAnalysis, plotView=configuration$plotView,
						maxY=configuration$maxY, minY=0, xaxisUnit=configuration$xaxisUnit, plotRegion='all', foldOut=GENOMEPLOT.MERGE.F, plotType=configuration$plotType)
				}

				#Export data
				if (make.bdg == TRUE){
					message("            *** Creating bedGraph file")

					if (nonBlind){
						mergebdgFile <- paste0(BDG.MERGE.F, "merge_condition_", mergeExpname, "_nonblind_WIN",wSize,".bedGraph")
					} else {
						mergebdgFile <- paste0(BDG.MERGE.F, "merge_condition_", mergeExpname, "_WIN", wSize ,".bedGraph")
					}
					gzMergeBDGName <- paste0(mergebdgFile, ".gz")

					if (file.exists(mergebdgFile) | file.exists(gzMergeBDGName)){
						error.msg <- paste0("               --- WARNING: merge_condition_", mergeExpname, " bedGraph file already exists.")
						message(error.msg)
						write(error.msg, log.path, append=TRUE)
					} else {
						if (unique(vpInfo$analysis) == "cis"){
							mergeReads.bdg <- merge[order(merge$pos)]
						} else {
							mergeReads.bdg <- merge[order(seqnames(merge), merge$pos)]
						}
						message("               --- Generation of bedGraph file.")
						exportBDG(gR=mergeReads.bdg, expName=mergeExpname, filename=mergebdgFile, vpPos=unique(vpInfo$pos), 
							vpChr=unique(vpInfo$chr), plotView=configuration$plotView, interval=FALSE)
					}
				}

				if (make.wig == TRUE){
					message("            *** Creating WIG file")

					if (nonBlind){
						mergeWigFile <- paste0(WIG.MERGE.F, "merge_condition_", mergeExpname, "_nonblind_WIN", wSize,".wig")
					} else {
						mergeWigFile <- paste0(WIG.MERGE.F, "merge_condition_", mergeExpname, "_WIN", wSize ,".wig")
					}
					gzMergeWigName <- paste0(mergeWigFile, ".gz")
					bwMergeName <- gsub(".wig", ".bw", mergeWigFile)

					if (file.exists(mergeWigFile) | file.exists(gzMergeWigName) & file.exists(bwMergeName)){
						error.msg <- paste0("               --- WARNING: merge_condition_", mergeExpname, " WIG and BIGWIG files already exist.")
						message(error.msg)
						write(error.msg, log.path, append=TRUE)
					} else if (file.exists(mergeWigFile) | file.exists(gzMergeWigName) & !file.exists(bwMergeName)){
						error.msg <- paste0("               --- WARNING: merge_condition_", mergeExpname, " WIG file already exists.")
						message(error.msg)
						write(error.msg, log.path, append=TRUE)
						if (file.exists(mergeWigFile)){
							do.wigToBigWig(wig_file=mergeWigFile, config_genomes=configuration$genomes, assemblyName=unique(vpInfo$genome))
						} else {
							do.wigToBigWig(wig_file=gzMergeWigName, config_genomes=configuration$genomes, assemblyName=unique(vpInfo$genome))
						}				
					} else {
						if (unique(vpInfo$analysis) == "cis"){
							mergeReads.wig <- merge[order(merge$pos)]
						} else {
							mergeReads.wig <- merge[order(seqnames(merge), merge$pos)]
						}
						message("               --- Generation of variableStep WIG file.")
						exportWig(gR=mergeReads.wig, expName=mergeExpname, filename=mergeWigFile, vpPos=unique(vpInfo$pos), vpChr=unique(vpInfo$chr),
							plotView=configuration$plotView, config_genomes=configuration$genomes, assemblyName=unique(vpInfo$genome), fixed=FALSE)
					}

					if (make.fixedStep & disabled == 0){
						if (nonBlind){
							mergeWigBinFile <- paste0(WIG.MERGE.F, "merge_condition_", mergeExpname, "_nonblind_BIN_", fixedStepWigBin, ".wig")
						} else {
							mergeWigBinFile <- paste0(WIG.MERGE.F, "merge_condition_", mergeExpname, "_BIN", fixedStepWigBin, ".wig")
						}

						gzMergeWigBinName <- paste0(mergeWigBinFile, ".gz")
						bwMergeBinName <- gsub(".wig", ".bw", mergeWigBinFile)

						if (file.exists(mergeWigBinFile) | file.exists(gzMergeWigBinName) & file.exists(bwMergeBinName)){
							error.msg <- paste("               --- WARNING: ", mergeExpname, " fixedStep WIG and BIGWIG file already exist.")
							message(error.msg)
							write(error.msg, log.path, append=TRUE)
						} else if (file.exists(mergeWigBinFile) | file.exists(gzMergeWigBinName) & !file.exists(bwMergeBinName)){
							error.msg <- paste0("               --- WARNING: ", mergeExpname, " fixedStep WIG file already exists.")
							message(error.msg)
							write(error.msg, log.path, append=TRUE)
							if (file.exists(mergeWigBinFile)){
								do.wigToBigWig(wig_file=mergeWigBinFile, config_genomes=configuration$genomes, assemblyName=unique(vpInfo$genome))
							} else {
								do.wigToBigWig(wig_file=gzMergeWigBinName, config_genomes=configuration$genomes, assemblyName=unique(vpInfo$genome))
							}
						} else {
							if (fixedStepWigBin < 1e3){
								mergeBinReads.wig <- paste0(RDS.BIN.MERGE.F, "merge_condition_", mergeExpname, "_bin_res_", fixedStepWigBin, "bp.rds")
							} else {
								mergeBinReads.wig <- paste0(RDS.BIN.MERGE.F, "merge_condition_", mergeExpname, "_bin_res_", (fixedStepWigBin/1e3), "kb.rds")
							}

							bin.GR.fixedWig <- doBins(file=mergeBinReads.wig,  expname=mergeExpname, bin=fixedStepWigBin, reads=merge, assemblyName=unique(vpInfo$genome), 
							config_genomes=configuration$genomes, report=reportAnalysis, vpInfo=vpInfo, log.path=log.path)

							message("               --- Generation of fixedStep WIG file.")
							exportWig(gR=bin.GR.fixedWig$bins, expName=mergeExpname, filename=mergeWigBinFile, vpPos=unique(vpInfo$pos), vpChr=unique(vpInfo$chr), 
								plotView=configuration$plotView, config_genomes=configuration$genomes, assemblyName=unique(vpInfo$genome), fixed=TRUE, bin=fixedStepWigBin)
						}
					}
				}

				if (tsv == TRUE){
					message("            *** Creating TSV file")

					if (nonBlind){
						mergeTsvFile <- paste0(TSV.MERGE.F, "merge_condition_", mergeExpname, "_nonblind_WIN", wSize,".tsv")
					} else {
						mergeTsvFile <- paste0(TSV.MERGE.F, "merge_condition_", mergeExpname, "_WIN", wSize ,".tsv")
					}
					gzMergeTSVName <- paste0(mergeTsvFile, ".gz")

					if (file.exists(mergeTsvFile) | file.exists(gzMergeTSVName)){
						error.msg <- paste0("               --- WARNING: merge_condition_", mergeExpname, " TSV file already exists.")
						message(error.msg)
						write(error.msg, log.path, append=TRUE)
					} else {
						if (unique(vpInfo$analysis) == "cis"){
							mergeReads.tsv <- merge[order(merge$pos)]
						} else {
							mergeReads.tsv <- merge[order(seqnames(merge), merge$pos)]
						}
						exportTSV(gR=mergeReads.tsv, filename=mergeTsvFile, merge=TRUE)
					}
				}

				# PeakC analysis
				if (peakC == TRUE){
					message("            *** PeakC analysis")
					mergeBedFile = paste0(PEAKC.BED_BDG.F, "merge_condition_", mergeExpname, "_peakC_peaks.bed")
					mergeBedFileRegion = gsub("peakC_peaks", "peakC_regions", mergeBedFile)
					mergeBedGraphPeakCFile = paste0(PEAKC.BED_BDG.F, "merge_condition_", mergeExpname, "_peakC_peaks.bedGraph")
					mergeBedGraphPeakCFileRegion = gsub("peakC_peaks", "peakC_regions", mergeBedGraphPeakCFile)
					mergePeakCRdsFile = paste0(PEAKC.RDS.F, "merge_condition_", mergeExpname, "_peakC_results.rds")

					if (file.exists(mergePeakCRdsFile)){
						error.msg <- paste("               --- WARNING: merge_condition_", mergeExpname, " peakC analysis already exists, continuing with the existing RDS file.")
						message(error.msg)
						write(error.msg, log.path, append=TRUE)
						mergeResPeakC <- readRDS(mergePeakCRdsFile)
						VPpos <- mergeResPeakC$vpPos
						VPChr <- mergeResPeakC$vpChr
					} else {
						num.exp <- length(expname)
						VPpos <- unique(vpInfo$pos)
						VPChr <- unique(vpInfo$chr)

						zoom <- GRanges(seqnames=VPChr, resize(IRanges(VPpos, VPpos), width=vpRegion, fix="center"))
						vpGR <- merge[unique(queryHits(findOverlaps(merge, zoom)))]

						mergePeakCDat <- list()
						for (f in 1:num.exp) {
							column <- paste0("normReads_rep", f)
							data <- as.data.frame(vpGR[, c("pos", column)])
							colnames(data)[7] <- "reads"
							mergePeakCDat[[f]] <- data.frame(pos=data$pos, reads=data$reads)
						}

						mergeResPeakC <- suppressWarnings(combined.analysis(data=mergePeakCDat, num.exp=num.exp, vp.pos=VPpos, wSize=wSize, 
							alphaFDR=alphaFDR, qWr=qWr, minDist=minDist))

						mergeResPeakC$vpPos <- VPpos
						mergeResPeakC$vpChr <- VPChr

						if(length(mergeResPeakC$peak)>0){
							mergeResPeakC$exportPeakGR <- getPeakCPeaks(resPeakC=mergeResPeakC, min.gapwidth=0)
						} else {
							mergeResPeakC$exportPeakGR <- NULL
						}

						#Save peakC results
						saveRDS(mergeResPeakC, file=mergePeakCRdsFile)
					}

					#Plot peakC results
					plotFile = paste0(PEAKC.PLOTS.F, "merge_condition_", mergeExpname,".pdf")
					pdf(file=plotFile)
					plot_C(data=mergeResPeakC, num.exp = 1, y.min = 0, y.max = 750)
					dev.off()

					#Export bed file
					mergeBedPeaks <- exportPeakCPeaks(resPeakC=mergeResPeakC, bedFile=mergeBedFile, 
						name=paste0("merge_condition_", mergeExpname,"peakC_peaks"), desc=NULL, includeVP=TRUE, min.gapwidth=0)
					mergeBedRegions <- exportPeakCPeaks(resPeakC=mergeResPeakC, bedFile=mergeBedFileRegion, 
						name=paste0("merge_condition_", mergeExpname,"peakC_region"), desc=NULL, includeVP=TRUE, min.gapwidth=minGapwidth)

					#export bedGraph
					merge2 <- merge
					ranges(merge2) <- merge2$pos
					peakCReads <-  GRanges(seqnames=mergeBedPeaks$chr, ranges=mergeBedPeaks$start)
					mergeReadsPeakC.bdg <- merge2[unique(queryHits(findOverlaps(merge2, peakCReads)))]
					
					peakCRegions <-  GRanges(seqnames=mergeBedRegions$chr, ranges=paste0(mergeBedRegions$start, "-", mergeBedRegions$end))
					score <- c()
					for (l in 1:length(peakCRegions)) {
						tmp <- peakCRegions[l]
						peaks <- as.data.frame(merge2[unique(queryHits(findOverlaps(merge2, tmp)))])
						score <- c(score, median(peaks$norm4C))
					}
					mergeRegionPeakC.bdg <- peakCRegions
					mergeRegionPeakC.bdg$norm4C <- score
					rm(tmp, peaks, score)
					
					exportBDG(gR=mergeReadsPeakC.bdg, expName=paste0(mergeExpname, "_peakC_peaks"), filename=mergeBedGraphPeakCFile, 
						vpPos=VPpos, vpChr=VPChr, plotView=configuration$plotView, interval=FALSE)
					exportBDG(gR=mergeRegionPeakC.bdg, expName=paste0(mergeExpname, "_peakC_regions"), filename=mergeBedGraphPeakCFileRegion,
						vpPos=VPpos, vpChr=VPChr, plotView=configuration$plotView, interval=TRUE)

					#Preparing data to DESeq2 (recovery of all positive position per viewpoint)
					if (DESeq2 == TRUE){
						VP <- paste0("chr", gsub("-", ":", ViewP))
						if (is.null(list_peakCMergePeaksPos_Per_VP[[VP]]) == TRUE){
							dataReads <- as.data.frame(mergeReadsPeakC.bdg)
							list_peakCMergePeaksPos_Per_VP[[VP]] <- paste0(dataReads$seqnames, ":", dataReads$start)
						} else {
							dataReads <- as.data.frame(mergeReadsPeakC.bdg)
							dataReadsVP <- paste0(dataReads$seqnames, ":", dataReads$start)
							list_peakCMergePeaksPos_Per_VP[[VP]] <- unique(c(list_peakCMergePeaksPos_Per_VP[[VP]], dataReadsVP))
						}
						if (is.null(list_peakCMergeRegionsPos_Per_VP[[VP]]) == TRUE){
							dataRegionsReads <- as.data.frame(mergeRegionPeakC.bdg)
							list_peakCMergeRegionsPos_Per_VP[[VP]] <- paste0(dataRegionsReads$seqnames, ":", dataRegionsReads$start, "-", dataRegionsReads$end)
						} else {
							dataRegionsReads <- as.data.frame(mergeRegionPeakC.bdg)
							dataRegionsReadsVP <- paste0(dataRegionsReads$seqnames, ":", dataRegionsReads$start, "-", dataRegionsReads$end)
							list_peakCMergeRegionsPos_Per_VP[[VP]] <- unique(c(list_peakCMergeRegionsPos_Per_VP[[VP]], dataRegionsReadsVP))
						}
					}
					rm(peakCReads, mergeReadsPeakC.bdg, peakCRegions,mergeRegionPeakC.bdg)				
				}
			}
		}
		message('\n')
	}

	if (peakC == TRUE & DESeq2 == TRUE){
		message("\n------ DESeq2 analysis")
		VP <- unique(paste0("chr", vpChr, ":", vppos))
		for (r in VP) {
			message(paste0("      >>> Analysis on experiments with the viewpoint ", r," <<<"))
			ViewP <- unlist(strsplit(r, ":"))
			subVPinfo <- VPinfo[which(VPinfo$vpchr == gsub("chr","", ViewP[1]) & VPinfo$vppos == ViewP[2]),]
			cond <- unique(subVPinfo$condition)

			if (length(cond) < 2){
				error.msg <- paste0("         ### WARNING: experiment(s) with the viewpoint ", r, " present only one condition: ", cond,".")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				next
			}
			DIR.OUT <- paste0(PEAKC.DESeq2.F, ViewP[1], "-", ViewP[2], "/")
			logDirs[[r]] <- ifelse(!dir.exists(DIR.OUT), dir.create(DIR.OUT), FALSE)

			matrixPeaksName <- paste0(DIR.OUT, "matrix_peaksPeakC_", ViewP[1], "-", ViewP[2], ".txt")
			if (file.exists(matrixPeaksName)){
				error.msg <- paste0("         ### WARNING: a peakC peaks matrix file for experiments with the viewpoint ", 
					r, "  already exists, continuing with the existing matrix.")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				matrix_peaks <- read.table(file=matrixPeaksName, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.name=1)
			} else {
				message("         ### Creation of matrix with peakC all significant peaks for DESeq2 analysis")
				positive_peaks <- GRanges(sort(list_peakCPeaksPos_Per_VP[[r]]))
				matrix_peaks <- make.DESeq2matrix(subVPinfo=subVPinfo, gR_positivePeaks=positive_peaks, pathToRDS=RDS.F, region=FALSE, merge=FALSE)
				write.table(matrix_peaks, file=matrixPeaksName, sep = "\t", quote= FALSE, row.names=TRUE, col.names = TRUE)
			}

			matrixRegionsName <- paste0(DIR.OUT, "matrix_regionsPeakC_", ViewP[1], "-", ViewP[2], ".txt")
			if (file.exists(matrixRegionsName)){
				error.msg <- paste0("         ### WARNING: a peakC region matrix file for experiments with the viewpoint ", 
					r, "  already exists, continuing with the existing matrix.")
				write(error.msg, log.path, append=TRUE)
				message(error.msg)
				matrix_regions <- read.table(file=matrixRegionsName, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.name=1)
			} else {
				message("         ### Creation of matrix with peakC all significant regions for DESeq2 analysis")
				positive_regions <- reduce(GRanges(sort(list_peakCRegionsPos_Per_VP[[r]])))			
				matrix_regions <- make.DESeq2matrix(subVPinfo=subVPinfo, gR_positivePeaks=positive_regions, pathToRDS=RDS.F, region=TRUE, merge=FALSE)			
				write.table(matrix_regions, file=matrixRegionsName, sep = "\t", quote= FALSE, row.names=TRUE, col.names = TRUE)
			}

			message("         ### Perfoming DESeq2 analysis on peakC significant peaks in at least one replicate")
			do.deseq2(matrix=matrix_peaks, metadata=subVPinfo, output=DIR.OUT, expname=paste0("peakC_peaks_", r), ctrlCondition=ctrl)

			message("         ### Perfoming DESeq2 analysis on peakC significant regions in at least one replicate\n")
			do.deseq2(matrix=matrix_regions, metadata=subVPinfo, output=DIR.OUT, expname=paste0("peakC_regions_", r), ctrlCondition=ctrl)


			if (replicates ==TRUE){
				matrixMergePeaksName <- paste0(DIR.OUT, "matrix_allReplicates_peaksPeakC_", ViewP[1], "-", ViewP[2], ".txt")
				if (file.exists(matrixMergePeaksName)){
					error.msg <- paste0("         ### WARNING: a merged peakC peaks matrix file for experiments with the viewpoint ", 
						r, "  already exists, continuing with the existing matrix.")
					write(error.msg, log.path, append=TRUE)
					message(error.msg)
					matrix_mergePeaks <- read.table(file=matrixMergePeaksName, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.name=1)
				} else {
					message("         ### Creation of matrix with peakC significant peaks (in all replicates) for DESeq2 analysis")
					positive_mergePeaks <- GRanges(sort(list_peakCMergePeaksPos_Per_VP[[r]]))
					matrix_mergePeaks <- make.DESeq2matrix(subVPinfo=subVPinfo, gR_positivePeaks=positive_mergePeaks,
						pathToRDS=RDS.F, region=FALSE, merge=FALSE)
					write.table(matrix_mergePeaks, file=matrixMergePeaksName, sep = "\t", quote= FALSE, row.names=TRUE, col.names = TRUE)
				}

				matrixMergeRegionsName <- paste0(DIR.OUT, "matrix_allReplicates_regionsPeakC_", ViewP[1], "-", ViewP[2], ".txt")
				if (file.exists(matrixMergeRegionsName)){
					error.msg <- paste0("         ### WARNING: a merged peakC region matrix file for experiments with the viewpoint ", 
						r, "  already exists, continuing with the existing matrix.")
					write(error.msg, log.path, append=TRUE)
					message(error.msg)
					matrix_mergeRegions <- read.table(file=matrixMergeRegionsName, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.name=1)
				} else {
					message("         ### Creation of matrix with peakC significant regions (in all replicates) for DESeq2 analysis")
					positive_mergeRegions <- reduce(GRanges(sort(list_peakCMergeRegionsPos_Per_VP[[r]])))		
					matrix_mergeRegions <- make.DESeq2matrix(subVPinfo=subVPinfo, gR_positivePeaks=positive_mergeRegions,
						pathToRDS=RDS.F, region=TRUE, merge=FALSE)		
					write.table(matrix_mergeRegions, file=matrixMergeRegionsName, sep = "\t", quote= FALSE, row.names=TRUE, col.names = TRUE)
				}

				message("         ### Perfoming DESeq2 analysis on peakC significant peaks in all replicates")
				do.deseq2(matrix=matrix_mergePeaks, metadata=subVPinfo, output=DIR.OUT, expname=paste0("allReplicates_peakC_peaks_", r), ctrlCondition=ctrl)

				message("         ### Perfoming DESeq2 analysis on peakC significant regions in all replicates")
				do.deseq2(matrix=matrix_mergeRegions, metadata=subVPinfo, output=DIR.OUT, expname=paste0("allReplicates_peakC_regions_", r), ctrlCondition=ctrl)
			}
		}
		message("\n")
	}

	#make report from all RDS files
	message("------ Writing final report\n")
	export.report(RDS.F=RDS.F, OUTPUT.F=OUTPUT.F, logname=LOGNAME)
}
