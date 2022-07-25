
#sudo docker run --rm --memory=30g --cpus=4 --mount type=bind,src=/home/dnanexus/,dst=/home/dnanexus/ -i -t rticode/deseq2:1.26.0_08a6163
#setwd("/home/dnanexus/")

#gunzip -c freeze6a.chr22.ssOnly.MAF0.001.MAC20.flt.vcf.gz | head -n 10000 | grep "#CHROM" | sed 's/#//' > header


#install.packages("survival")
library("survival")

args <- commandArgs(TRUE)
loop <- TRUE
while (loop) {
	if (args[1] == "--inFile") {
		#vcf
		inFile = args[2]
		#inFile="freeze6a.chr22.ssOnly.MAF0.001.MAC20.flt.vcf.gz" 
	}
	if (args[1] == "--chr") {
		chr = args[2]
		#chr = 22
	}
	if (args[1] == "--pheno") {
		pheno = args[2]
		#pheno = "ssOnly_n1563_coxAge_10PCs_26Aug2020.txt" 
	}

	if (length(args) > 1) {
		args = args[2:length(args)]
	} else {
		loop=FALSE
	}
}
setwd("/home/dnanexus/")

print(paste("inFile =",inFile))
print(paste("chr =",chr))
print(paste("pheno =",pheno))

print("reading in pheno and geno files")

master <- read.table(pheno,sep="\t",header=T,stringsAsFactors = F)
#tmp_vcf<-readLines(gzfile(inFile))
#gunzip -c freeze6a.chr22.ssOnly.MAF0.001.MAC20.flt.vcf.gz | head -n 10000 | grep -P "#CHROM" | sed 's/^#//' > header
tmp_vcf_data<-read.table(gzfile(inFile), stringsAsFactors = FALSE)

#get colnames from tmp_vcf
header<-read.table("header",sep="\t",header=F,stringsAsFactors=F)
#vcf_names<-unlist(strsplit(header,"\t"))
names(tmp_vcf_data)<-header[1,] 


# cox regression function
cox.me <- function(x) {
	res.cox <- coxph(Surv(cox.age, Stroke_final) ~ EV1 + 
										EV2 + EV3 + EV4 + EV5 + 
										EV6 +EV7 + EV8 + EV9 + EV10
										+ x, data=df)
	summary(res.cox)$coefficients[11,]
}
#genotype converter
# 0|0 = 0
# 0|1 or 1|0 = 1
# 1|1 = 2
convert.genotype <- function(x) {
	g <- ifelse(x == "0|0",0, 
			ifelse(x == "1|0",1,
				ifelse(x == "0|1",1,
					ifelse(x == "1|1",2,NA))))
	return(g)
}


# load in a chunk of vcf (rows=snp, col=samps)
# convert to 0,1,2 and merge with pheno
# run cox

chunk_size <- 1000
chunks<-ceiling(nrow(tmp_vcf_data)/chunk_size)

out<-NULL
for (i in 0:(chunks-1) ) {
	b <- i * chunk_size + 1
	e <- b + chunk_size
	if (e > nrow(tmp_vcf_data) ) {
		e <- nrow(tmp_vcf_data)
	}
	pct <- ceiling((i/chunks)*100)
	print(paste("chunk",i,"of",chunks,"(",pct,"%)"))
	
	# convert genotypes and save in temp matrix
	g <- as.data.frame(apply(tmp_vcf_data[b:e,10:ncol(tmp_vcf_data)],1,convert.genotype))
	af=gsub("AF=","",sapply(strsplit(tmp_vcf_data[b:e,"INFO"],";"),function(l) l[[4]]))
	af=gsub("\\.","d",af)
	colnames(g)<-paste0(tmp_vcf_data$CHROM[b:e],":",tmp_vcf_data$POS[b:e],":",tmp_vcf_data$REF[b:e],":",tmp_vcf_data$ALT[b:e],":",af)
	g$NWD<-rownames(g)
	#add alt freq
	
	# merge with pheno
	df<-merge(master,g,by="NWD")
	
	# cox regression on chunk
	res <- t(sapply(df[,14:ncol(df)],cox.me))
	out<-rbind(out,res)
}

out<-as.data.frame(out)
out$CHR <- sapply(strsplit(rownames(out),"\\."),function(l) l[[1]])
out$POS <- sapply(strsplit(rownames(out),"\\."),function(l) l[[2]])
out$REF <- sapply(strsplit(rownames(out),"\\."),function(l) l[[3]])
out$ALT <- sapply(strsplit(rownames(out),"\\."),function(l) l[[4]])
out$ALT_AF <- sapply(strsplit(rownames(out),"d"),function(l) l[[2]])
out$ALT_AF<-paste0("0.",out$ALT_AF)



colnames(out) <- c("BETA","HR","SE_BETA","Z","P","CHR","POS","REF","ALT","ALT_AF")
out<-out[,c("CHR","POS","REF","ALT","ALT_AF","BETA","HR","SE_BETA","Z","P")]
outFile = paste0("cox_results_ssOnly.chr",chr,".txt")
write.table(out, file=outFile, sep=" ", col.names=T, row.names=F, quote=F)