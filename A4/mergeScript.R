install.packages("reshape2")
#merged_bed <- read.table("merged_hg38.bed",header = FALSE,stringsAsFactors=FALSE, quote="") #read.table() will always return a data frame

merged_bed <- read.csv("C:/Users/bkirs/Documents/merged_hg38.bed",header = FALSE, sep="",stringsAsFactors=FALSE, quote="")
colnames(merged_bed) <- c('read_id','chr_hg38','start_hg38','end_hg38','strand_hg38','chr_hg19','start_hg19','end_hg19','strand_hg19') 
#read_id <- c('i1','i2','i3','i4')
#chr_hg38 <- c('1','1','1','2')
#start_hg38 <- c('10','20','30','40')
#end_hg38 <-  c('110','120','130','140')
#strand_hg38 <-  c('+','+','-','+')
#chr_hg19 <- c('1','2','3','5')
#start_hg19 <- c('10','22','30','38')
#end_hg19 <- c('100','120','130','140')
#strand_hg19 <- c('+','+','+','+')
#merged_bed <- data.frame(read_id,chr_hg38,start_hg38,end_hg38,strand_hg38,chr_hg19,start_hg19,end_hg19,strand_hg19)

read_statistics <- group_by(merged_bed, read_id)

cols2compare <- c(2,6)
#type.a <- c("hg38", "redo")
#type.b <- c("hg19", "noDet")

type.a <- "hg38"
type.b <- "hg19"

reads.per.chr <- function(merged_bed, cols2compare=c(2,6), type.a=c("hg38", "redo"), type.b=c("hg19", "noDet")){
  
  canonical_chromosomes <- paste0("chr", 1:22) #this line creates a character vector comprised of UCSC-style chromosome names (i.e. "chrN")
  # it contains chromosomes 1 through 22, but excludes chrX, chrY, chrM, and unplaced or unlocalized sequences. 
  
  chr_subset <- merged_bed[,c(cols2compare[1])]
  table_chrs1 <- table(chr_subset)
  
  chr_subset <- merged_bed[,c(cols2compare[2])]
  table_chrs2 <- table(chr_subset)
  
  
  compare.df <- data.frame(column1=table_chrs1[names(table_chrs1) %in% canonical_chromosomes],
                           column2=table_chrs2[names(table_chrs2) %in% canonical_chromosomes])
  
  ccompare.df <- compare.df[,c(1,2,4)]
  colnames(compare.df) <- c("Chr",paste0(type.a, "_reads"), paste0(type.b, "_reads"))
  
  compare.df <- melt(compare.df)
  
  return(compare.df)
  
}
