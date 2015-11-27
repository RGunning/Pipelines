args <- '--inputfiles /Volumes/lustre110/Pipeline/Sailfish/output/run2/job0/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job1/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job2/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job3/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job4/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job5/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job6/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job7/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job8/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job9/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job10/sailfishtranscriptome/quant_bias_corrected.genes.sf /Volumes/lustre110/Pipeline/Sailfish/output/run2/job11/sailfishtranscriptome/quant_bias_corrected.genes.sf --sex F F F M M M F F F M M M --strain B B B B B B B B B B B B --cell B B B T T T T T T B B B
--fastq FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE --filenames 10966_8#3.bam 10966_8#2.bam 10966_8#1.bam 11049_5#18.bam 11049_5#17.bam 11049_5#16.bam 11049_4#6.bam 11049_4#5.bam 11049_4#4.bam 11048_5#15.bam 11048_5#14.bam 11048_5#13.bam'
hh <- paste(unlist(args),collapse=' ')

listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- lapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
})
options.names <- lapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)
print(options.args)
par(mfrow=c(1,1))
edgerfdr <- 0.001
deseq2padj <- 0.001

designmat <- cbind(options.args$filenames,options.args$strain,options.args$cell,options.args$sex,options.args$fastq)
designmat<-designmat[which(options.args$strain=="B" & options.args$fastq=="FALSE"),]
inputfiles <- options.args$inputfiles[which(options.args$strain=="B"& options.args$fastq=="FALSE")]

# Reformat designmatrix
#designmat <- t(designmat)
strain    <- factor(designmat[,2])
cell      <- factor(designmat[,3])
sex    <- factor(designmat[,4])
colData <- data.frame(strain=strain,cell=cell,sex=sex)


# Initialise blank count matrix and design matrix
file <- inputfiles[1]
x <- read.table(file,header=FALSE,as.is=TRUE)
Count.Matrix <- data.frame(Row.names=x[,1])
TPM.Matrix <- data.frame(Row.names=x[,1])
FPKM.Matrix <- data.frame(Row.names=x[,1])

for (i in 1:length(inputfiles)){
  file <- inputfiles[i]
  print(file)
  # Read table
  x <- read.table(file,header=TRUE,as.is=TRUE)
  Count.Matrix <- merge(Count.Matrix,   x[,c(1,5)], by.x= "Row.names", by.y= 1)
  TPM.Matrix <- merge(TPM.Matrix,   x[,c(1,3)], by.x= "Row.names", by.y= 1)
  FPKM.Matrix <- merge(FPKM.Matrix,   x[,c(1,4)], by.x= "Row.names", by.y= 1)
  
  names(Count.Matrix)[i+1] <- names(TPM.Matrix)[i+1] <- paste(designmat[i,1],designmat[i,2],designmat[i,3],designmat[i,4],sep = "_")
}

rownames(Count.Matrix) <- Count.Matrix$Row.names
rownames(TPM.Matrix) <- TPM.Matrix$Row.names
rownames(FPKM.Matrix) <- FPKM.Matrix$Row.names

Count.Matrix$Row.names <- TPM.Matrix$Row.names <- FPKM.Matrix$Row.names <- NULL
TPMmeanB <- rowMeans(TPM.Matrix[,which(cell=="B")], na.rm = FALSE, dims = 1)
TPMmeanT <- rowMeans(TPM.Matrix[,which(cell=="T")], na.rm = FALSE, dims = 1)
FPKMmeanB <- rowMeans(FPKM.Matrix[,which(cell=="B")], na.rm = FALSE, dims = 1)
FPKMmeanT <- rowMeans(FPKM.Matrix[,which(cell=="T")], na.rm = FALSE, dims = 1)
TPMmean <- data.frame(TPMmeanB, TPMmeanT,FPKMmeanB,FPKMmeanT)
rownames(TPMmean) <- rownames(TPM.Matrix)

LMG <- read.table("~/github/D1/LMG-table.tsv",header = T, sep="\t",as.is = F,na.strings = "\\N")
#LMG[LMG=="\\N"]<-NA
All <- merge(TPMmean,LMG,by.x='row.names',by.y='gene_id2_GM3')

B0_4m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.0.0.txt",as.is=T)
B0_36m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.0.1.txt",as.is=T)
B0_27ac<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.0.2.txt",as.is=T)
B1_4m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.1.0.txt",as.is=T)
B1_36m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.1.1.txt",as.is=T)
B1_27ac<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.1.2.txt",as.is=T)
B100_4m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.100.0.txt",as.is=T)
B100_36m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.100.1.txt",as.is=T)
B100_27ac<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/genes.100.2.txt",as.is=T)

T0_4m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.0.0.txt",as.is=T)
T0_36m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.0.1.txt",as.is=T)
T0_27ac<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.0.2.txt",as.is=T)
T1_4m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.1.0.txt",as.is=T)
T1_36m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.1.1.txt",as.is=T)
T1_27ac<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.1.2.txt",as.is=T)
T100_4m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.100.0.txt",as.is=T)
T100_36m<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.100.1.txt",as.is=T)
T100_27ac<-read.table("/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/genes.100.2.txt",as.is=T)

write.table(merge(B0_4m,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B0_4m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B0_36m,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B0_36m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B0_27ac,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B0_27ac",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B1_4m,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B1_4m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B1_36m,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B1_36m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B1_27ac,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B1_27ac",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B100_4m,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B100_4m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B100_36m,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B100_36m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(B100_27ac,All[,c(1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/B100_27ac",row.names=F,col.names=T,sep=",",quote = F)

write.table(merge(T0_4m,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),    "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T0_4m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T0_36m,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),   "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T0_36m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T0_27ac,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),  "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T0_27ac",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T1_4m,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),    "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T1_4m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T1_36m,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),   "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T1_36m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T1_27ac,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),  "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T1_27ac",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T100_4m,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),  "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T100_4m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T100_36m,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'), "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T100_36m",row.names=F,col.names=T,sep=",",quote = F)
write.table(merge(T100_27ac,All[,c(7,1:5,89:102)],by.x=1, by.y='Row.names'),"/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryT/T100_27ac",row.names=F,col.names=T,sep=",",quote = F)

library(VennDiagram)
z<-list(H3K4me3_B=B0_4m[[1]],H3K4me3_T=T0_4m[[1]],H3K36me3_B=B0_36m[[1]],H3K36me3_T=T0_36m[[1]],H3K27ac_B=B0_27ac[[1]])
#get.venn.partitions(z)
x<-calculate.overlap(z)
plot.new()
grid.draw(venn.diagram(z,NULL))

venn.diagram(list(H3K4me3_B=B0_4m[[1]],H3K4me3_T=T0_4m[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K4me3_0expressed.png",
             imagetype="png",main="H3K4me3 not expressed",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))
venn.diagram(list(H3K4me3_B=B1_4m[[1]],H3K4me3_T=T1_4m[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K4me3_1expressed.png",
             imagetype="png",main="H3K4me3 1-100TPM",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))
venn.diagram(list(H3K4me3_B=B100_4m[[1]],H3K4me3_T=T100_4m[[1]]),
              "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K4me3_100expressed.png",
              imagetype="png",main="H3K4me3 >100TPM",
              category=c("B cell","T cell"),
              scaled=TRUE,
              fill=c("Green","Blue"),
              cat.pos=c(-125,125))

venn.diagram(list(H3K36me3_B=B0_36m[[1]],H3K36me3_T=T0_36m[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K36me3_0expressed.png",
             imagetype="png",main="H3K36me3 not expressed",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))
venn.diagram(list(H3K36me3_B=B1_36m[[1]],H3K36me3_T=T1_36m[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K36me3_1expressed.png",
             imagetype="png",main="H3K36me3 1-100TPM",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))
venn.diagram(list(H3K36me3_B=B100_36m[[1]],H3K36me3_T=T100_36m[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K36me3_100expressed.png",
             imagetype="png",main="H3K36me3 >100TPM",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))

venn.diagram(list(H3K27ace3_B=B0_27ac[[1]],H3K27ace3_T=T0_27ac[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K27ace3_0expressed.png",
             imagetype="png",main="H3K27ace3 not expressed",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))
venn.diagram(list(H3K27ace3_B=B1_27ac[[1]],H3K27ace3_T=T1_27ac[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K27ace3_1expressed.png",
             imagetype="png",main="H3K27ace3 1-100TPM",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))
venn.diagram(list(H3K27ace3_B=B100_27ac[[1]],H3K27ace3_T=T100_27ac[[1]]),
             "/Volumes/lustre110/Pipeline/ChIPseqAnalysis/output/run1/job0/ngsavgcategoryB/H3K27ace3_100expressed.png",
             imagetype="png",main="H3K27ace3 >100TPM",
             category=c("B cell","T cell"),
             scaled=TRUE,
             fill=c("Green","Blue"),
             cat.pos=c(-125,125))

library(Vennerable)
z<-list(H3K4me3_B=B0_4m[[1]],H3K4me3_T=T0_4m[[1]],H3K36me3_B=B0_36m[[1]],H3K36me3_T=T0_36m[[1]],H3K27ac_B=B0_27ac[[1]],H3K27ac_T=T0_27ac[[1]])
Vstem= Venn(z)
Vstem3=Vstem[,c("H3K4me3_B","H3K4me3_T","H3K36me3_B","H3K36me3_T")]
plot(Vstem,type="battle",doEuler=TRUE)

