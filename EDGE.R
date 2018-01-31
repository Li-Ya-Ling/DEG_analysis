#!/usr/local/bin/Rscript
library('edgeR')
library('getopt')

# getting parameters 获取输入参数
spec = matrix(c(
        'help',     'h',  0, "logical",
        'infile',   'i',  1, "character",
        'outdir',   'o',  0, "character",
        'sample1',  's1', 1, "character",
        'sample2',  's2', 1, "character",
        'bcv',      'b',  0, "numeric",
    ),
    byrow=TRUE, ncol=4,
);
opt = getopt(spec);

# 设置参数提示信息
print_usage <- function(spec=NULL) {
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat(
        "Usage example:
        Rscript EDGE.R --infile A_VS_B.count.txt --outdir od --sample1 A --sample2 B

        Options:
        --help          -h	NULL		get this help
        --infile        -i	character	the input file [forced]
        --outdir        -o	character	the output dir [forced]
        --sample1       -s1	character	the sample 1 name; or the samples of group1, such as "Sample1_Sample2_Sample3"
        --sample2       -s2	character	the sample 2 name; or the samples of group2, such as "Sample4_Sample5_Sample6"
	--bcv           -b      numeric		Biological cofficient of variation, while only 1 sample in each group, (0,1), Default: 0.1
        \n"
    )
    q(status=1);
}

# 设置 bcv 默认值, 0.1，通常如果是实验控制的好的人类数据，那么选择BCV=0.4，比较好的模式生物选择BCV=0.1，技术重复的话选择BCV=0.1。
if ( is.null(opt$bcv)) { opt$bcv = 0.1 }

# 获取脚本所在目录
getScriptPath <- function() {
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}
scr_path <- getScriptPath()

# 加载脚本所在目录下的 dge_plot_funcs2.r 
plot2_source <- paste (scr_path,"/dge_plot_funcs2.r",sep="")
source (plot2_source)

# 读取输入文件
ALL<- read.delim(as.vector(opt$infile), row.names = 1, header=TRUE, sep="\t")

# 获取差异分析样本分组信息
v1<-strsplit(as.character(opt$sample1),"_")
v2<-strsplit(as.character(opt$sample2),"_")
v1<-as.vector(v1[[1]])
v2<-as.vector(v2[[1]])

# 提取差异分析样本对应的数据
counts <- data.frame(ALL[,c(v1,v2)])
keep<-rownames(counts)
counts<-counts[keep,]
rownames(counts)<-rownames(ALL[keep,])

# 计算rpkm表达量
geneLength<-ALL[keep,length(ALL)]
rpkm<-rpkm(counts, gene.length = geneLength)
rownames(rpkm)<-rownames(counts)



et<-0
if (length(v1) == 1 | length(v2) == 1) {
# 两组内均只有一个样本时
    # 生成edgeR能够识别的object
    y <- DGEList(counts=counts, group=1:2)
    # 寻找差异表达基因
    et <- exactTest(y, dispersion=opt$bcv^2)
}
else {
    # 设置分组
    group<-factor(c(rep(1, length(v1)), rep(2, length(v2))))
    # 生成edgeR能够识别的object
    y <- DGEList(counts=counts, group=group)
    # 对表达值进行标准化以及估计disperse值
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    # 寻找差异表达基因
    et <- exactTest(y)
}

# 进行fdr矫正
result <- topTags(et,dim(et)[1])
# 获取上/下调信息
x <- data.frame(result[,c(-2)])
regulated<-''
for(i in 1:dim(x)[1]) {
    if (x[i,1]>0) regulated[i]<-'up'  else regulated[i]<-'down'
}
# 生成表达量结果数据
x<-x[,c(2,3,1)]
x<-cbind(x,regulated)
counts=counts[row.names(counts)%in%row.names(x),]
counts = round(counts/(rep(1, times = nrow(counts)) %*% t(colSums(counts)))*1000000,2)
# 生成输出的结果数据
sx <- rownames(counts)
x<-cbind(counts[sx,],x[sx,])
x=data.frame("ID"=row.names(x),x)
colnames(x)=c("#ID",colnames(x)[-1])
colnames(x)[ncol(x)-1]="log2FC"
# 输出两组样本间的差异表达分析结果，结果包括每个样本的基因表达量，每个基因的差异表达倍数，pvalue，FDR，上下调信息。
out1=paste(opt$outdir,"/",opt$sample1,"_vs_",opt$sample2,".final.xls",sep="")
write.table(x,file=out1,row.names = FALSE,sep = "\t", col.names = TRUE,quote=F)
