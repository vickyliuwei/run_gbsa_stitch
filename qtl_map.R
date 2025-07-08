# 安装必要包
if (!require("vcfR")) install.packages("vcfR")
if (!require("qtl")) install.packages("qtl")
if (!require("tidyverse")) install.packages("tidyverse")

library(vcfR)
library(qtl)
library(tidyverse)

# 1. 读取VCF文件
vcf <- read.vcfR("F2_population.vcf")

# 2. 提取基因型信息
gt <- extract.gt(vcf, element = "GT")

# 3. 数据过滤
# 过滤缺失率过高的标记
missing_rate <- apply(gt, 1, function(x) sum(x == "./.")/length(x))
gt_filtered <- gt[missing_rate < 0.2, ]

# 过滤次要等位基因频率(MAF)过低的标记
maf <- apply(gt_filtered, 1, function(x) {
  alleles <- unlist(strsplit(x, "/"))
  alleles <- alleles[alleles != "."]
  if(length(alleles) == 0) return(NA)
  min(table(alleles))/length(alleles)
})
gt_filtered <- gt_filtered[maf > 0.2 & !is.na(maf), ]

# 4. 确定亲本基因型 (假设前两个样本是亲本P1和P2)
p1 <- gt_filtered[, 1]
p2 <- gt_filtered[, 2]

# 筛选亲本间多态的标记
polymorphic <- which(p1 != p2 & p1 %in% c("0/0", "1/1") & p2 %in% c("0/0", "1/1"))
gt_poly <- gt_filtered[polymorphic, ]

# 5. 转换基因型格式为qtl包所需格式
# 创建交叉对象
cross <- list()
cross$geno <- list()
cross$pheno <- data.frame(id = colnames(gt_poly)[-c(1,2)])

# 转换基因型为AB格式 (A=亲本P1型, B=亲本P2型)
convert_geno <- function(x, p1, p2) {
  ifelse(x == p1, "a", 
         ifelse(x == p2, "b", "h"))
}

geno_data <- matrix("-", nrow = nrow(gt_poly), ncol = ncol(gt_poly)-2)
for(i in 1:4031) {
  geno_data[i, ] <- convert_geno(f_geno[i,], p1[i], p2[i])
}

# 添加染色体和位置信息
chr_pos <- str_split_fixed(rownames(gt_poly), "_", 2)
chr <- chr_pos[, 1]
pos <- as.numeric(chr_pos[, 2])

# 构建遗传图谱数据
for(ch in unique(chr)) {
  idx <- which(chr == ch)
  cross$geno[[ch]] <- list(data = t(geno_data[idx, ]),
                           map = pos[idx])
  class(cross$geno[[ch]]) <- "A"
}

class(cross) <- c("f2", "cross") 
cross <- calc.genoprob(cross, step=1, error.prob=0.001)

# 6. 构建连锁群
# 计算标记间的重组率
cross <- est.rf(cross)

# 形成连锁群
lg <- formLinkageGroups(cross, max.rf=0.35, min.lod=6)
cross <- formLinkageGroups(cross, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)

# 7. 估计遗传距离
cross <- est.map(cross, error.prob=0.001)
out.hk <- scanone(mydata,method = "hk")
# 8. 可视化结果
# 绘制连锁图
plotMap(cross, show.marker.names=FALSE)

# 绘制重组率矩阵
plotRF(cross)

# 输出遗传图谱
map <- pull.map(cross)
write.csv(map, "genetic_map.csv")
file_list <- dir("D:\\jy\\qtl_map\\bin\\")
for(i in file_list){
  data <- fread(paste0("D:\\jy\\qtl_map\\bin\\",i))
  data <- data[-c(1,2),]
  chr <- c("",rep(1,ncol(data)-1))
  chr <- t(as.data.frame(chr))
  chr <- as.data.frame(chr)
  geno <- rbind(data[1,],chr)
  geno <- rbind(geno,data[c(2:nrow(data)),])
  geno[1,1] <- "id"
  name <- gsub(".txt","",i)
  write.table(geno,file = paste0("D:\\jy\\qtl_map\\bin\\",name,".csv"),row.names = F,quote = F,col.names = F,sep = ",")
}
file_list <- dir("D:\\jy\\qtl_map\\bin\\")
nf <- file_list[grep(".csv$",file_list)]
library(qtl)
library(qtl2)
for(i in nf){
  mydata <- read.cross(format="csvs", genfile=paste0("D:\\jy\\qtl_map\\bin\\",i),phefile = "D:\\jy\\qtl_map\\phe.csv", genotypes = c("a", "b", "h"),crosstype="f2")
  mydata_qtl2=convert2cross2(mydata)
  pr <- calc_genoprob(mydata_qtl2, error_prob=0.0001, cores=4,quiet=F)
  kinship=calc_kinship(pr)
  out_pg <- scan1(pr, mydata_qtl2$pheno, kinship)
  operm <- scan1perm(pr, mydata_qtl2$pheno,kinship=kinship, n_perm=1000, cores=10)
  name <- gsub(".csv","",i)
  h <- summary(operm,alpha = 0.05)
  h <- as.numeric(h)
  png(paste0("D:\\jy\\qtl_map\\bin\\qtlbin\\",name,".png"))
  plot(out_pg, mydata_qtl2$gmap, lodcolumn=1, col="tomato",)
  abline(h = h, col = "red", lty = 2, lwd = 2)
  dev.off()
  c <- out_pg[,1]
  p <- which(c >= h)
  if(length(p) != 0){
    outlier <- rownames(out_pg)[p]
    write.table(outlier,file = paste0("D:\\jy\\qtl_map\\bin\\qtlbin\\peak",name,".txt"),col.names = F,row.names = F,sep = "\t",quote = F)
  }
}

 

