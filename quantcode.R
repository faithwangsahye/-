# GWAS
install.packages("qqman")
if(!requireNamespace("BiocManager",quietly = TRUE))
install.packages("BiocManager")
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
install.packages("qqman")

getwd()
setwd("D:/RNAseq/Experiment/quant/lab3")
workdir <- getwd()

myY <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt" ,head = FALSE)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3)

#View(gwasResults)
head(gwasResults)
str(gwasResults)
#View(snpsOfInterest)
head(snpsOfInterest)
vignette("qqman")

#绘制曼哈顿图
manhattan(gwasResults)
manhattan(gwasResults, 
          main = "Manhattan Plot", 
          ylim =c(0, 10), 
          cex =0.6, 
          cex.axis =0.9, 
          col = c("blue4", "orange3"), 
          suggestiveline = F,
          genomewideline = F,
          chrlabs =c(1:20, "P", "Q"))
#一条染色体的曼哈顿图
manhattan(subset(gwasResults, CHR == 1))
#对3号染色体上一些snp位点进行高亮标注。
manhattan(gwasResults, highlight = snpsOfInterest)
# 将高亮显示和限制在单个染色体上，并使用xlim图形参数放大感兴趣的区域(位置在200-500之间)
manhattan(subset(gwasResults, CHR == 3), highlight = snpsOfInterest, xlim =c(200,500), main = "Chr3")
# 根据snp的p值对其进行注释,只对每个超过p值阈值的染色体的顶部SNP进行注释
manhattan(gwasResults, annotatePval =0.01)                                                     
#标注所有满足阈值的snp:
manhattan(gwasResults, annotatePval =0.005, annotateTop =FALSE)
# 添加统计测验
gwasResults <- transform(gwasResults, zscore = qnorm(P/2,lower.tail =FALSE))
head(gwasResults)
# 绘制一个新的图形
manhattan(gwasResults, 
          p = "zscore",
          logp =FALSE, 
          ylab ="Z-score", 
          genomewideline =FALSE,  
          uggestiveline =FALSE, 
          main = "Manhattan plot of Z-scores")

qq(gwasResults$P)
qq(gwasResults$P,
   main = "Q-Q plot of GWAS p-values", 
   xlim = c(0, 7), 
   ylim =c(0, 12), 
   pch =18, 
   col = "blue4", 
   cex =1.5, 
   las = 1)

## ====================


getwd()
setwd("D:/RNAseq/Experiment/quant")
workdir <- getwd()

mdpdata = read.table("mdp_numeric.txt",header = T)
mdpdata[,1]

X.raw = mdpdata[,-1]
X = X.raw

#Set missing values
mr = .2
n = nrow(X)
m = ncol(X)
dp = m*n
uv = runif(dp)
#画直方图
hist(uv)
missing = uv<mr
length(missing)
#uv<0.9
missing[1:10]
index.m = matrix(missing,n,m)
dim(index.m)
X[index.m] = NA
X.raw[1:5,1:5]
X[1:5,1:5]

#Define Stochastic Impute funciton
StochasticImpute=function(X){
  
  n=nrow(X)
  m=ncol(X)
  fn=colSums(X,na.rm=T)  # sum of genotypes for all individuals
  fc=colSums(floor(X/3+1),na.rm=T) #count number of non missing individuals
  fa=fn/(2*fc) #Frequency of allele"2"
  
  for(i in 1:m){
    
    index.a=runif(n)<fa[i]
    index.na=is.na(X[,i])
    index.m2=index.a  &  index.na
    index.m0=!index.a  &  index.na
    X[index.m2,i]=2
    X[index.m0,i]=0
  }
  return(X)
}

#Impute
XI= StochasticImpute(X)

#Correlation
accuracy.r=cor(X.raw[index.m], XI[index.m])

#Proportion of match
index.match=X.raw==XI
index.mm=index.match&index.m
accuracy.m=length(X[index.mm])/length(X[index.m])
accuracy.r
accuracy.m

#Impute with StochasticImpute
XI= StochasticImpute(X)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
install.packages("impute")
library(impute)

#Imputeand calculate correlation
XI= StochasticImpute(X)
X.knn= impute.knn(as.matrix(t(X)),k=10)

accuracy.r.si=cor(X.raw[index.m], XI[index.m])
accuracy.r.knn=cor(X.raw[index.m], t(X.knn$data)[index.m])
accuracy.r.si
accuracy.r.knn

getwd()
system("java -jar beagle.22Jul22.46e.jar unphased = test.bgl out=result")

system("java -jar beagle.22Jul22.46e.jar unphased = test.bgl out=D:/RNAseq/Experiment/quant/result")
