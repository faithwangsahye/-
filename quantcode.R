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