getwd()
setwd("D:/RNAseq/Experiment/R_experiment")
workdir <- getwd()

##实验1##

# 1、标量，向量，矩阵，数组，有序因子，数据框，列表
sca <- 10
vec <- c(1,2,3)
matr <- matrix(1:9,nrow = 3,ncol = 3)#三行三列
arr <- array(1:8,dim = c(2,2,2))#两组两行两列
of <- ordered(c("A","B"), levels = c("A", "B"))
View(of)
data_frame <- data.frame(
  name = c("J", "A", "B"),
  age = c(18, 19, 20),
  height = c(170, 165, 180)
)#按列的方向进行填充
list <- list(
  name = "J",
  age = 18,
  scores = c(80, 90, 85)
)
list$name
data_frame$name
#list和df格式都可以用$进行数据选择

# 2、导入数据
data = read.table("D:/RNAseq/Experiment/R_experiment/gene.xlsx")
data1 = read.table("D:/RNAseq/Experiment/R_experiment/manager.txt")
View(data)

# 3、从不同的数据结构取值并赋给新变量
datav1 <- data$V1
name <- data_frame$name
age <- list$age

# 4、seq和rep的用法
#seq()生成序列,rep()重复向量
seq(1,10,by = 2)
seq(1,10,length.out = 2)
rep(1:10,times = 2)#整体重复
rep(1:3, each = 3) #分别依次重复

# 5、保存文件
write.table(data_frame, file.path(workdir, 'dataframe'), sep='\t')


##实验2##

# 1、读入数据,以第一列为行名-header=TRUE
ex2_data = read.table("D:/RNAseq/Experiment/R_experiment/manager.txt",
                      header=TRUE )

# 2、将99定义为NA
ex2_data[5,5] <- NA

# 3、创建新变量
ex2_data$agecat <- ifelse(ex2_data$age>75,"elder",
                          ifelse(ex2_data>55,"middle","young"))

# 4、数据框重命名
library(dplyr)
ex2_data <- rename(ex2_data,birthday = date)
View(ex2_data)
  #查看所有列名
names(ex2_data)
colnames(ex2_data)
colnames(ex2_data) <- c("manager","birthday","country","gender",
                        "age","one","two","three","q4","q5","agecat")
  #查看所有行名
rownames(ex2_data)
rownames(ex2_data) <- c("a","b","c","d","e")
  #重命名行名or列名
row.names(probe_data)[5] <- "Rename_as"

# 5、日期值输出
weekdays(as.Date(2004-05-7))

# 6、产生新变量：平均数、方差
probe_data <- read.table("D:/RNAseq/Experiment/R_experiment/exp_probe.csv",
                         sep = "," , 
                         head = T , 
                         row.name = "probe_ID")
probe_data$average <- apply(probe_data[2:85],1,mean)
probe_data$varience <- apply(probe_data[2:85],1,var)
  #降序重排：根据variance，从大到小降序排列
ord <- order(probe_data$varience,decreasing = T)    #是否降序？Ture
  #对行排列，则空出列位置
probe_data_ord <- probe_data[ordered,]    
view(probe_data_ord)
#apply函数,直接按照平均数排序
apply(probe_data, 1,mean)

# 7、产生数据框新变量sum（分数总和）和average（平均分）
data2 = read.table("D:/RNAseq/Experiment/R_experiment1/data.txt",head = T)
#产生dataframe新变量
probe_data$sum <- Sum(matrix(1:20 , 4))
probe_data$average <- average(matrix(1:20 , 4))


##实验3##

getwd()
setwd("D:/RNAseq/Experiment/R_experiment/R_experiment3")
workdir <- getwd()
data1 = read.table(paste0(workdir, "/manager.txt"), sep="\t", header=T)

#(一)	数据集合并、子集选取
#1、剔除变量q3 与q5
data1 = data1[,!(names(data1) %in% c("q3","q5"))]
data1 = data1[,which(!(names(data1) %in% c("q3","q5")))]

#2、选择所有35岁以上的男性
data2 = data1[which(data1$age >= 35 & data1$gender == c("M")),]
data2 = subset(data1 , data1$age >= 35 & data1$gender == c("M"))
data2$q4 <- data2$gender
data2

#3、数据框格式读入，数据集列合并
exp_probe1 = read.table(paste0(workdir,"/exp-probe.csv"), sep=",", header=T)
exp_probe2 = read.table(paste0(workdir,"/exp-probe.2.csv"), sep=",", header=T)
probe_gene = read.table(paste0(workdir,"/probe-gene.csv"), sep=",", header=T)

gene.probe <- merge(exp_probe1,probe_gene)
library(dplyr)
probe.exp = bind_rows(exp_probe1,exp_probe2)

#(二)	数学、统计和字符处理函数
#1、
data1 = read.table(paste0(workdir, "/manager.txt"), sep="\t", header=T)
#2、
mean(data1$q1)
sd(data1$q1,na.rm = TRUE)#标准差
median(data1$q1)
var(data1$q1,na.rm = TRUE)#方差

#3、分位数
quantile(data1$q2,0.95)
quantile(data1$q2,0.99)

#4、最大值、最小值
min(data1$q5)
max(data1$q5)
  #返回各自的位置
which.min(data1$q5)
which.max(data1$q5)

#5、累积和，累积乘积
cumsum(data1$age)
cumprod(data1$age)

#6、滞后差分
diff(data1$q1)

#7、标准化
scale(data1$q3)
mean((scale(data1$q3)*5)+10)

#8、返回每行q1与q2 的最大值、最小值
data1[which.max(data1$q1),]
data1[which.min(data1$q1),]
data1[which.max(data1$q2),]
data1[which.min(data1$q2),]

#9、标准正态分布x为1.96左侧的曲线下面积
#mean=0,标准差=1的正态分布当x=1.96左侧的累积分布函数
pnorm(1.96 , mean = 0 , sd = 1 , lower.tail = T)
#均值为500，标准差为100的正态分布的0.9分为点值
#概率密度函数
dnorm(0.9 , mean = 500 , sd = 100)

#10、生成50个均值为50，标准差为10的正态随机数
rnorm(50 , mean = 50 , sd = 10)

#11、产生两个一模一样的正态分布随机数
set.seed(1)
rnorm(1,0,1)

#12、字符串替换、分隔
Mystr = c("hello,every,body,how,are,you") 
Mystrsplit = strsplit(Mystr,",")
sub("body" , "rookie" , Mystrsplit)
gsub("body" , "rookie" , Mystrsplit)
library(stringr)
str_replace(Mystrsplit,"body" , "rookie")#分隔

#13、对任意小数，均保留指定的小数点后有效位数
options(digits = 5)
a <- pi
a


##实验四##
install.packages("sampling")
library(sampling)
library(MASS)
data("Insurance")
Insurance
nrow(Insurance)
mydata <- 1:999
#1.	随机抽样sample
sam1 <- sample(nrow(Insurance),10,replace = T)
sam2 <- sample(mydata,10,replace = T)

#2.	分层抽样strata
#"srswor"【不放回】,"srswr"【有放回】,"poisson"【泊松分布】,"systematic"【系统抽样】
#size指定每个层的抽样大小
strata(Insurance,stratanames = "District",
       size = c(2,3,1,2),
       method = "srswor")

#3.	整群抽样cluster
#size用于指定抽样群组数量
cluster(Insurance, 
        clustername = "District", 
        size = 2, 
        method = "srswor", 
        description=T)

#1.	水仙花数是指一个三位数，其各位数字的立方和等于该数本身。如153=13 + 53 + 33，请统计水仙花数的个数;
#a（a ！= 0）【(1:9)】、b、c【(0:9)】、num、count
#for a,b,c in c(1:9),c(0:9),c(0:9),if a^3 + b^3 + c^3 = num,count ++

count = 0
for (a in 1:9) {
  for (b in 0:9) {
    for (c in 0:9) {
      num = a*100 + b*10 + c
      judg = a^3 + b^3 + c^3
      if(num == judg){
        count = count +1
      }
    }
  }
}
print(count)


#2.	请用循环计算表达矩阵exp_probe.csv中每个探针的平均数;
setwd("D:/RNAseq/Experiment/R_experiment/R_experiment4")
exprob <- read.table("exp-probe.csv",header = T, sep = ",",row.names = 1)
View(exprob)
nrow(exprob)#600个探针
exprob[1,]
rowMeans(exprob[1,])
for(i in 600){
  means <- rowMeans(exprob[i,])
  exprob$means <- means
}
print(exprob)
exprob$means


#3.	随机产生35bp的DNA序列，将T替换为U;
# 生成随机DNA序列
jianji <- c("A", "T", "C", "G")
samp <- sample(jianji, 35, replace = TRUE)
samp
#替换
for (i in 1:35){
  if(samp[i] == "T"){
    samp[i] <- "U"
  }
}

randomseq_result <- paste(samp,collapse = "")
randomseq_result


##实验5##
getwd()
setwd("D:/RNAseq/Experiment/R_experiment")
workdir <- getwd()

install.packages("tidyverse")
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)

data <- mtcars
data$wt  #weight - 横坐标
data$disp  #displacement - 气泡大小
data$mpg  #miles per gallon - 纵坐标
data$cyl  #cylinders - 颜色
#View(data)

ggplot(data = data) + 
  geom_point(mapping = aes(x = wt, 
                           y = mpg, 
                           color = factor(cyl), 
                           size = disp
                           ), 
             alpha = 1/2) +
  scale_color_discrete() +
  labs(
    title = "Auto mileage by weight and horsepower",
    subtitle = "Motor Trend US Magazine(1973-74 models)",
    x = "Weight（1000lbs）",
    y = "Miles/(US) gallon",
    color = "Cylinders",
    size = "Engine displacement"
  )



Dose <- c(20,30,40,45,60)
y1 <- c(16,20,27,40,60)
y2 <- c(15,18,25,31,40)
#主标题字体为斜体，绿色，增大3倍；其它字体增大2倍；Y1为红色，Y2为蓝色
par(cex.main = 3, col.main = "green", font.main = 3,
    cex.lab = 2, cex.axis = 1)
plot(Dose, y1, type = "b", 
     col = "red",pch = 15, 
     ylim = c(0, 100),
     xlim = c(0,80),
     main = "Clinical Trials for two Drugs", 
     xlab = "Dose", ylab = "Drug Responses",
     xaxt = "n")
# 定义横坐标标签和刻度
axis(side = 1, at = Dose, labels = c("Dose1", "Dose2", "Dose3", "Dose4", "Dose5"))
#添加第二条线
lines(Dose,y2,type = "b", 
      col = "blue",pch = 17,lty = 2)
# 添加图例
legend("topleft", inset = 0.05, legend = c("y1", "y2"), title = "Drug type",
       col = c("red", "blue"), 
       lty = c(1, 2), pch = c(15, 17), xpd = T)
#title("Drug type", cex.main = 0.8, col.main = "black", adj = 0, line = -1)
#添加下标题
title(xlab = "Fig1.This is hypothetical data", line = 4.5)



【1】
ball_function <- function(m,n){
  round <- c()
  for (i in 1:n) {
    counts <- 0
    turn <- 0
    while (counts <= m) {
      ball <- sample(c("blue","red"), 1, T)
      if(ball == "blue"){
        counts <- counts + 1
      }
      turn <- turn + 1
    }
    print(paste("第", i, "次需要抽取:", turn))
    round <- c(round,turn)
  }
  average <- sum(round)/n
  
return(average)
}
ball_function(100,5)

【2】
jiou_function <- function(array1,array2){
  
  ji <- 0
  ou <- c()
  
  for (i in 1:length(array1)) {
    if(array1[i] %% 2 == 1){
      ji <- ji +1
    }
    else{
      ou <- c(ou,array1[i])
    }
    i <- i + 1
  }
  
  for (i in 1:length(array2)) {
    if(array2[i] %% 2 == 1){
      ji <- ji +1
    }
    else{
      ou <- c(ou,array2[i])
    }
    i <- i + 1
  }
  
  print(ou)
  return(ji)
}

jiou_function(array1 = array(6:9), array2 = array(1:5))

【3】

randomlist <- sample(0:9,20,T)
#length(randomlist)

for (i in 1:length(randomlist)-1) {
  for (u in (i+1):length(randomlist)) {
    if(u <= length(randomlist) && randomlist[i] > randomlist[u]){
      switch <- randomlist[i]
      randomlist[i] <- randomlist[u]
      randomlist[u] <- switch
    }
  }
}
print(randomlist)
