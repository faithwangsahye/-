#2023.10.20-2023.12.5
#FaithHuangSiHui实验课记录

##实验一##

data = read.table("D:/RNAseq/Experiment/Biostat/1/习题data.txt")
d2 = matrix(as.matrix(data),ncol = 1)
bi = seq(25.5,47.5,by = 3)
table = hist(d2,breaks = bi,plot = F)
table$counts#不要图，将直方图转换成列表

hist(d2,breaks = bi,axes = F)# 绘制次数分布图，先不绘制坐标轴
axis(2,at = seq(0,60,5))#y轴
axis(1,at = seq(20,50,by = 5))#x轴

library(tidyverse)
library(ggplot2)

df_count <- data.frame(
  "组限" = c("25.5-28.5","28.5-31.5","31.5-34.5","34.5-37.5","37.5 40.5","40.5-43.5","43.5-46.5")
)
df_count$"中点值" <- seq(27,45,by = 3)
df_count$"次数" <- table$counts

df_count1 <- data.frame(
  "组限" = c("25.5-28.5","28.5-31.5","31.5-34.5","34.5-37.5","37.5 40.5","40.5-43.5","43.5-46.5")
)
df_count1$"次数" <- table$counts

plot(df_count$中点值 , df_count$次数,type = "b")

#算术平均
summary(data)
summary(data$V1)
mean(data$V1)
#标准差
sd(data$V1)
#标准误
sd(data$V1)/sqrt(5)
#变异系数
sd(data$V1)/mean(data$V1)*100
#二项分布概率
y<- 0:4
dbinom(y,10,0.4)#（for i in y）A事件发生的概论为0.4，做10次实验，发生i次的概论
#二项分布累积概率
pbinom(y,10,0.4)

x <- c(88, 78, 67, 69, 62, 100, 73, 45, 70, 60, 93, 97, 84, 82, 81, 73, 68, 76, 77, 92)
tmp<-table(x) #计算出x中每个值出现的次数
index<-which.max(tmp) #找出最多次数的索引
tmp[index] #输出对应的数据及次数

#习题4.7
#(1)
a <- c(2.38,2.41,2.50,2.47,2.41,2.38,2.26,2.32,2.41)
t.test(a,mu = 2.5)
b1 <- pnorm(20,16,2) - pnorm(10,16,2)
#result : 0.9759

#(2)
b2 <- pnorm(12,16,2)
#result:0.0228
b3 <- 1 - pnorm(20,16,2)
#result:0.0228

#(3)
b4 <- qnorm(0.75,16.2)
b5 <- qnorm(0.25,16.2)
#15.52551~16.87449

#(4)
qnorm(0.975,16.2)
qnorm(0.025,16.2)
#14.240041~18.15996

##DONE##




##实验二##
#单个样本平均数 的假设测验
#u测验
#已知某地区的当地水稻原品种平均株高为100cm, 但总体方差未知；现测量了60株【生成60个随机数】水稻新品种(由R随机产生, 样本均数和方差容易算得), 问新品种所属总体与原品种总体之间是否有显著差异？
a <- round(rnorm(60, mean = 120, sd = 15))
u <- (mean(a)-100) / (sd(a)/sqrt(60))
u
#计算确切的P值
p1 <- pnorm(-u)*2
p1

#t测验
weight <- c(35.6,37.6,33.4,35.1,32.7,36.8,35.9,34.6)
t.test(weight, mu=34)#mu是数据的真实均值
#RESULT：p-value = 0.07485 > 0.05 ，确认是零假设u=u0=34

#两个样本平均数相比较 的假设测验
#t测验
y1 <- c(400,420,435,460,425)
y2 <- c(450,440,445,445,420)
t.test(y1,y2,var.equal = T)
#var.equal=TURE指定样本之间是等方差的，alternative=指定单侧检验
# 配对t检验
t.test(y1,y2,paired=TRUE) 
#p-value = 0.3126 > 0.05，接受零假设，没有显著差异
#二项资料：单个样本百分数
#二项资料：两个样本百分数相比较

##DONE##




##实验三##
#1.1

#因素：
grow <- c(42.9,45.4,46.4,42.8,41.5,46.5,47.0,
          42.5,43.2,41.1,43.1,41.6,42.9,
          19.9,20.3,21.5,24.4,23.7,21.5,21.8,23.41)
#水平 重复数分别为：7，6，8
pz <- factor(rep(c("1", "2", "3"), times = c(7, 6, 8)))
#方差分析表
result1 <- aov(grow ~ pz)
summary(result1) #F value：407.4
#差异绘图
plot(grow ~ pz)
#多重比较结果
pairwise.t.test(grow, pz, p.adjust.method = "bonferroni")


#1.2
#比较资料转换前后方差分析的差别
#转换前
per <- c(0.8,4.0,9.8,6.0,
         3.8,1.9,56.2,79.8,
         0.0,0.7,66.0,7.0,
         6.0,3.5,10.3,84.6,
         1.7,3.2,9.2,2.8)
pz1 <- factor(rep(1:4, times = 5))

result2 <- aov(per ~ pz1)
summary(result2)
#反正弦转换后：
#百分化-正平方根
scaper <- sqrt(per/100)
chanper <- asin(scaper)

result3 <- aov(chanper ~ pz1)
summary(result3)

#2
#指出哪些是因素【稻田和水样】、每个因素有几个水平【各3个水样】、
#每个水平下有几次重复【每水样分析2次】, 列出方差分析表、多重比较结果
y <- c(1.1,1.3,1.2,
       1.2,1.1,1.0,
       1.3,1.3,1.4,
       1.4,1.5,1.2,
       1.8,2.1,2.2,
       2.0,2.0,1.9)
daotian <- gl(3,6,labels = c("A1","A2","A3"))
shuiyang <- rep(gl(3,1,labels = c("B1","B2","B3")),6) 
data.frame(daotian, shuiyang, y)
#差异图
op <- par(mfrow = c(1,2))
plot(y ~ daotian + shuiyang)
#方差分析表
result4 <- aov(y ~ daotian + shuiyang + daotian:shuiyang)
summary(result4)
#关于稻田因素的多重比较结果
pairwise.t.test(y, daotian, p.adjust.method = "bonferroni")

##DONE##




##中心极限##

#原始总体 population
normal_data <- rnorm(n, mean, sd)
uniform_data <- runif(n, min, max)
binomial_data <- rbinom(n, size, prob)  # 成功概率
poisson_data <- rpois(n, lambda)  # lambda平均发生率
t_data <- rt(n, df) # df自由度

sample_size <- c(1,2,4,10,20,50)

num_iterations <- 6 # 迭代次数 
for (i in 1:num_iterations){
  size <- sample(sample_size,1)
  sample(population , size = size , replace = TRUE)
  cat("n = ", i, ": ", resampled_data, "\n")
}



dataa <- rnorm(100, 68, 7)

sample_size <- c(1,2,4,10,20,50)
num_iterations <- 6 # 迭代次数 

for (i in 1:num_iterations){
  size <- sample(sample_size,1
                 )
  resampled_data <- sample(dataa , size = size , replace = TRUE)
  cat("n = ", i, ": ", resampled_data, "\n")
}





##实验四-卡方##

#适合性测验：第 1 题、第 2 题、第 6 题
#chisq.test 函数
#矫正：chi2.adjusted <- sum((abs(O-E)-0.5)^2/E)
#O 观察值，p 理论值/理论比例，E 期望值E <- sum(O) * p
#1-pchisq(chi2.adjusted,1)

chisq.test(c(3437,3482),p=c(0.5,0.5))
O <- c(3437,3482)
p <- c(0.5,0.5)
E <- sum(O) * p
E
chi2.adjusted <- sum((abs(O-E)-0.5)^2/E)
chi2.adjusted
1-pchisq(chi2.adjusted,1)

#【1】
adjustment_test <- function(O,P){
  E <- sum(O) * P
  chi2.adjusted <- sum((abs(O-E)-0.5)^2/E)
  print(chi2.adjusted)
  adj <- 1-pchisq(chi2.adjusted,1)
  print(adj)
}
#(1)
chisq.test(c(134,36), p=c(3/4,1/4))
adjustment_test(O = c(134, 36), P = c(3/4, 1/4))
#(2)
chisq.test(c(240,120), p=c(3/4,1/4))
adjustment_test(O = c(240,120), P = c(3/4, 1/4))
#(3)
chisq.test(c(76,56), p=c(0.5,0.5))
adjustment_test(O = c(76,56), P=c(0.5,0.5))
#(4)
chisq.test(c(240,13), p=c(15/16,1/16))
adjustment_test(O = c(240,13), P = c(15/16,1/16))

#【2】
chisq.test(c(348,115,157),p=c(9/16,3/16,4/16))

#【6】
chisq.test(c(132,42,38,14),p=c(9/16,3/16,3/16,1/16))


#独立性测验：第 5 题、第 7 题、第 8 题
#chisq.test(four_table) chisq.test 函数默认对四格表进行了连续性矫正
#2×C 表独立性检验
table23 <- matrix(c(29,68,96,22,199,2),nrow = 2, byrow = T)
chisq.test(table23)
#R×C 表独立性检验
table33 <- matrix(c(146,7,7,183,8,13,152,14,16),nrow = 3, byrow = T)
chisq.test(table33)

#【5】
apple <- matrix(c(150,168,14,26),nrow = 2, byrow =T)
apple
chisq.test(apple)

#【7】
plants <- matrix(c(142,51,3,13,404,2,2,17,176),nrow = 3, byrow =T)
chisq.test(plants)

#【8】
potatoes <- matrix(c(19,75,12,76,6,108,15,107),nrow = 4, byrow =T)
potatoes
chisq.test(potatoes)

##DONE##




##实验四：常用试验设计及方差分析##
#第十一章第1题，第3题，第5题，第6题，第9题；第十二章，第3、4、7题。

#【11.5】
y<-c(24.2,21.4,26.4,26.4,20.4,
     21.0,26.9,19.6,26.8,25.8,
     17.7,24.7,26.7,23.5,25,
     25.8,25.5,21.8,20.7,23.6,
     22.4,22.5,25.6,29.3,24.4)
fertilizer2 <- factor(c("B","E","A","C","D",
                          "D","A","E","B","C",
                          "E","B","C","D","A",
                          "A","C","D","E","B",
                          "C","D","B","A","E"))
row_block <- c(rep(1:5,each=5))
col_block <- c(rep(1:5,5))
data5 <- data.frame(y, fertilizer2)
data5$col_block <- factor(col_block) #col列标
data5$row_block <- factor(row_block) #row 行标
head(data5)
fit5 <- aov(y~ fertilizer2+col_block+row_block, data=data5)
summary(fit5)

#【11.6】单因素随机区组设计
#四个区组四次重复
data6 <- c(6.2,5.8,7.2,5.6,6.9,7.5,
           6.6,6.7,6.6,5.8,7.2,7.8,
           6.9,6.0,6.8,5.4,7.0,7.3,
           6.1,6.3,7.0,6.0,7.4,7.6)
varities <- rep(gl(6,1,
                   labels = c("A","B","C","D","E","F")),4)
blocks <- factor(rep(1:4,each=6))
xiaomai <- data.frame(varities,blocks,data6)
fit <- aov(data6~varities+blocks)
summary(fit)


data6ps <- c(6.2,6.6,6.9,6.1,
           5.8,6.7,6.0,6.3,
           7.2,6.6,6.8,7.0,
           5.6,5.8,5.4,6.0,
           6.9,7.2,7.0,7.4,
           7.5,7.8,7.3,7.6)
varities2 <- gl(6,4,labels = c("A","B","C","D","E","F"))
blocks2 <- factor(rep(c(1:4),6))
xiaomai2 <- data.frame(data6ps,varities2,blocks2)
xiaomai2
fit2 <- aov(data6ps ~ varities2 + blocks2 + varities2:blocks2)
summary(fit2)


#【11.9】
x119 <- c(10,12,17,14,12,
          8,13,15,14,10,
          6,8,13,17,10,
          8,11,11,15,16)
y119 <- c(18,36,40,21,42,
          17,38,36,23,36,
          14,28,35,24,38,
          15,30,29,20,52)
blocks119 <- gl(4,5,labels = c(1:4))
varities119 <- factor(rep(c("A","B","C","D","E"),4))
#x变数
yumi119x <- data.frame(x119,blocks119,varities119)
yumi119x
fit119x <- aov(x119 ~ varities119 + blocks119)
summary(fit119x)
#y变数
yumi119y <- data.frame(y119,blocks119,varities119)
yumi119y
fit119y <- aov(y119 ~ varities119 + blocks119)
summary(fit119y)

cov(yumi119x, yumi119y)


#【12.3】
y2 <- c(12,13,14,15,13,16,14,13,16,12,14,14,
        16,14,14,15,12,13,16,13,13,15,13,17,
        13,15,11,14,17,14,12,15,15,13,15,13)
pinzhong <- factor(c("a1","a2","a3","a4","a2","a4","a3","a1","a4","a1","a3","a2",
                     "a4","a1","a2","a3","a1","a2","a4","a3","a2","a3","a1","a4",
                     "a2","a3","a1","a2","a4","a3","a2","a4","a3","a1","a4","a1"))

bozhong <- factor(c("b1","b2","b3","b2","b1","b3","b2","b3","b1","b2","b1","b3",
                    "b2","b3","b1","b3","b2","b3","b1","b2","b2","b1","b1","b3",
                    "b3","b1","b2","b1","b3","b2","b2","b1","b3","b3","b2","b1"))

quzu <- c(rep(1:3,each=12))
col_block <- c(rep(1:12,3))

data113 <- data.frame(y, bozhong)
data5$col_block <- factor(col_block)
data5$row_block <- factor(row_block)
head(data5)

##实验六##
#【8.4】
x <- c(0,5,10,15,20,25,30)
y <- c(0.00,0.11,0.23,0.34,0.46,0.57,0.71)
data1 <- data.frame(x,y)
data1
lm(y~x,data = data1)
plot(data1$x, data1$y,
     xlab = "density", ylab = "corelation")
abline(fit)
fit <- lm(y~x,data = data1)
summary(fit)

#[8.5]
x <- c(13,25,27,23,26,1,15)
y <- c(50,55,50,47,51,29,48)
data2 <- data.frame(x,y)
r <- cor(data2$x, data2$y)
lm(y ~ x, data = data2)
fit <- lm(y ~ x, data = data2)
summary(fit)
sqrt(sum(residuals(fit)^2)/(7-2))
data3 <- data.frame(x = 5, y = 37.0279)
predict(fit, newdata = data2, 
        interval = "prediction", 
        level = 0.95)
