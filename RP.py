#!/usr/bin/env python
# coding: utf-8
#FaithHuangSihui 模式识别与预测项目实现
# # LAB 2

# In[ ]:


import numpy as np # 导入numpy包，并重命名为np

def file2matrix(filename, KMAX, bpTable):#filename是要处理的文件名，KMAX是最大的k值，bpTable是一个空的碱基对频率表。
    fr = open(filename) # 打开文件
    arrayOLines = fr.readlines() # 读取所有内容
    fr.close() # 及时关闭文件

    numberOfLines = len(arrayOLines) # 得到文件行数
    returnMat = np.zeros((numberOfLines, 16*(KMAX+1))) # 为返回的结果矩阵开辟内存
    lineNum = 0

    for line in arrayOLines:
        line = line.strip() # 删除空白符，包括行尾回车符
        listFromLine = line.split(': ') # 以': '为分隔符进行切片
        nt_seq = list(listFromLine[1]) # 取出核酸序列并转换成list
        del(nt_seq[70:72]) # 删除位于第71，72位的供体位点
        
        kSpaceVec = []

        for k in range(KMAX+1): # 计算不同k条件下的kSpace特征
            bpFreq = bpTable.copy() # bpTable是一个字典型变量，一定要用字典的copy函数.计算不同k值（从0到KMAX）下的碱基对频率（bpFreq）

            for m in range(len(nt_seq)-k-1): # 扫描序列，并计算不同碱基对的频率。迭代核酸序列(nt_seq)，更新bpFreq字典中的频率
                bpFreq[nt_seq[m]+nt_seq[m+1+k]] += 1 # 序列的子串会自动在字典中寻找对应的key，很神奇！否则要自己写if语句匹配
            bpFreqVal = list(bpFreq.values()) # 取出bpFreq中的值并转换成list
            kSpaceVec.extend(np.array(bpFreqVal)/(len(nt_seq)-k-1)) # 每个k下的特征，需除以查找的所有子串数。bpFreqVal中的值除以总的子串数，以获得频率，计算得到的频率（bpFreqVal）添加到kSpaceVec列表中

        returnMat[lineNum,:] = kSpaceVec
        lineNum += 1
    return returnMat, lineNum

if __name__ == '__main__':
    filename = 'EI_true1.seq'
    KMAX = 4
    bpTable = {}
    for m in ('A','T','C','G'):
        for n in ('A','T','C','G'):
            bpTable[m+n] = 0

    kSpaceMat, SeqNum = file2matrix(filename, KMAX, bpTable)
    outputFileName = 'EI_true1_kSpace.txt'
    np.savetxt(outputFileName, kSpaceMat, fmt='%g', delimiter=',')
    print('The number of sequences is %d. Matrix of features is saved in %s' % (SeqNum, outputFileName))


# # LAB 3

# #### 对之前的代码进行了一些修改和添加，可以从命令行接收参数，并在处理过程中输出进度信息。

# In[ ]:


import numpy as np # 导入numpy包，并重命名为np
import sys # 导入sys包，用于从命令行传递参数给python程序

def file2matrix(filename, bpTable, KMAX=2): # 为KMAX提供默认参数(updated)
    fr = open(filename) # 打开文件
    arrayOLines = fr.readlines() # 读取所有内容
    del(arrayOLines[:4]) # 删除头4行（updated, 避免了运行程序之前，另外使用sed删除头4行）
    fr.close() # 及时关闭文件
    
    numberOfLines = len(arrayOLines) # 得到文件行数
    returnMat = np.zeros((numberOfLines, 16*(KMAX+1))) # 为返回的结果矩阵开辟内存
    
    lineNum = 0
    for line in arrayOLines:
        line = line.strip() # 删除空白符，包括行尾回车符
        listFromLine = line.split(': ') # 以': '为分隔符进行切片
        nt_seq = list(listFromLine[1]) # 取出核酸序列并转换成list
        del(nt_seq[70:72]) # 删除位于第71，72位的供体位点
        
        kSpaceVec = []
        for k in range(KMAX+1): # 计算不同k条件下的kSpace特征
            bpFreq = bpTable.copy() # bpTable是一个字典型变量，一定要用字典的copy函数，Python函数参数使用的址传递

            for m in range(len(nt_seq)-k-1): # 扫描序列，并计算不同碱基对的频率
                sub_str = nt_seq[m]+nt_seq[m+1+k] # 提出子串(updated)
                if sub_str in bpFreq.keys(): # 如果子串在bpFreq中有对应的key，才统计频次(updated, NOTE:在供体虚假位点序列中存在非正常碱基)
                    bpFreq[sub_str] += 1 # 序列的子串会自动在字典中寻找对应的key，很神奇！否则要自己写if语句匹配
            bpFreqVal = list(bpFreq.values()) # 取出bpFreq中的值并转换成list
            kSpaceVec.extend(np.array(bpFreqVal)/(len(nt_seq)-k-1)) # 每个k下的特征，需除以查找的所有子串数

        returnMat[lineNum,:] = kSpaceVec#将kSpaceVec赋值给returnMat矩阵的当前行（lineNum）
        lineNum += 1#追踪当前处理的行数
        if (lineNum % 1000) == 0:
            print('Extracting k-spaced features: %d sequences, done!' % lineNum)#如果lineNum是1000的倍数，输出进度信息。
    return returnMat, lineNum

if __name__ == '__main__':
    filename = sys.argv[1]
    outputFileName = sys.argv[2]
    KMAX = int(sys.argv[3])
    bpTable = {}
    for m in ('A','T','C','G'):
        for n in ('A','T','C','G'):
            bpTable[m+n] = 0
    
    kSpaceMat, SeqNum = file2matrix(filename, bpTable, KMAX)
    np.savetxt(outputFileName, kSpaceMat, fmt='%g', delimiter=',')
    print('The number of sequences is %d. Matrix of features is saved in %s' % (SeqNum, outputFileName))


# ## 构建训练集、测试集

# #### 命令行参数1、trueSiteFileName 2、falseSiteFileName 3、trainingSetFileName 4、testSetFileName

# In[ ]:


import numpy as np
import sys
from random import sample # 导入sample函数，用于从虚假位点数据中随机抽取样本
from sklearn.model_selection import train_test_split # 用于产生训练集、测试集

trueSiteFileName = sys.argv[1]#将命令行传递给脚本的第一个参数赋值给trueSiteFileName变量
falseSiteFileName = sys.argv[2]

trueSitesData = np.loadtxt(trueSiteFileName, delimiter = ',') # 载入true位点数据
numOfTrue = len(trueSitesData)
falseSitesData = np.loadtxt(falseSiteFileName, delimiter = ',') # 载入false位点数据
numOfFalse = len(falseSitesData)

randVec = sample(range(numOfFalse), len(trueSitesData)) # 随机产生true位点样本个数的随机向量。sample函数生成一个随机向量randVec，长度与真实位点数据的数量相同
falseSitesData = falseSitesData[randVec,] # 以随机向量从false位点数据中抽取样本

Data = np.vstack((trueSitesData, falseSitesData)) # 按行将true位点与false位点数据组合
Y = np.vstack((np.ones((numOfTrue,1)),np.zeros((numOfTrue,1)))) # 包含了真实位点数据和虚假位点数据对应的标签（1表示真实位点，0表示虚假位点）
testSize = 0.3 # 测试集30%，训练集70%
X_train, X_test, y_train, y_test = train_test_split(Data, Y, test_size = testSize, random_state = 0)

trainingSetFileName = sys.argv[3]
testSetFileName = sys.argv[4]
testSetFileName = sys.argv[4]
#将训练集数据和标签组合后保存到训练集文件中。数据和标签在水平方向上进行组合，即标签作为第一列，数据作为后续列。
np.savetxt(trainingSetFileName, np.hstack((y_train, X_train)), fmt='%g', delimiter=',') # 将Y与X以列组合后，保存到文件
np.savetxt(testSetFileName, np.hstack((y_test, X_test)), fmt='%g', delimiter=',')
print('Generate training set(%d%%) and test set(%d%%): Done!' % ((1-testSize)*100, testSize*100))


# ## KNN

# In[ ]:


import numpy as np
from sklearn import neighbors # 导入KNN包
import sys

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集，在命令行指定文件名
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

n_neighbors = int(sys.argv[3]) # 在命令行指定邻居数
weights = 'uniform' # 每个邻居的权重相等
clf = neighbors.KNeighborsClassifier(n_neighbors, weights=weights) # 创建一个KNN的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of KNN: %g%% (%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))


# ## LR

# In[ ]:


import numpy as np
from sklearn import linear_model # 导入线性模型包
import sys

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

maxIterations = int(sys.argv[3]) # 在命令行指定最大迭代次数
clf = linear_model.LogisticRegression(max_iter=maxIterations) # 创建一个LR的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of LR: %g%% (%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))


# ## DT

# In[ ]:


import numpy as np
from sklearn import tree # 导入Decision Trees包
import sys
import graphviz # 导入Graphviz包

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

clf = tree.DecisionTreeClassifier() # 创建一个DT的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of DT: %g%% (%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))

# Export the tree in Graphviz format
graphFileName = sys.argv[3] # 从命令行指定图文件名称
dotData = tree.export_graphviz(clf, out_file=None)
graph = graphviz.Source(dotData)
graph.render(graphFileName)
print('The tree in Graphviz format is saved in "%s.pdf".' % graphFileName)


# ## NB

# In[ ]:


import numpy as np
from sklearn import naive_bayes # 导入NB包
import sys

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集，在命令行指定文件名
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

clf = naive_bayes.GaussianNB() # 创建一个NB的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of NB: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))


# ## SVM svc

# #### 包含网格寻优

# In[ ]:


import numpy as np
from sklearn import svm # 导入svm包
import sys
from sklearn import preprocessing # 导入数据预处理包
from sklearn.model_selection import GridSearchCV # 导入参数寻优包
from random import sample

#数据载入
train = np.loadtxt(sys.argv[1], delimiter=',') 
test = np.loadtxt(sys.argv[2], delimiter=',') 

# 从train中随机抽200样本用于后续建模
train = train[sample(range(len(train)), 200),] 
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

isScale = int(sys.argv[3]) # 建模前，是否将每个特征归一化到[-1,1]
kernelFunction = sys.argv[4] # {‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’}, default=’rbf’
isChooseCG = int(sys.argv[5]) # 是否寻找最优参数：c, g

if isScale:
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(-1,1))
    trX = min_max_scaler.fit_transform(trX)
    teX = min_max_scaler.transform(teX)

if isChooseCG:
    numOfFolds = int(sys.argv[6]) # 命令行指定寻优过程的交叉验证次数
    C_range = np.power(2, np.arange(-5,15,2.0)) # 指定C的范围
    gamma_range = np.power(2, np.arange(3,-15,-2.0)) # 指定g的范围
    parameters = dict(gamma=gamma_range, C=C_range) # 将c, g组成字典，用于参数的grid遍历
    
    clf = svm.SVC(kernel=kernelFunction) # 创建一个SVC的实例
    grid = GridSearchCV(clf, param_grid=parameters, cv=numOfFolds) # 创建一个GridSearchCV实例
    grid.fit(trX, trY) # grid寻优c, g
    print("The best parameters are %s with a score of %g" % (grid.best_params_, grid.best_score_))
    clf = svm.SVC(kernel=kernelFunction, C=grid.best_params_['C'], gamma=grid.best_params_['gamma'])
else:
    clf = svm.SVC(kernel=kernelFunction)#如果不需要寻找最优参数，则直接创建一个SVC实例，并指定核函数为命令行参数中指定的值。
    
clf.fit(trX, trY) # 训练模型
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of SVC: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))


# ## 运行时间较长的将命令写到脚本中再用qsub提交任务到集群
#!/bin/bash
#$ -S /bin/bash
#$ -N mySVC
#$ -j y
#$ -cwd

# SVC分类器：在命令行指定训练集、测试集，规格化，线性核，参数寻优，10次交叉
echo '------ scale: 1; kernel: linear; chooseCG: 1; numOfCV: 10 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 linear 1 10
echo

# 规格化，线性核，不参数寻优
echo '------ scale: 1; kernel: linear; chooseCG: 0 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 linear 0
echo

# 规格化，径向基核(rbf)，参数寻优，10次交叉
echo '------ scale: 1; kernel: rbf; chooseCG: 1; numOfCV: 10 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 rbf 1 10
echo

# 规格化，径向基核(rbf)，不参数寻优
echo '------ scale: 1; kernel: rbf; chooseCG: 0 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 rbf 0
echo
# # LAB 5

# In[75]:


import numpy as np
import sys, os

workdir = 'D:\RNAseq\Experiment\PatternR&P\QSAR'

#原始数据载入
#AA531FileName = 'AA531properties.txt'
fr = open('D:\RNAseq\Experiment\PatternR&P\QSAR\AA531properties.txt')
arrayOLines = fr.readlines()
del(arrayOLines[0]) # 删除head行
fr.close()

#数据类型转成float
AA531Dict = {}
for line in arrayOLines:
    line = line.strip()
    listFromLine = line.split('\t')
    AA = listFromLine[0]
    properties = [float(i) for i in listFromLine[1:]]
    AA531Dict[AA] = properties


# In[76]:


#原始数据载入
#AASeqFileName = 'ACEtriPeptidesSequencesActivities.txt'
fr = open('D:\RNAseq\Experiment\PatternR&P\QSAR\ACEtriPeptidesSequencesActivities.txt')
arrayOLines = fr.readlines() 
fr.close()

seqLength = 3
lineNum = 0

numberOfLines = len(arrayOLines)
#定义一个0矩阵用来存结果，numberOfLines行，531*seqLength列
returnMat = np.zeros((numberOfLines, 531*seqLength))
Y = np.zeros((numberOfLines, 1))#只有一列的0矩阵

for line in arrayOLines:
        line = line.strip() # 删除空白符，包括行尾回车符
        listFromLine = line.split('\t') # 以'\t'为分隔符进行切片
        AASeq = listFromLine[0] # 取出氨基酸序列
        Y[lineNum] = float(listFromLine[1]) # 取出活性值Y
        
        feaVec = []
        for AA in AASeq: # 遍历listFromLine，将每个氨基酸替换为相应的531个理化属性
            if AA in AA531Dict.keys(): # 如果序列中的氨基酸在AA531Dict中有对应的key，才进行替换
                feaVec.extend(AA531Dict[AA])
            else: # 否则以0替换
                print('Warning: nonregular amino acid found! Coding "%s" in "%s"(seqId: %d) with 531 zeros.' % (AA, AASeq, lineNum))
                feaVec.extend([0.0]*531)
                Y[lineNum] = -1
        #把returnMat的第lineNum行 换成新的feaVec
        returnMat[lineNum,:] = np.array(feaVec)
        lineNum += 1


# In[77]:


outputFileName = 'result.txt'
np.savetxt(outputFileName, np.hstack((Y, returnMat)), fmt='%g', delimiter='\t')
print('The number of sequences is %d. Matrix of features is saved in %s' % (lineNum, outputFileName))

import numpy as np
import sys

# 1. 将AA531properties.txt做成字典
# def makeAA531Dict(filename):
AA531FileName = 'AA531properties.txt'
    fr = open(AA531FileName) # 打开文件
    arrayOLines = fr.readlines() # 读取所有内容
    del(arrayOLines[0]) # 删除head行
    fr.close() # 及时关闭文件

    AA531Dict = {}
    for line in arrayOLines:
        line = line.strip()
        listFromLine = line.split('\t')
        AA = listFromLine[0]
        properties = [float(i) for i in listFromLine[1:]] # 从文件读取的数值默认是字符串类型，需要转换为浮点型
        AA531Dict[AA] = properties
    #return AA531Dict

# 2. 肽序列表征
# def file2matrix(filename, seqLength, AA531Dict):
AASeqFileName = 'ACEtriPeptidesSequencesActivities.txt'
seqLength = 3
    fr = open(AASeqFileName) # 打开文件
    arrayOLines = fr.readlines() # 读取所有内容
    fr.close() # 及时关闭文件

    numberOfLines = len(arrayOLines) # 得到文件行数
    returnMat = np.zeros((numberOfLines, 531*seqLength)) # 为返回的结果矩阵开辟内存
    Y = np.zeros((numberOfLines, 1))
    lineNum = 0

    for line in arrayOLines:
        line = line.strip() # 删除空白符，包括行尾回车符
        listFromLine = line.split('\t') # 以'\t'为分隔符进行切片
        AASeq = listFromLine[0] # 取出氨基酸序列
        Y[lineNum] = float(listFromLine[1]) # 取出活性值Y
        
        feaVec = []
        for AA in AASeq: # 扫描序列，将每个氨基酸替换为相应的531个理化属性
            if AA in AA531Dict.keys(): # 如果序列中的氨基酸在AA531Dict中有对应的key，才进行替换
                feaVec.extend(AA531Dict[AA])
            else: # 否则以0替换
                print('Warning: nonregular amino acid found! Coding "%s" in "%s"(seqId: %d) with 531 zeros.' % (AA, AASeq, lineNum))
                feaVec.extend([0.0]*531)
                Y[lineNum] = -1

        returnMat[lineNum,:] = np.array(feaVec)
        lineNum += 1
    #return Y, returnMat, lineNum

# 3. 将结果写入文件
#if __name__ == '__main__':
#    AASeqFileName = sys.argv[1]
#    AA531FileName = sys.argv[2]
#    seqLength = int(sys.argv[3])
#    outputFileName = sys.argv[4]
#    AA531Dict = makeAA531Dict(AA531FileName)
#    Y, AA531Mat, SeqNum = file2matrix(AASeqFileName, seqLength, AA531Dict)

outputFileName = 'result.txt'
np.savetxt(outputFileName, np.hstack((Y, returnMat)), fmt='%g', delimiter='\t')
print('The number of sequences is %d. Matrix of features is saved in %s' % (lineNum, outputFileName))