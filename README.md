# 专业实验课记录
1、生物统计学Biostatistics
2、R语言
## 模式识别与预测Pattern Recognize and Predicate
### 序列表征/数值化：以剪接位点识别为例
剪接位点识别
权重矩阵模型(WMM)、加权数组方法、GeneSplicer、NNsplice
极度不平衡的分类任务
数据：HS(3)D数据集：所有虚假剪接位点序列中，随机抽取2796/2880条虚假供体/受体位点序列，与所有真实剪接位点序列构建正负样本均衡的数据集。
DNA序列的k-spaced碱基对组分特征表征/数值化：如KMAX = 4，分别计算k = 0, 1, ..., KMAX时的组分特征。对于任意一条核酸序列，其k-spaced碱基对组分特征维数为：16x(KMAX+1)
### 剪接位点识别分类器：KNN, Logistic Regression, Decision Tree、NB、SVC
训练集与测试集：供体真实/虚假位点序列表征。供体真实位点与虚假位点序列的k-space组分特征。
序列表征/数值化：以定量构效关系QSAR建模为例
数据：
- ACE_tri-peptides_150数据集的肽序列及其活性(lg(IC(50)))
- ACEtriPeptidesSequencesActivities.txt
- 符合条件的生理生化属性共531个
- 每条ACE三肽序列1593个特征
- AA531properties.txt
### QSAR
从分子结构中提取、总结分子结构的信息与规律
有机小分子（如抑制剂等）与生物大分子（如受体、酶等）间的相互作用；可用于高效生物活性分子化合物筛选、药物的环境毒性评价、新药设计与合成
### ACE抑制剂
具有降高血压活性的多肽，没有降压药物的毒副作用
ACE抑制肽的ACE抑制能力：与其分子质量、与其氨基酸序列、其立体空间构象 之间存在高度相关性
ACE抑制肽的抑制类型：与ACE抑制活性、构效关系存在一定相关性
