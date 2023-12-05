import numpy as np
import sys

#数据载入和类型转成float
def makeAA531Dict(filename):
    AA531FileName = 'AA531properties.txt'
    fr = open(AA531FileName)
    arrayOLines = fr.readlines()
    del(arrayOLines[0]) # 删除head行
    fr.close()

    AA531Dict = {}
    for line in arrayOLines:
        line = line.strip()
        listFromLine = line.split('\t')
        AA = listFromLine[0]
        properties = [float(i) for i in listFromLine[1:]]
        AA531Dict[AA] = properties
    return AA531Dict

#肽序列表征
def file2matrix(filename, seqLength, AA531Dict):
    #数据载入
    AASeqFileName = 'ACEtriPeptidesSequencesActivities.txt'
    fr = open(AASeqFileName)
    arrayOLines = fr.readlines()
    fr.close()

    seqLength = 3
    lineNum = 0
    numberOfLines = len(arrayOLines) #文件行数

    '''定义一个0矩阵用来存结果，numberOfLines行，531*seqLength列'''
    returnMat = np.zeros((numberOfLines, 531*seqLength))
    Y = np.zeros((numberOfLines, 1))#只有一列的0矩阵

    for line in arrayOLines:
        line = line.strip()
        listFromLine = line.split('\t')
        AASeq = listFromLine[0]
        Y[lineNum] = float(listFromLine[1])
        
        feaVec = []
        '''遍历listFromLine，将每个氨基酸替换为相应的531个理化属性,如果序列中的氨基酸在AA531Dict中有对应的key，才进行替换,否则以0替换'''
        for AA in AASeq:
            if AA in AA531Dict.keys():
                feaVec.extend(AA531Dict[AA])
            else:
                print('Warning: nonregular amino acid found! Coding "%s" in "%s"(seqId: %d) with 531 zeros.' % (AA, AASeq, lineNum))
                feaVec.extend([0.0]*531)
                Y[lineNum] = -1
        '''把returnMat的第lineNum行 换成新的feaVec'''
        returnMat[lineNum,:] = np.array(feaVec)
        lineNum += 1
    return Y, returnMat, lineNum

#将结果写入文件
if __name__ == '__main__':
    AASeqFileName = sys.argv[1]
    AA531FileName = sys.argv[2]
    seqLength = int(sys.argv[3])
    outputFileName = sys.argv[4]
    AA531Dict = makeAA531Dict(AA531FileName)
    Y, AA531Mat, SeqNum = file2matrix(AASeqFileName, seqLength, AA531Dict)

outputFileName = 'result.txt'
np.savetxt(outputFileName, np.hstack((Y, returnMat)), fmt='%g', delimiter='\t')
print('The number of sequences is %d. Matrix of features is saved in %s' % (lineNum, outputFileName))