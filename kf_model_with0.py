# -* coding: utf-8 -*-
"""
Created on Fri Dec  2 21:38:05 2022

@author: 81411
"""


import warnings
import collections
warnings.filterwarnings("ignore")
import argparse
import json
import matplotlib
import matplotlib.pyplot as plt
import math
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from torch import cuda
import sys, os
import random
import numpy as np
from sklearn import metrics
import gc
import csv
from scipy.stats import pearsonr
import scipy 
import scipy.sparse as sp
import logging
from six.moves import xrange
from collections import OrderedDict
import sys
import pdb
from sklearn import metrics
import torch.nn.functional as F
from torch.autograd import Variable

parser = argparse.ArgumentParser(description='RNN_LSTM')
parser.add_argument('--lr', type=float, default=0.0001, help='initial learning rate') #初始学习率
parser.add_argument('--model_name', type=str, default='RNN_LSTM', help='model of prediction') 
parser.add_argument('--clip', type=float, default=1,help='gradient clipping')
parser.add_argument('--epochs', type=int, default=100, help='upper epoch limit')
parser.add_argument('--batch_size', type=int, default=10, help='')
parser.add_argument('--dropout', type=float, default=0.5, help='dropout applied to layers (0 = no dropout) if n_layers LSTM > 1')
parser.add_argument('--tissue', type=str, default='Ileum', help='Tissue name')
parser.add_argument('--species', type=str, default='pig', help='Animal Species')
parser.add_argument('--save_root', type=str, default='/BIGDATA2/scau_xlyuan_1/CJL/keti/2022.11test/res/model_res/RNN_LSTM', help='where to save')
parser.add_argument('--data_root', type=str, default='/BIGDATA2/scau_xlyuan_1/CJL/keti/2022.11test/res/test20230107', help='data location')
parser.add_argument('--gpuid', type=int, default=0, help='CUDA gpu')
parser.add_argument('--gpu', type=int, default=0, help='CUDA gpu')
parser.add_argument('--n_features', type=int, default=5, help='number of histone modifications')
parser.add_argument('--n_bins', type=int, default=100, help='number of bins')
parser.add_argument('--bin_rnn_size', type=int, default=32, help='bin rnn size')  #单次训练用的样本数,通常为2^N,如32、64、128
parser.add_argument('--num_layers', type=int, default=1, help='number of layers')
parser.add_argument('--unidirectional', action='store_true', help='bidirectional/undirectional LSTM')
parser.add_argument('--save_attention_maps',action='store_true', help='set to save validation beta attention maps')
parser.add_argument('--attentionfilename', type=str, default='beta_attention.txt', help='where to save attnetion maps')
parser.add_argument('--test_on_saved_model',action='store_true', help='only test on saved model')
parser.add_argument('--num',default='p348', help='sample number')
args = parser.parse_args()

torch.manual_seed(1)
#model名称设置
model_name = '' #填模型名称
model_name += (args.species)+('_')+(args.tissue)+('_')+(args.num)+('_')
model_name+=args.model_name
args.bidirectional=not args.unidirectional
print(model_name)
#数据处理
#基因表达量字典
  
def ExpDict(gene_exp):
	# get expression value of each gene from cell*.expr.csv
	exp_dict={}
	with open(gene_exp) as fi:
		for line in fi:
			geneID,geneExpr=line.split(',')
			exp_dict[str(geneID)]=float(geneExpr)
	fi.close()
	return(exp_dict)

#组蛋白修饰数据字典
def loadData(features_file,windows,exp_dict):
	with open(features_file) as fi:
		csv_reader=csv.reader(fi)
		data=list(csv_reader)
		ncols=(len(data[0]))  #列数
	fi.close()
	nrows=len(data)
	ngenes=nrows/windows
	nfeatures=ncols-1
	print("Number of genes: %d" % ngenes)
	print("Number of entries: %d" % nrows)
	print("Number of features: %d" % nfeatures)

	count=0
	attr=collections.OrderedDict()

	for i in range(0,nrows,windows):
		hm1=torch.zeros(windows,1)
		hm2=torch.zeros(windows,1)
		hm3=torch.zeros(windows,1)
		hm4=torch.zeros(windows,1)
		hm5=torch.zeros(windows,1)
	#	hm6=torch.zeros(windows,1)
		for w in range(0,windows):
			hm1[w][0]=float(data[i+w][1])
			hm2[w][0]=float(data[i+w][2])
			hm3[w][0]=float(data[i+w][3])
			hm4[w][0]=float(data[i+w][4])
			hm5[w][0]=float(data[i+w][5])
		#	hm6[w][0]=int(data[i+w][6])
		geneID=str(data[i][0].split("_")[0])

		attr[count]={
			'geneID':geneID,
		    	'expr':exp_dict[geneID],
		    	'hm1':hm1,
		    	'hm2':hm2,
			'hm3':hm3,
            		'hm4':hm4,
            		'hm5':hm5,
            	#	'hm6':hm6,
		}
		count+=1

	return attr

class FeaturesData(Dataset):
	# Dataset class for loading data
	def __init__(self,attr,transform=None):  #经loadData函数处理后，得到的细胞数据，包含了基因，表达量和5种组蛋白的数据
		self.features=attr   #attr的数据
	def __len__(self):
		return len(self.features)    ## 返回attr数据集的行数
	def __getitem__(self,i):   ## 这里i可以理解为从数据集中提取一个基因的结果（基因ID，表达量，5种组蛋白的修饰数据，共两百个bin）
		final_data=torch.cat((
			self.features[i]['hm1'],
                        self.features[i]['hm2'],
                        self.features[i]['hm3'],
                        self.features[i]['hm4'],
                        self.features[i]['hm5']),1)
                  #      self.features[i]['hm6']),1) # torch.cat(inputs, dim=1为横拼) → Tensor
		label=self.features[i]['expr']
		geneID=self.features[i]['geneID']  ## 根据attr数据集，得到基因ID
		sample={'geneID':geneID,     ## 基因ID
			'Features':final_data,  ## 细胞1中，5种组蛋白修饰数据拼接后的矩阵
			'expr':label,         ## 即getlable中的log_fold_change=math.log((fold_change),2)，基因表达量差
               }  ## 即getlabel函数中的label1=math.log((float(c1)+1.0),2)和label2=math.log((float(c2)+1.0),2)
		return sample 

def load_data(args):
	'''
	Loads data into a 3D tensor for each of the 4 splits.

	'''
	print("==>loading train data")
	gene_dict=ExpDict(args.data_root+"/"+"last_exp_tpm.csv")
	train_dict=loadData(args.data_root+"/"+"training_tpm.csv",args.n_bins,gene_dict)
	train_inputs = FeaturesData(train_dict)
	print("==>loading valid data")
	valid_dict=loadData(args.data_root+"/"+"validation_tpm.csv",args.n_bins,gene_dict)
	valid_inputs = FeaturesData(valid_dict)
	print("==>loading test data")
	test_dict=loadData(args.data_root+"/"+"test_tpm.csv",args.n_bins,gene_dict)
	test_inputs = FeaturesData(test_dict)
	print("==>loading all gene data")
	all_dict=loadData(args.data_root+"/"+"all_gene_feature_tpm.csv",args.n_bins,gene_dict)
	all_inputs = FeaturesData(all_dict)
#gene dict 使用的是一样的，表达量
#batch_size：每次输入数据的行数，suffer：迭代训练时是否将数据洗牌
#对数据进行 batch 的划分,在训练模型时使用到此函数，用来把训练数据分成多个小组 .
	Train = torch.utils.data.DataLoader(train_inputs, batch_size=args.batch_size, shuffle=True)  
	Valid = torch.utils.data.DataLoader(valid_inputs, batch_size=args.batch_size, shuffle=False)
	Test = torch.utils.data.DataLoader(test_inputs, batch_size=args.batch_size, shuffle=False)
	Allgene = torch.utils.data.DataLoader(all_inputs, batch_size=args.batch_size, shuffle=False)
	return Train,Valid,Test,Allgene

#模型构建

def batch_product(iput,mat2): #输入 xx和bin_context_vector
		result = None
		for i in range(iput.size()[0]):
			op = torch.mm(iput[i], mat2) #torch.mm是两个矩阵相乘，即两个二维的张量相乘
			op = op.unsqueeze(0)  #对数据维度进行扩充。给指定位置加上维数为一的维度
			if(result is None): 
				result = op
			else:
				result = torch.cat((result,op),0)
		return result.squeeze(2)  #a.squeeze(N) 就是去掉a中指定的维数为一的维度

#注意力机制
#nn.Module是Pytorch封装的一个类，是搭建神经网络时需要继承的父类，以下class继承nn.Module的所有特性
class rec_attention(nn.Module): 
	# attention with bin context vector per HM and HM context vector
	def __init__(self,hm,args):
		super(rec_attention,self).__init__() #父类初始化，等价于nn.Module.__init__()
		self.num_directions=2 if args.bidirectional else 1  #如果是双向lstm，num_direction=2,否则为1
		if (hm==False):  #bin级注意力,双向
			self.bin_rep_size=args.bin_rnn_size*self.num_directions #32*2 
		else:    #HM级注意力
			self.bin_rep_size=args.bin_rnn_size   #32 
#生成单精度浮点类型的tensor,转换成可以训练的类型 parameter
		self.bin_context_vector=nn.Parameter(torch.Tensor(self.bin_rep_size,1),requires_grad=True) 
#nn.Softmax(dim=n),不同维度归一化,得到一个与输入具有相同维度和形状的张量，其值在 [0, 1] 范围内
		self.softmax=nn.Softmax(dim=1)  

		self.bin_context_vector.data.uniform_(-0.1, 0.1) #就是个初始化，uniform_指的是均匀分布

	def forward(self,iput):
		alpha=self.softmax(batch_product(iput,self.bin_context_vector))
		[batch_size,source_length,bin_rep_size2]=iput.size()
		repres=torch.bmm(alpha.unsqueeze(2).view(batch_size,-1,source_length),iput)
		return repres,alpha


class recurrent_encoder(nn.Module): 
	# modular LSTM encoder
	def __init__(self,n_bins,ip_bin_size,hm,args): #定义初始化方法
		super(recurrent_encoder,self).__init__() 
		self.bin_rnn_size=args.bin_rnn_size #32
		self.ipsize=ip_bin_size #1
		self.seq_length=n_bins #200

		self.num_directions=2 if args.bidirectional else 1
		if (hm==False):
			self.bin_rnn_size=args.bin_rnn_size
		else:
			self.bin_rnn_size=args.bin_rnn_size // 2
		self.bin_rep_size=self.bin_rnn_size*self.num_directions


		self.rnn=nn.LSTM(self.ipsize,self.bin_rnn_size,num_layers=args.num_layers,dropout=args.dropout,bidirectional=args.bidirectional)

		self.bin_attention=rec_attention(hm,args)
	def outputlength(self):
		return self.bin_rep_size
	def forward(self,single_hm,hidden=None):

		bin_output, hidden = self.rnn(single_hm,hidden)
		bin_output = bin_output.permute(1,0,2)
		hm_rep,bin_alpha = self.bin_attention(bin_output)
		return hm_rep,bin_alpha



class raw_d(nn.Module):
	def __init__(self,args):
		self.n_features=args.n_features    #2
		self.n_bins=args.n_bins  #160
		self.ip_bin_size=1
		super(raw_d,self).__init__() #nn.Module初始化
		self.rnn_hms=nn.ModuleList() #torch.nn.ModuleList：可以把任意 nn.Module 的子类加到这个 list 里面
		for i in range(self.n_features):  #[0, 1, 2, 3, 4]
			self.rnn_hms.append(recurrent_encoder(self.n_bins,self.ip_bin_size,False,args)) #在ModuleList后面添加网络层
		self.opsize = self.rnn_hms[0].outputlength() #bin_rep_size=self.bin_rnn_size*self.num_directions
		self.hm_level_rnn_1=recurrent_encoder(self.n_features,self.opsize,True,args)
		self.opsize2=self.hm_level_rnn_1.outputlength()
		self.diffopsize=2*(self.opsize2)
		self.fdiff1_1=nn.Linear(self.opsize2,1)

	def forward(self,iput): ##输入组蛋白数据

		iput=iput
		bin_a=None
		level1_rep=None
		[batch_size,_,_]=iput.size()

		for hm,hm_encdr in enumerate(self.rnn_hms):
			hmod=iput[:,:,hm].contiguous() #取3维中的第一维，contiguous()使tensor变量在内存中的存储变得连续
			hmod=torch.t(hmod).unsqueeze(2) #转置，升维度

			op,a= hm_encdr(hmod)
			if level1_rep is None:
				level1_rep=op
				bin_a=a
			else:
				level1_rep=torch.cat((level1_rep,op),1) #(10,5,64)
				bin_a=torch.cat((bin_a,a),1) 
		level1_rep=level1_rep.permute(1,0,2) #块，行，列 0，1，2  转置 (5,10,64)
		final_rep_1,hm_level_attention_1=self.hm_level_rnn_1(level1_rep)
		final_rep_1=final_rep_1.squeeze(1)
		prediction_m=((self.fdiff1_1(final_rep_1)))
		return prediction_m,hm_level_attention_1, bin_a

#train model
model_dir = args.save_root
attentionmapfile=model_dir+'/'+args.attentionfilename

model = raw_d(args)
print(model)
if(args.test_on_saved_model==False):
	print("==>initializing a new model") #测试不是在已有模型上进行，初始化新模型
	for p in model.parameters():   
		p.data.uniform_(-0.1,0.1) 
        
if torch.cuda.device_count() >= 1:  #可用GPU数量
	torch.cuda.manual_seed_all(1)   #为所有的GPU设置种子
	dtype = torch.cuda.FloatTensor  #数据类型为“GPU的单精度浮点数张量”
	cuda.set_device(args.gpuid)     #训练时指定GPU显卡
	model.type(dtype)               
	print('Using GPU '+str(args.gpuid))
else:
	dtype = torch.FloatTensor  

optimizer = optim.Adam(model.parameters(), lr = args.lr)
Loss = nn.MSELoss(size_average=True).type(dtype)

#####Train
def train(Train):
    model.train() 
    exp_targets = torch.zeros(Train.dataset.__len__(),1)  #创建1560行1列全为0的tensor
    exp_predictions = torch.zeros(exp_targets.size(0),1)  #1560
    all_attention_bin=torch.zeros(Train.dataset.__len__(),(args.n_features*args.n_bins))  #1560行320列0的tensor
    all_attention_hm=torch.zeros(Train.dataset.__len__(),args.n_features)  
    num_batches = int(math.ceil(Train.dataset.__len__()/float(args.batch_size))) #取整数20(一共取20次数据去训练)
    gene_ids=[None]*Train.dataset.__len__()   #一列200个None,基因id               
    per_epoch_loss = 0 
    for idx, Sample in enumerate(Train): #开始迭代,10个基因为1个sample的迭代输入 0<=idx<=19
        if(idx%100==0):    #除100是否有余数
            print('TRAINING ON BATCH:',idx)
        start,end = (idx*args.batch_size), min((idx*args.batch_size)+args.batch_size, Train.dataset.__len__())
        optimizer.zero_grad()
        # get HM profiles
        iput = Sample['Features'] 
        # get targets: both differential and cell specific expression
        batch_exp_targets=(Sample['expr']).float().unsqueeze(1)
        exp_targets[start:end,0] = batch_exp_targets[:,0]
        
        gene_ids[start:end]=Sample['geneID']
        batch_size = iput.size(0)
#注意力权重
        batch_exp_predictions,batch_beta,batch_alpha = model(iput.type(dtype)) 
        all_attention_bin[start:end]=batch_alpha.data #α是bin级注意力权重
        all_attention_hm[start:end]=batch_beta.data   #β是HM级注意力权重
        loss = Loss(batch_exp_predictions,batch_exp_targets.type(dtype))

        exp_predictions[start:end] = batch_exp_predictions.data.cpu()
        per_epoch_loss += loss.item()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), args.clip)
        optimizer.step()
    per_epoch_loss=per_epoch_loss/num_batches
    return exp_predictions,exp_targets,all_attention_bin,all_attention_hm,per_epoch_loss,gene_ids



#test
def test(MyValid):
    model.eval()  #验证模式，保证BN层使用全部训练数据的均值和方差，利用所有的网络连接来训练更新参数。

    exp_targets = torch.zeros(MyValid.dataset.__len__(),1)  #创建660行1列全为0的tensor
    exp_predictions = torch.zeros(exp_targets.size(0),1)   #200行1列
    all_attention_bin=torch.zeros(MyValid.dataset.__len__(),(args.n_features*args.n_bins))
    all_attention_hm=torch.zeros(MyValid.dataset.__len__(),args.n_features)

    num_batches = int(math.ceil(MyValid.dataset.__len__()/float(args.batch_size)))
    all_gene_ids=[None]*MyValid.dataset.__len__()
    per_epoch_loss = 0
    for idx, Sample in enumerate(MyValid):
        if(idx%668==0):
            print('TESTING ON BATCH:',idx)
        start,end = (idx*args.batch_size), min((idx*args.batch_size)+args.batch_size, MyValid.dataset.__len__())
        optimizer.zero_grad()
		# get HM profiles
        iput = Sample['Features']
		# get targets: both differential and cell specific expression
        batch_exp_targets=(Sample['expr']).float().unsqueeze(1)
        exp_targets[start:end,0] = batch_exp_targets[:,0]
        
        all_gene_ids[start:end]=Sample['geneID']
        batch_size = iput.size(0)

        batch_exp_predictions,batch_beta,batch_alpha = model(iput.type(dtype))  #α是bin级注意力权重，β是HM级注意力权重
        all_attention_bin[start:end]=batch_alpha.data
        all_attention_hm[start:end]=batch_beta.data
        loss = Loss(batch_exp_predictions,batch_exp_targets.type(dtype))

        exp_predictions[start:end] = batch_exp_predictions.data.cpu()
        per_epoch_loss += loss.item()
    per_epoch_loss=per_epoch_loss/num_batches
    return exp_predictions,exp_targets,all_attention_bin,all_attention_hm,per_epoch_loss,all_gene_ids

def compute_metrics(predictions, targets):

	pred=predictions.numpy()
	targets=targets.numpy()

	R2,p=scipy.stats.pearsonr(np.squeeze(targets),np.squeeze(pred))
	MSE=metrics.mean_squared_error(targets, pred)
	return MSE, R2
#数据导入
print('==>processing data')
Train,Valid,Test,Allgene = load_data(args) 

#######################result
best_valid_loss = 10000000000
best_valid_MSE=100000
best_valid_R2=-1
if(args.test_on_saved_model==False):
	for epoch in range(0, args.epochs):
		print('=---------------------------------------- Training '+str(epoch+1)+' -----------------------------------=')
		exp_predictions,exp_targets,alpha_train,beta_train,train_loss,train_geneID = train(Train)
		train_MSE, train_R2 = compute_metrics(exp_predictions,exp_targets)
		np.savetxt(args.save_root+"/"+"train_predictions.txt",np.array(exp_predictions))
		np.savetxt(args.save_root+"/"+"train_targets.txt",np.array(exp_targets))
		np.savetxt(args.save_root+"/"+"alpha_train.txt",np.array(alpha_train))
		np.savetxt(args.save_root+"/"+"beta_train.txt",np.array(beta_train))
		np.savetxt(args.save_root+"/"+"train_geneID.txt",np.array(train_geneID),fmt='%s')
		exp_predictions,exp_targets,alpha_valid,beta_valid,valid_loss,gene_ids_valid = test(Valid)
		valid_MSE, valid_R2 = compute_metrics(exp_predictions,exp_targets)
		np.savetxt(args.save_root+"/"+"valid_predictions.txt",np.array(exp_predictions))
		np.savetxt(args.save_root+"/"+"valid_targets.txt",np.array(exp_targets))
		np.savetxt(args.save_root+"/"+"valid_geneID.txt",np.array(gene_ids_valid),fmt='%s')
		np.savetxt(args.save_root+"/"+"alpha_valid.txt",np.array(alpha_valid))
		np.savetxt(args.save_root+"/"+"beta_valid.txt",np.array(beta_valid))
		if(valid_R2 >= best_valid_R2):
			 #save best epoch -- models converge early
			best_valid_R2=valid_R2
			torch.save(model,model_dir+"/"+model_name+'_R2_model.pt')

		print("train_MSE:",train_MSE)
		print("Epoch:",epoch)
		print("train_MSE:",train_MSE)
		print("train R2:",train_R2)
		print("valid R2:",valid_R2)
		print("best valid R2:", best_valid_R2)
		print("train_loss:",train_loss)
		print("valid_loss:",valid_loss)

	print("finished training!!")
	print("best validation R2:",best_valid_R2)
	print("testing")

	model=torch.load(model_dir+"/"+model_name+'_R2_model.pt')

	exp_predictions,exp_targets,alpha_test,beta_test,test_loss,gene_ids_test = test(Test)
	test_MSE, test_R2 = compute_metrics(exp_predictions,exp_targets)
	np.savetxt(args.save_root+"/"+"test_predictions.txt",np.array(exp_predictions))
	np.savetxt(args.save_root+"/"+"test_targets.txt",np.array(exp_targets))
	np.savetxt(args.save_root+"/"+"alpha_test.txt",np.array(alpha_test))
	np.savetxt(args.save_root+"/"+"beta_test.txt",np.array(beta_test))
	np.savetxt(args.save_root+"/"+"test_geneID.txt",np.array(gene_ids_test),fmt='%s')
	all_predictions,all_targets,alpha_all,beta_all,all_loss,gene_ids_all = test(Allgene)
	all_MSE, all_R2 = compute_metrics(all_predictions,all_targets)
	np.savetxt(args.save_root+"/"+"all_predictions.txt",np.array(all_predictions))
	np.savetxt(args.save_root+"/"+"all_targets.txt",np.array(all_targets))
	np.savetxt(args.save_root+"/"+"all_gene_id.txt",np.array(gene_ids_all),fmt='%s')
	np.savetxt(args.save_root+"/"+"alpha_all.txt",np.array(alpha_all))
	np.savetxt(args.save_root+"/"+"beta_all.txt",np.array(beta_all))
	
	print("test R2:",test_R2)
	print("test_MSE:",test_MSE)
	print("all_MSE:",all_MSE)
	print("all R2:",all_R2)
#	if(args.save_attention_maps==False):
#		attentionfile=open(attentionmapfile,'w')
#		attentionfilewriter=csv.writer(attentionfile)
#		beta_test=beta_test.numpy()
#		for i in range(len(gene_ids_test)):
#			gene_attention=[]
#			gene_attention.append(gene_ids_test[i])
#			for e in beta_test[i,:]:
#				gene_attention.append(str(e))
#			attentionfilewriter.writerow(gene_attention)
#		attentionfile.close()


else:
	model=torch.load(model_dir+"/"+model_name+'_R2_model.pt')
	exp_predictions,exp_targets,alpha_test,beta_test,test_loss,gene_ids_test = test(Test)
	test_MSE, test_R2 = compute_metrics(exp_predictions,exp_targets)
	np.savetxt(args.save_root+"/"+"test_predictions.txt",np.array(exp_predictions))
	np.savetxt(args.save_root+"/"+"test_targets.txt",np.array(exp_targets))
	np.savetxt(args.save_root+"/"+"alpha_test.txt",np.array(alpha_test))
	np.savetxt(args.save_root+"/"+"beta_test.txt",np.array(beta_test))
	np.savetxt(args.save_root+"/"+"test_geneID.txt",np.array(gene_ids_test),fmt='%s')
	all_predictions,all_targets,alpha_all,beta_all,all_loss,gene_ids_all = test(Allgene)
	all_MSE, all_R2 = compute_metrics(all_predictions,all_targets)
	np.savetxt(args.save_root+"/"+"all_predictions.txt",np.array(all_predictions))
	np.savetxt(args.save_root+"/"+"all_targets.txt",np.array(all_targets))
	np.savetxt(args.save_root+"/"+"all_gene_id.txt",np.array(gene_ids_all),fmt='%s')
	np.savetxt(args.save_root+"/"+"alpha_all.txt",np.array(alpha_all))
	np.savetxt(args.save_root+"/"+"beta_all.txt",np.array(beta_all))
	
	print("test R2:",test_R2)
	print("test_MSE:",test_MSE)
	print("all_MSE:",all_MSE)
	print("all R2:",all_R2)

#	if(args.save_attention_maps):
#		attentionfile=open(attentionmapfile,'w')
#		attentionfilewriter=csv.writer(attentionfile)
#		beta_test=beta_test.numpy()
#		for i in range(len(gene_ids_test)):
#			gene_attention=[]
#			gene_attention.append(gene_ids_test[i])
#			for e in beta_test[i,:]:
#				gene_attention.append(str(e))
#			attentionfilewriter.writerow(gene_attention)
#		attentionfile.close()



