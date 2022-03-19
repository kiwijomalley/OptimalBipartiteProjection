# To binarize the network based on certain threshold
# Originally written by Chuankai An
# Adapted by Xin Ran on 10/22/2019
# edited by Xin Ran on 10/30/2019, turns out reading in the dataset using pandas will be faster\
# edited by James O'Malley 08/18/2021 to correct filepath
# 08/29/2021 - James updated filepaths for Sound project

#Need to run separately for undirected and directed networks?
#Was this done previously for the CVD factorial-design networks

#import random
import pandas as pd
import networkx as nx
import numpy as np
import math
import decimal
# changed on 11/24/2019
#path = "/projects/active/14593/idata/core_c/core_c_jomalley/"
path = "/projects/active/54054/idata/jomalley/"
#path = "/projects/active/14593/idata/Project3_HospClosure_VQI/xran/CAn_network/"

# fp=open(path + "undirected_network_2010_new.txt",'r')

df = pd.read_csv(path + "undirected_network_2018a.txt", sep = ";", header = None, index_col = False)
#df = pd.read_csv(path + "directed_network_2018a.txt", sep = ";", header = None, index_col = False)
print('hey is loading the dataset')\
df = df.iloc[:, :-1] # drop the last column due to the deliminator used
wc = df.shape[0] # number of rows\


# generate the threshold
# why only pick 1 percent of the data uniformly
#size=100000
#size = 10
ratio = 0.2 # remove lowest 20% weights
prob = 0.01
size = int(wc * prob)
sampleID = np.random.choice(range(wc), size = size, replace = False)
sample = df.iloc[sampleID, :]

# get the threshold for each weight, sorted - ascending by default
threshold = list(sample.iloc[:, 1:].apply(lambda x: sorted(x)[int(size * ratio)], axis = 0))
print(threshold)

# [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 0.25, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]



# sample=[[0 for col in range(size)] for row in range(20)]
# sample = pd.DataFrame(index = range(20),columns = range(size))
# index = 0

# regenerate the threshold by using the cutoff at the 20% strength
# fixed sample probability

# number of rows, w/o replacement, 100000 out of it
# generate a random number for them (randomly order the rows with equal probability), then sorted and take the first 3
# wi/0 predifiend sample size, just use random number generator 

# total number of rows, prob = sample size/ toral number of rows, pick each row with this prob

# c1 = 0
# while 1:
#     line = fp.readline()
#     c1 += 1
#     if not line:
#         break
#     if c1 % 10000 == 0:
#         print(c1)
#     # if c1 == 100:
#     #     break
#     values = line.split(";")[1:-1] # there's /n in the end
#     tmp = random.randint(1, 100)
#     if tmp == 1: 
#         sample.iloc[:, index] = list(map(lambda x: float(x), values)) 
#         index = index + 1
#         if index == size:
#             break
# fp.close()
# print(sample.iloc[:, 0:100])


# threshold=[0 for col in range(20)]
# for i in range(20):
#     tmp = sorted(sample.iloc[i, :])
#     threshold[i] = tmp[int(size*ratio)] # the cutoff of weight below 20%, for each weight type
# print(threshold)


# calculate the network density under each weight
# col = ["edge"]
# col.extend(["w" + str(i) for i in range(20)])
# df = pd.read_csv(path + "undirected_network_2010_new.txt", sep = ";", header = None, index_col = False)
# print('hey is loading the dataset')
# df = df.iloc[:, :-1] # drop the last column due to the deliminator used
# df.columns = col # assign columns names

non_zero_edge = [sum(df.iloc[:, i] > 0) for i in range(1, len(df.columns))]

nodes = list(map(lambda x: x.split('_'), df.iloc[:, 0]))
nodes = [item for sub in nodes for item in sub] # flattern the list
print("number of total nodes is %d" %len(nodes))
nodes = set(nodes) # unique nodes 
print("number of unique nodes is %d" %len(nodes))

# def LogDen(non_zero_edge, nodes):
#     return math.log(2) + math.log(non_zero_edge) - math.log(len(nodes)) - math.log(len(nodes) - 1)

# dens1 = [math.exp(LogDen(non_zero_edge[i], nodes)) for i in range(20)] # are all zeros, too small # use the log or split them

def Dec(num):
    return decimal.Decimal(num)
# undirected network
dens2 = [round(2 * Dec(non_zero_edge[i]) / (Dec(len(nodes)) * Dec(len(nodes) - 1)), 8) for i in range(20)] # display more decimals
# directed network
#dens2 = [round(Dec(non_zero_edge[i]) / (Dec(len(nodes)) * Dec(len(nodes) - 1)), 8) for i in range(20)] # display more decimals
#print(dens1)
print(dens2)

# densDF = pd.DataFrame(dens1)
# densDF.to_csv(path + "undirected_network_2010_new_density_v1.csv", header = False, index = False)
densDF2 = pd.DataFrame(dens2)
densDF2.to_csv(path + "undirected_network_2018a_density.txt", header = False, index = False)
#densDF2.to_csv(path + "directed_network_2018a_density.txt", header = False, index = False)


# use dens2, the second approach



# apply the threshold
def binarize(col, den, th, ratio, i):
    if i%5 != 0:
        if ratio <= den:
            col = list(map(lambda x: 0 if x <= th else 1, col))
        else:
            col = list(map(lambda x: 1 if x > 0 else 0, col))
    else:
        col = list(map(lambda x: 0 if x == 0 else 1, col))
    return col


for i in range(1, len(df.columns)):\
    print(i)
    df.iloc[:, i] = binarize(df.iloc[:, i], dens2[i - 1], threshold[i - 1], ratio, i)

df.to_csv(path + "binary_undir_network_2018a.txt", header = False, index = False) # 0225, just change it to .txt
#df.to_csv(path + "binary_dir_network_2018a.txt", header = False, index = False) # 0225, just change it to .txt}