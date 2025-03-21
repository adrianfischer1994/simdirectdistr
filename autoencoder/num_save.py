import torch
import math

from autoencoder import ModelAE

modelS = []

for label in range(0,10):
    modelS.append(torch.load('result\\modelS' + str(label)))
    modelS[label].eval()


file_read = open("gen_sph.txt", "r")
file_write = open("num.txt", "w")

firstline = True
k = 0

for line in file_read.readlines():
    
    x = line.split()
    for i in range(0,len(x)):
        x[i] = float(x[i])

    tt = modelS[math.floor(k/10)].decode(torch.tensor(x)).reshape(28, 28).tolist()

    if firstline == False:
        _=file_write.write("\n ")
    firstline = False
    
    for i in range(0, 28):
        for j in range(0, 28):
            _=file_write.write("%s " % tt[i][j])
        if i != 27:
            _=file_write.write("\n ")
            
    k = k + 1

file_read.close()
file_write.close()
