#Parts of this code are taken from https://github.com/nicola-decao/s-vae-pytorch/


import os
import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data
from torchvision import datasets, transforms
from collections import defaultdict

train_dat = datasets.MNIST('./data', train=True, download=True,transform=transforms.ToTensor())
test_dat = datasets.MNIST('./data', train=False, download=True,transform=transforms.ToTensor())

train_loader = []
test_loader = []

for label in range(0, 10):
    labels = [label]
    indices = [idx for idx, target in enumerate(train_dat.targets) if target in labels]
    train_loader.append(torch.utils.data.DataLoader(torch.utils.data.Subset(train_dat, indices),batch_size=64,shuffle=True))
    test_loader.append(torch.utils.data.DataLoader(torch.utils.data.Subset(train_dat, indices),batch_size=64,shuffle=False))

class ModelAE(torch.nn.Module):
    
    def __init__(self, h_dim, z_dim, activation=nn.ReLU()):
        """
        ModelAE initializer
        :param h_dim: dimension of the hidden layers
        :param z_dim: dimension of the latent representation
        :param activation: callable activation function
        """
        super(ModelAE, self).__init__()
        
        self.z_dim, self.activation = z_dim, activation
        
        # 2 hidden layers encoder
        self.fc_e0 = nn.Linear(784, h_dim * 2)
        nn.init.uniform_(self.fc_e0.weight,a=-1/math.sqrt(784),b=1/math.sqrt(784))
        self.fc_e1 = nn.Linear(h_dim * 2, h_dim)
        nn.init.uniform_(self.fc_e1.weight,a=-1/math.sqrt(h_dim * 2),b=1/math.sqrt(h_dim * 2))
        self.fc_sphere = nn.Linear(h_dim, z_dim)
        nn.init.uniform_(self.fc_sphere.weight,a=-1/math.sqrt(h_dim),b=1/math.sqrt(h_dim))

            
        # 2 hidden layers decoder
        self.fc_d0 = nn.Linear(z_dim, h_dim)
        nn.init.uniform_(self.fc_d0.weight,a=-1/math.sqrt(z_dim),b=1/math.sqrt(z_dim))
        self.fc_d1 = nn.Linear(h_dim, h_dim * 2)
        nn.init.uniform_(self.fc_d1.weight,a=-1/math.sqrt(h_dim),b=1/math.sqrt(h_dim))
        self.fc_logits = nn.Linear(h_dim * 2, 784)
        nn.init.uniform_(self.fc_logits.weight,a=-1/math.sqrt(h_dim * 2),b=1/math.sqrt(h_dim * 2))


    def encode(self, x):
        x = self.activation(self.fc_e0(x))
        x = self.activation(self.fc_e1(x))
        z = self.fc_sphere(x)
        z = z - z.mean()
        z = z / z.norm(dim=-1, keepdim=True)
        return z
        
    def decode(self, z):
        x = self.activation(self.fc_d0(z))
        x = self.activation(self.fc_d1(x))
        x = F.sigmoid(self.fc_logits(x))
        return x
        
    def forward(self, x): 
        z = self.encode(x)
        x_ = self.decode(z)       
        return z, x_
    

def train(model, optimizer, train_loader_numb):
    for i, (x_mb, y_mb) in enumerate(train_loader_numb):

        optimizer.zero_grad()
        
        # dynamic binarization
        x_mb = (x_mb > torch.distributions.Uniform(0, 1).sample(x_mb.shape)).float()

        _, x_mb_ = model(x_mb.reshape(-1, 784))

        loss_recon = nn.BCEWithLogitsLoss(reduction='none')(x_mb_, x_mb.reshape(-1, 784)).sum(-1).mean()

        loss = loss_recon

        loss.backward()
        optimizer.step()
            
            
def test(model, optimizer, test_loader_numb):
    recon_loss = []
    for x_mb, y_mb in test_loader_numb:
        
        # dynamic binarization
        x_mb = (x_mb > torch.distributions.Uniform(0, 1).sample(x_mb.shape)).float()
        
        _, x_mb_ = model(x_mb.reshape(-1, 784))
        
        recon_loss.append(float(nn.BCEWithLogitsLoss(reduction='none')(x_mb_,x_mb.reshape(-1, 784)).sum(-1).mean().data))

    return np.mean(-recon_loss)


if __name__ == "__main__":
    
    # file in which results and models are saved
    filename = 'result'
    os.makedirs(filename)
    
    # hidden dimension and dimension of latent space
    H_DIM = 128
    Z_DIM = 3

    for label in range(2, 3):

        # hyper-spherical AE
        modelS = ModelAE(h_dim=H_DIM, z_dim=Z_DIM + 1)
        optimizerS = optim.Adam(modelS.parameters(), lr=1e-3)
        
        # warm up
        warm_up_iter = 1500
        for i in range(0, warm_up_iter):
            train(modelS, optimizerS,train_loader[label])

        # training

        # maximum number of training iterations
        max_epoch = 1000

        # look-ahead
        thres = 50
        
        epochs = 0
        loss = []
        for i in range(0, max_epoch):
            train(modelS, optimizerS, train_loader[label])
            epochs = epochs + 1
            loss.append(test(modelS, optimizerS, test_loader[label]))
            if epochs > thres:
                lower = max(1, epochs - thres)
                if max(loss[lower: epochs]) < loss[lower - 1]:
                    break

        torch.save(modelS,filename + '/modelS' + str(label))
        
        ll_ = []
        res_ = []

        for x_mb, y_mb in train_loader[label]:
            for i, dat in enumerate(x_mb):
                dat = dat[0]
                if y_mb[i] == label:
                    ll_.append(dat)
                    z_ ,_ = modelS(dat.reshape(-1, 784))
                    res_.append(z_.tolist()[0])


        with open(filename + '/results' + str(label) + '.txt', 'w') as fp:
            for item in res_:
                for i in range(0,Z_DIM+1):
                    _=fp.write("%s " % item[i])
                _=fp.write("\n ")

    with open(filename + '/param.txt', 'w') as fp:
        _=fp.write("H-DIM: " + str(H_DIM) + "\nZ-DIM: " + str(Z_DIM) + "\nWarm-up: " + str(warm_up_iter) + "\nThreshold for iterations: " + str(thres))           



            
