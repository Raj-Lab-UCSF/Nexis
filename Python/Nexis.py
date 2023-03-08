import os
import numpy as np
import scipy as sp
import scipy.io
from scipy.linalg import expm
import pandas as pd

class loading_matfiles:
    def __init__(self,datadir_=''):
        if (datadir_==''):
            curdir = os.getcwd()
            subdir = 'raw_data_mouse'
            datadir_ = os.path.join(curdir,subdir)
        self.datadir = datadir_ # Directory to load dependences from
    
    def load_file(self,datafile):
        datafile = os.path.join(self.datadir,datafile)
        datamat = scipy.io.loadmat(datafile)
        return datamat

class run_Nexis:
    def __init__(self,C_,U_,seed_vec_,t_vec_,volcorrect_=0,w_dir_=0,datadir_=''):
        self.C = C_ # Connectivity matrix, nROI x nROI
        self.U = U_ # Matrix or vector of cell type or gene expression, nROI x nTypes
        self.seed_vec = seed_vec_ # Binary vector indicated seed location, nROI x 1 
        self.t_vec = t_vec_ # Vector of time points to output model predictions, 1 x nt
        self.volcorrect = volcorrect_ # Binary flag indicating whether to use volume correction
        self.w_dir = w_dir_ # Binary flag indicating whether to use directionality or not 
        if (datadir_==''):
            curdir = os.getcwd()
            subdir = 'raw_data_mouse'
            datadir_ = os.path.join(curdir,subdir)
        self.datadir = datadir_ # Directory to load dependences from

    def forward_sim(self,A_,t_,x0_):
        y_ = np.zeros([np.shape(A_)[0],len(t_)])
        for i in list(range(len(t_))):
            ti = t_[i]
            y_[:,i] = np.dot(expm(A_*ti),np.squeeze(x0_))
        return y_

    def simulate_nexis(self, parameters):
        """
        Returns a matrix, Y, that is nROI x nt representing the modeled Nexis pathology
        given the provided parameters. alpha, beta, and gamma should be nonnegative scalars;
        s should be bounded between 0 and 1; b and p should be nCT-long vectors
        """
        # Define parameters
        ntypes = np.size(self.U,axis=1)
        alpha = parameters[0] # global connectome-independent growth
        beta = parameters[1] # global diffusivity rate 
        gamma = parameters[2] # seed rescale value
        s = parameters[3] # directionality (0 = anterograde, 1 = retrograde)
        if self.w_dir==0:
            s = 0.5
        else:
            s = parameters[3] # directionality (0 = anterograde, 1 = retrograde)
        b = np.transpose(parameters[4:(ntypes+4)]) # cell-type-dependent spread modifier
        p = np.transpose(parameters[(ntypes+4):]) # cell-type-dependent growth modifier
        
        # Define starting pathology x0
        x0 = gamma * self.seed_vec
        
        # Define diagonal matrix Gamma containing spread-independent terms
        s_p = np.dot(self.U,p)
        Gamma = np.diag(s_p) + (alpha * np.eye(len(s_p)))

        # Define Laplacian matrix L
        C_dir = (1-s) * np.transpose(self.C) + s * self.C
        coldegree = np.sum(C_dir,axis=0)
        L_raw = np.diag(coldegree) - C_dir
        s_b = np.dot(self.U,b)
        s_b = np.reshape(s_b,[len(s_b),1])
        S_b = np.tile(s_b,len(s_b)) + np.ones([len(s_b),len(s_b)])
        L = np.multiply(L_raw,np.transpose(S_b))
        
        # Apply volume correction if applicable
        if self.volcorrect:
            regionfile = os.path.join(self.datadir,'regionvoxels.mat')
            volmat = scipy.io.loadmat(regionfile)
            voxels = volmat['voxels']
            voxels_2hem = np.vstack((voxels,voxels)) / 2
            inv_voxels_2hem = np.diag(np.squeeze(voxels_2hem)**(-1))
            L = np.mean(voxels_2hem) * np.dot(inv_voxels_2hem,L)

        # define system dydt = Ax
        A = Gamma - (beta * L)

        # solve analytically
        y = self.forward_sim(A,self.t_vec,x0)
        return y


        # def network_transfer(params,w):
        #     tau_e = params[0]
        #     tau_i = params[1]
        #     gii = params[2]
        #     gei = params[3]
        #     # pw_scale = params[4]
        #     pw_scale = 1
        #     gee = 1


        #     # Cortical model
        #     Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
        #     Fi = np.divide(1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)

        #     Hed = (1 + (Fe * Fi * gei)/(tau_e * (1j * w + Fi * gii/tau_i)))/(1j * w + Fe * gee/tau_e + (Fe * Fi * gei)**2/(tau_e * tau_i * (1j * w + Fi * gii / tau_i)))
            
        #     Hid = (1 - (Fe * Fi * gei)/(tau_i * (1j * w + Fe * gee/tau_e)))/(1j * w + Fi * gii/tau_i + (Fe * Fi * gei)**2/(tau_e * tau_i * (1j * w + Fe * gee / tau_e)))

        #     Htotal = Hed + Hid
            
        #     return pw_scale*Htotal
        

        # freq_mdl = []
        # for freq in self._fvec:
        #     _w = 2 * np.pi * freq
        #     freq_model = network_transfer(parameters, _w)
        #     freq_mdl.append(freq_model)

        # freq_mdl = np.transpose(np.asarray(freq_mdl))
        
        # # freq_out = functions.mag2db(np.abs(freq_mdl))
        
        # return np.abs(freq_mdl)
