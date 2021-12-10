import numpy as np

class Nexis:
    def __init__(self,C,U,seed_vec,t_vec,volcorrect=0,w_dir=0):
        self.C = C # Connectivity matrix, nROI x nROI
        self.U = U # Matrix or vector of cell type or gene expression, nROI x nTypes
        self.seed_vec = seed_vec # Binary vector indicated seed location, nROI x 1 
        self.t_vec = t_vec # Vector of time points to output model predictions, 1 x nt
        self.volcorrect = volcorrect # Binary flag indicating whether to use volume correction
        self.w_dir = w_dir # Binary flag indicating whether to use directionality or not 
    
    def simulate(self, parameters):
        # Define parameters
        ntypes = np.size(self.U,1)
        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        s = parameters[3]
        b = parameters[4:(ntypes+4)]
        p = parameters[(ntypes+4):]

        # Define diagonal matrix Gamma containing spread-independent terms
        



        def network_transfer(params,w):
            tau_e = params[0]
            tau_i = params[1]
            gii = params[2]
            gei = params[3]
            # pw_scale = params[4]
            pw_scale = 1
            gee = 1


            # Cortical model
            Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
            Fi = np.divide(1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)

            Hed = (1 + (Fe * Fi * gei)/(tau_e * (1j * w + Fi * gii/tau_i)))/(1j * w + Fe * gee/tau_e + (Fe * Fi * gei)**2/(tau_e * tau_i * (1j * w + Fi * gii / tau_i)))
            
            Hid = (1 - (Fe * Fi * gei)/(tau_i * (1j * w + Fe * gee/tau_e)))/(1j * w + Fi * gii/tau_i + (Fe * Fi * gei)**2/(tau_e * tau_i * (1j * w + Fe * gee / tau_e)))

            Htotal = Hed + Hid
            
            return pw_scale*Htotal
        

        freq_mdl = []
        for freq in self._fvec:
            _w = 2 * np.pi * freq
            freq_model = network_transfer(parameters, _w)
            freq_mdl.append(freq_model)

        freq_mdl = np.transpose(np.asarray(freq_mdl))
        
        # freq_out = functions.mag2db(np.abs(freq_mdl))
        
        return np.abs(freq_mdl)
