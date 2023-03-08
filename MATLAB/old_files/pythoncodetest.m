matdir = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';
load([matdir filesep 'eNDM_mousedata.mat'],'Networks');
load([matdir filesep 'KaufmanDiamond_datasets_dat&seed.mat']);
load([matdir filesep 'regional_Gene_Data.mat'],'regvgene_mean');
C = Networks.ret/max(Networks.ret(:));
tpts = tpts.DS6;
data = data426.DS6;
seed = seed426.DS6;
U_null = zeros(426,1);
alpha_glob = 0.21;
beta_glob = 2.09;
gamma_glob = 0.1;
s = 0.5;
[y_glob] = eNDM_general_dir((seed * gamma_glob),tpts,...
    C,U_null,alpha_glob,beta_glob,s,0,0,0,'analytic',0);
U_trem2 = regvgene_mean(:,3578); U_trem2 = U_trem2/mean(U_trem2);
alpha_trem2 = 0.2688;
beta_trem2 = 1.575;
gamma_trem2 = 0.0854;
p_trem2 = -0.3908;
b_trem2 = 0.1749;
[y_trem2] = eNDM_general_dir((seed * gamma_trem2),tpts,...
    C,U_trem2,alpha_trem2,beta_trem2,s,0,b_trem2,p_trem2,'analytic',0);
R_glob = (corr(y_glob(:),data(:),'rows','complete'))^2
R_trem2 = (corr(y_trem2(:),data(:),'rows','complete'))^2
