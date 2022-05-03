clear; clc;
matdir = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';
load([matdir filesep 'eNDM_mousedata.mat'],'Networks');
load([matdir filesep 'Zeisel_CTMaps.mat'],'Zeisel_ng1360');
load([matdir filesep 'classkey_zeisel.mat'],'classkey_zeisel');
load([matdir filesep 'mouse_aSynData_426.mat']);
C = Networks.ret/max(Networks.ret(:));
%% Human, no cell types
% Nexis:global, no directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
U_null = zeros(426,1);
alpha = 0.109;
beta = 5.085;
gamma = 0.423;
b = 0;
p = 0;
s = 0.5;
wdir = 0;
color = [0 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U_null,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{},'cell',0,wdir,color,'s',1);

% Nexis:global, with directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
U = zeros(426,1);
alpha = 0.1;
beta = 4.17;
gamma = 0.4;
b = 0;
p = 0;
s = 0.72;
wdir = 1;
color = [0 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{},'cell',0,wdir,color,'s',1);
%% Human HBNOR
% Nexis:HBNOR, no directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'HBNOR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.04;
beta = 5.65;
gamma = 0.36;
b = -1.05;
p = 0.37;
s = 0.5;
wdir = 0;
color = [1 0 0];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

% Nexis:HBNOR, with directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'HBNOR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.05;
beta = 4.65;
gamma = 0.36;
b = -1.05;
p = 0.37;
s = 0.73;
wdir = 1;
color = [1 0 0];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

%% Human MBDOP1
% Nexis:MBDOP1, no directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'MBDOP1';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.1;
beta = 4.91;
gamma = 0.39;
b = -0.82;
p = 0.45;
s = 0.5;
wdir = 0;
color = [1 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

% Nexis:MBDOP1, with directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'MBDOP1';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.08;
beta = 4.02;
gamma = 0.38;
b = -0.29;
p = 0.8;
s = 0.72;
wdir = 1;
color = [1 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

%% Human MBDOP2
% Nexis:MBDOP2, no directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'MBDOP2';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.13;
beta = 5.6;
gamma = 0.4;
b = 13.88;
p = -4.72;
s = 0.5;
wdir = 0;
color = [0.5 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

% Nexis:MBDOP2, with directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'MBDOP2';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.13;
beta = 4.32;
gamma = 0.42;
b = 44.9;
p = -11.98;
s = 0.81;
wdir = 1;
color = [0.5 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

%% Human HBADR
% Nexis:HBADR, no directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'HBADR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.101;
beta = 5.329;
gamma = 0.39;
b = -0.845;
p = 0.014;
s = 0.5;
wdir = 0;
color = [0 0.5 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);

% Nexis:HBADR, with directionality, human
ts = tpts.human;
data = data426.human;
data = data/nansum(data(:,1));
seed = seed426.human;
type = 'HBADR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.098;
beta = 4.358;
gamma = 0.37;
b = -0.429;
p = 0.051;
s = 0.709;
wdir = 1;
color = [0 0.5 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'human',{type},'cell',0,wdir,color,'s',1);