clear; clc;
matdir = '/Users/justintorok/Documents/MATLAB/Nexis/raw_data_mouse';
load([matdir filesep 'eNDM_mousedata.mat'],'Networks');
load([matdir filesep 'Zeisel_CTMaps.mat'],'Zeisel_ng1360');
load([matdir filesep 'classkey_zeisel.mat'],'classkey_zeisel');
load([matdir filesep 'mouse_aSynData_426.mat']);
C = Networks.ret/max(Networks.ret(:));
%% Mouse, no cell types
% Nexis:global, no directionality, mouse
ts = tpts.mouse;
data = data426.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
U_null = zeros(426,1);
alpha = 0.19;
beta = 6.87;
gamma = 0.43;
b = 0;
p = 0;
s = 0.5;
wdir = 0;
color = [0 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U_null,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{},'cell',0,wdir,color,'o',1);

% Nexis:global, with directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
U = zeros(426,1);
alpha = 0.18;
beta = 5.59;
gamma = 0.41;
b = 0;
p = 0;
s = 0.74;
wdir = 1;
color = [0 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{},'cell',0,wdir,color,'o',1);
%% Mouse HBNOR
% Nexis:HBNOR, no directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'HBNOR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.14;
beta = 7.39;
gamma = 0.36;
b = -1.03;
p = 0.29;
s = 0.5;
wdir = 0;
color = [1 0 0];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

% Nexis:HBNOR, with directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'HBNOR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.15;
beta = 6.12;
gamma = 0.37;
b = -0.75;
p = 0.19;
s = 0.75;
wdir = 1;
color = [1 0 0];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

%% Mouse MBDOP1
% Nexis:MBDOP1, no directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'MBDOP1';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.14;
beta = 5.98;
gamma = 0.39;
b = 1.22;
p = 1.89;
s = 0.5;
wdir = 0;
color = [1 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

% Nexis:MBDOP1, with directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'MBDOP1';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.15;
beta = 5.51;
gamma = 0.39;
b = 0.21;
p = 1.21;
s = 0.73;
wdir = 1;
color = [1 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

%% Mouse MBDOP2
% Nexis:MBDOP2, no directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'MBDOP2';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.19;
beta = 7.39;
gamma = 0.42;
b = 23.94;
p = 0.53;
s = 0.5;
wdir = 0;
color = [0.5 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

% Nexis:MBDOP2, with directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'MBDOP2';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.2;
beta = 6.06;
gamma = 0.44;
b = 38.04;
p = -7.05;
s = 0.8;
wdir = 1;
color = [0.5 0 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

%% Mouse HBADR
% Nexis:HBADR, no directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'HBADR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.154;
beta = 6.947;
gamma = 0.404;
b = -0.375;
p = 0.349;
s = 0.5;
wdir = 0;
color = [0 0.5 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);

% Nexis:HBADR, with directionality, mouse
ts = tpts.mouse;
data = data/nansum(data(:,1));
seed = seed426.mouse;
type = 'HBADR';
cellind = ismember(classkey_zeisel,type);
U = Zeisel_ng1360(:,cellind);
U = (U - min(U))/(max(U) - min(U));
alpha = 0.148;
beta = 5.709;
gamma = 0.382;
b = 0.014;
p = 0.388;
s = 0.727;
wdir = 1;
color = [0 0.5 1];
[y] = eNDM_general_dir(seed*gamma,ts,C,U,alpha,beta,s,0,b,p,'analytic',1);
CorrelationPlotter_single(y,data,ts,'mouse',{type},'cell',0,wdir,color,'o',1);