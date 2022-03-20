% C = Networks.ret;
C = Networks.ret/max(Networks.ret(:));
[y2,a2] = eNDM_general_dir((seed426.mouse * 0.423113229089271),[1,3,6,12],...
    C,zeros(426,1),0.109202364102107,5.08516589545512,0.5,0,0,0,'analytic',1);
