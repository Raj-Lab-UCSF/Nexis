function [flow_mat] = FlowCalculator(X_,C_,beta_,isndm,U_,b_)
% This function will calculate the flows between each pair of regions using
% the NDM and eNDM frameworks

%isndm is a binary flag - where replaced by 1 means the model is for the
%standard ndm and if 0 then it means the model is for eNDM

nroi = size(C_,1);
%number of time points
nt = size(X_,2);
%initialize the flow matrix
flow_mat= zeros(nroi,nroi,nt);
% flow_mat_endm = flow_mat_ndm;

for k = 1:nt
    for i = 1:nroi
        for j = 1:nroi
            if ~isnan(X_(i,k)) && ~isnan(X_(j,k))
                if isndm
                    flow_mat(i,j,k) = beta_*C_(i,j)*(X_(i,k) - X_(j,k));
                else
                    flow_mat(i,j,k) = beta_*C_(i,j)*(X_(i,k)*(1 + dot(b_,U_(i,:)))-...
                        X_(j,k)*(1 + dot(b_,U_(j,:))));
                end
            end
        end
    end
end
% for k = 1:nt
%     for i = 1:nroi
%         flow_mat(i,:,k) = beta_*C_(i,:)'.*(X_(i,k) - X_(:,k));
%     end
% end
end