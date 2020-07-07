%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: solve regularized LSP with sequential thresh.
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz


function [W,G,b,resid,dW] = wsindy_pde_RGLS(lambda,gamma,Theta_pdx,lhs_ind,axi)

b = zeros(size(Theta_pdx,1),length(lhs_ind));
W = zeros(size(Theta_pdx,2)-1,length(lhs_ind));
resid = zeros(size(Theta_pdx,1),length(lhs_ind));
dW = cell(length(lhs_ind)+1,1);
G = cell(length(lhs_ind)+1,1);

for k=1:length(lhs_ind)
    b(:,k) = Theta_pdx(:,lhs_ind(k));
    G{k} = Theta_pdx(:,[1:lhs_ind(k)-1 lhs_ind(k)+1:end]);
    W(:,k) = sparsifyDynamics(G{k}, b(:,k), lambda, 1, gamma);
    resid(:,k) = b(:,k) - G{k} * W(:,k);
    dW{k+1} = W(:,k)-axi(:,k);
    dW{1}(k) = norm(dW{k+1},inf);
end

end
