function [W,G,b,resid,dW,its_all,thrs_EL,M,lambda_hat,lossvals,ET_wsindy,tags_pde_G,lib_list_G] = wsindy_pde_solve(lambda,gamma,Theta_pdx,lhs_ind,true_nz_weights,M_full,maxits,tags_pde,lib_list, sparsity_scale)
if and(length(lambda)==1,all(lambda>=0))
    [W,G,b,resid,dW,its_all,thrs_EL,M] = wsindy_pde_RGLS(lambda,gamma,Theta_pdx,lhs_ind,true_nz_weights,M_full,maxits);
    lambda_hat = lambda;
    lossvals = [];
else
    if sparsity_scale ==0
        [W,G,b,resid,dW,its_all,lossvals,thrs_EL,M] = wsindy_pde_RGLS_seq(lambda,gamma,Theta_pdx,lhs_ind,true_nz_weights,M_full,maxits);
    elseif sparsity_scale ==1
        [W,G,b,resid,dW,its_all,lossvals,thrs_EL,M] = wsindy_pde_RGLS_seq2(lambda,gamma,Theta_pdx,lhs_ind,true_nz_weights,M_full);
    end
    lambda_hat = lossvals(min(end,4),lossvals(min(end,4),:)>0);
    lambda_hat = lambda_hat(end);
end
ET_wsindy = toc;

tags_pde_G = tags_pde(~ismember(1:length(tags_pde),lhs_ind));
lib_list_G = lib_list(~ismember(1:length(tags_pde),lhs_ind),:);
end