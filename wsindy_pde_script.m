%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: script for recoverying PDE systems
%%%%%%%%%%%% 
%%%%%%%%%%%% pde_num selects a PDE system from the list pde_names
%%%%%%%%%%%% noise_ratio sets the signal-to-noise ratio (L2 sense)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

pde_num = 3;
noise_ratio=0.04;
filename = [];
toggle_plot_basis_fcn = 0;

% Load clean data
pde_names = {'KS.mat','NLS.mat','Sine_Gordon.mat','rxn_diff.mat','Nav_Stokes.mat'};
load(['datasets/',pde_names{pde_num}])
clc;

% Add white noise
n = length(U_exact);
dims = size(squeeze(U_exact{1}));
dim = length(dims);

if noise_ratio>0
    U_obs = cell(n,1);
    stdvs = zeros(1,n);
    noise = cell(1,n);
    snr = zeros(1,n);
    for k=1:n
        stdvs(k) = rms(U_exact{k}(:))^2;
    end
    for j=1:n
        rng('shuffle');
        sigma = noise_ratio*sqrt(stdvs(j));
        noise{j} = normrnd(0,sigma,size(U_exact{j}));
        snr(j) = norm(noise{j}(:))/norm(U_exact{j}(:));
        U_obs{j} = U_exact{j} + noise{j};
        disp(['SNR=',num2str(snr(j))]); 
    end
else
    U_obs = U_exact;
end

% Identify PDE
tic,
[pdx_list,tags_pde,lib_list,lhs_ind,axi] = get_lib_tags(n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_pt,custom_remove,custom_add,axi_tags);
[Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds] = get_testfcn_weights(dims,x,t,max_dx,max_dt,supp_phi_x,supp_phi_t,s_x,s_t);
Theta_pdx = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,supp_phi_x,supp_phi_t,dx,dt,sub_inds,dim);
[W,G,b,resid,dW] = wsindy_pde_RGLS(lambda,gamma,Theta_pdx,lhs_ind,axi);
ET_wsindy = toc;

% Print results to filename, or command window if filename = []
print_results(W,G,b,resid,dW,filename,dims,polys,trigs,max_dx,max_dt,lambda,gamma,lhs_ind,tags_pde,supp_phi_x,supp_phi_t,p_x,p_t,s_x,s_t,ET_wsindy)
if toggle_plot_basis_fcn
    plot_basis_fcn(Cfs_x,Cfs_t,supp_phi_x,dx,supp_phi_t,dt,max_dx,max_dt,pdx_list,[1000 5 1400 700]);
end
