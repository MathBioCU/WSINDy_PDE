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

%% Load data

% clc; 
% clear all; 
% close all;

pde_num = 1;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','Sine_Gordon.mat','rxn_diff.mat','Nav_Stokes.mat','porous.mat','AC_temp'};
% load(['datasets/',pde_names{pde_num}])

%% Subsample data (if desired)

U_obs = U_exact;
xs_obs = xs;
dims = size(U_obs{1});
dim = length(dims);

coarsen_data = [ones(dim,2) dims'];
% coarsen_data(end,end) = floor(coarsen_data(end,end)/2); %%% uncomment to use only first 1/2 of timepoints
[xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims);

%% Add noise

sigma_NR = 0;
noise_dist = 0; 
noise_alg = 0;
rng('shuffle');
rng_seed = rng().Seed; 
 
rng(rng_seed);
[U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,rng_seed,0);

%% Set hyperparameters

%---------------- weak discretization

% m_x = min(10,floor((length(xs_obs{1})-1)/2));
% m_t = min(10,floor((length(xs_obs{end})-1)/2));
% s_x = max(floor(length(xs_obs{1})/25),1);
% s_t = max(floor(length(xs_obs{end})/25),1);
% phi_class = 1;
% tau = 10^-10;
% tauhat = 0; 
tauhat_inds = 1;
% toggle_scale = 0;

%---------------- model library

% max_dx = 6;
% max_dt = 1;
% polys = 0:4;
% trigs = [];
use_all_dt = 0;
% use_cross_dx = 0;

% custom_add = [];
% custom_remove = [];

%---------------- find test function hyperparameters from tau,tauhat

if tauhat > 0
    [m_x,m_t,p_x,p_t,sig_est,corners_all] = findcorners(U_obs(tauhat_inds),xs_obs,tau,tauhat,max_dx,max_dt,phi_class);
    tols = [-p_x -p_t];
else
    tols = [tau tau];
end

% Build Library

n = length(U_obs);
dims = size(U_obs{1});
dim = length(dims);

[axi,tags_pde,lib_list,pdx_list,lhs_ind,Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx] = wsindy_pde_fun(U_obs,xs,true_nz_weights,...
    lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_dt,use_cross_dx,...
    toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class);

%% Solve Sparse Regression Problem

%---------------- regularized least squares solve
lambda = 10.^(linspace(-4,0,50));
%lambda = [];
gamma = 0;
maxits = 5;
sparsity_scale = 0;                     % 0, enforce sparsity on original data; 1, enforce on rescaled data

[W,G,b,resid,dW,its_all,thrs_EL,M,lambda_hat,lossvals,ET_wsindy,tags_pde_G,lib_list_G] = wsindy_pde_solve(lambda,gamma,Theta_pdx,lhs_ind,axi,M_full,maxits,tags_pde,lib_list,sparsity_scale);

% Display results

print_loc = 1;
toggle_plot_basis_fcn = 0; 
toggle_plot_sol = 3; 
toggle_plot_loss = 0; 

if print_loc~=0
    str_wsindy  = print_results(W,G,resid,dW,print_loc,dims,polys,trigs,max_dx,max_dt,lambda_hat,gamma,lhs_ind,tags_pde,m_x,m_t,p_x,p_t,sub_inds,scales,ET_wsindy,its_all,sigma_NR,sigma,axi);
end

if toggle_plot_basis_fcn
    figure(1);
    plot_basis_fcn(Cfs_x,Cfs_t,m_x,dx,m_t,dt,max_dx,max_dt,pdx_list,[1000 5 1400 700],scales(end-nstates:end));
end

if toggle_plot_sol~=0
    figure(2); set(gcf,'units','normalized','outerposition',[0.7 0.5 0.3 0.45])
    set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    xlabel('$x$','interpreter','latex','fontsize',14)
    ylabel('$t$','interpreter','latex','fontsize',14)
    if dim==2
        surf(1:length(xs_obs{1}),1:length(xs_obs{2}),U_obs{min(end,toggle_plot_sol)}', 'EdgeColor','none')
        colorbar
        view([0 90])       
        zlabel('$u$','interpreter','latex','fontsize',14)
        set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    elseif dim==3
        for j=1:length(xs_obs{3})
            surf(1:length(xs_obs{1}),1:length(xs_obs{2}),squeeze(U_obs{min(end,toggle_plot_sol)}(:,:,j))', 'EdgeColor','none')
            colorbar
            view([0 90])       
            drawnow
        end
    end
end

if ~isempty(lossvals) && toggle_plot_loss
    figure(3); set(gcf,'units','normalized','outerposition',[0.7 0.05 0.3 0.45])
    loglog(lossvals(2,:),lossvals(1,:),'o-')
    xlabel('$\lambda$','interpreter','latex','fontsize',14)
    ylabel('$\mathcal{L}$','interpreter','latex','fontsize',14)
    set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    xtick =10.^linspace(log10(lossvals(2,1)),log10(lossvals(2,end)),4);
    xticks(xtick)
    xticklb = num2str(log10(xtick)');
    xticklabels(strcat('$10^{',xticklb(:,1:min(4,end)),'}$'))
end

