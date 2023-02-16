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

clc;
% close all;
clear all;

pde_num = 3;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','Sine_Gordon.mat','rxn_diff.mat','Nav_Stokes.mat','porous.mat','sod.mat'};
dr = ['datasets/',pde_names{pde_num}];
% dr = ['/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/',pde_names{pde_num}];

try
    %%% don't reload data if already loaded
    U_obs = U_exact;
    xs_obs = xs;
catch
    load(dr);
    U_obs = U_exact;
    xs_obs = xs;
end

%%% select subset of equations
eq = 1:length(U_obs);
lhs = lhs(unique(min(eq,end)),:);
true_nz_weights = true_nz_weights(unique(min(eq,end)));

dims = size(U_obs{1});
dim = length(dims);
n = length(U_obs);

%% Subsample data (if desired)

coarsen_data = repmat([0 1 1],dim,1);

%%% set row d of coarsen_data to [initial_frac inc final_frac] to subsample dth coordinate to 
%%% start at index initial_frac*L, where L is the number of points in dth coordinate
%%% end at index final_frac*L,
%%% skip every inc gridpoint.

coarsen_data(1:dim-1,:) = repmat([0 1 1],dim-1,1); 


[xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims);
dims = cellfun(@(x) length(x), xs_obs);

%% Add noise

sigma_NR = 0;
noise_dist = 0; 
noise_alg = 0;
rng(1);
rng_seed = rng().Seed;
 
rng(rng_seed);
[U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,rng_seed,0);

%% Set hyperparameters 

use_presets = 0;
if ~use_presets
    %---------------- weak discretization
    %%% phi_class = 1 for piecewise polynomial test function, 2 for Gaussian
    phi_class = 1;

    %%% set query point spacing by maximum row restriction
    K_max = 5000;

    %%% set reference test function parameters using spectrum of data:
    tauhat = 2;      %%% if tauhat<=0, explicit vals for m_x,m_t,p_x,p_t used. 
    tau = 10^-10;

    %%% set query point spacing manually
    s_x = 4; 
    s_t = 4; 

    %%% set reference test function parameters explicitly:
    m_x = 36;
    m_t = 34;
    p_x = 8;
    p_t = 9;
    
    %%% toggle rescale state variables and spatiotemporal coordinates
    toggle_scale = 2;
    
    %---------------- model library
    max_dx = 3+max(max(true_nz_weights{1}(:,n+1:n+dim-1))); 
    polys = 0:3+max(max(true_nz_weights{1}(:,1:n)));
    max_dt = max(lhs(:,end));
    trigs = 1:max(abs(reshape(imag(true_nz_weights{1}),[],1)));
    use_all_dt = 0;
    use_cross_dx = any(sum(logical(true_nz_weights{1}(:,n+1:n+dim-1)),2)>1);
    custom_add = [];
    custom_remove = {};%{@(mat) mat(:,3)>1};
end
% lhs = [1 1 0 0 1];
% true_nz_weights = {};

%% Build Linear System

%---------------- find test function hyperparams using Fourier spectrum of U
if tauhat > 0
    tauhat_inds = 1;
    [m_x,m_t,p_x,p_t,sig_est,corners] = findcorners(cellfun(@(x) x.^1, U_obs(tauhat_inds), 'uni',0),xs_obs,tau,tauhat,max_dx,max_dt,phi_class);
else
    m_x = min(m_x,floor((length(xs_obs{1})-1)/2));
    m_t = min(m_t,floor((length(xs_obs{end})-1)/2));
end
tols = [-p_x -p_t];

%---------------- build linear system
if K_max>0
    s_x = max(ceil((length(xs_obs{1})-2*m_x)/K_max^(1/dim)),1);
    s_t = max(ceil((length(xs_obs{end})-2*m_t)/K_max^(1/dim)),1);
end

[axi,tags_pde,lib_list,pdx_list,lhs_ind,Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx] = wsindy_pde_fun(U_obs,xs_obs,true_nz_weights,...
    lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_dt,use_cross_dx,...
    toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class);

%% Solve Sparse Regression Problem

lambda = 10.^(linspace(-4,0,100));
gamma = 0;
maxits = Inf;

%%% sparsity_scale =  0 enforces sparsity on original data; = 1 enforces on rescaled data
sparsity_scale = 0;

[W,G,b,resid,dW,its_all,thrs_EL,M,lambda_hat,lossvals,ET_wsindy,tags_pde_G,lib_list_G] = wsindy_pde_solve(lambda,gamma,Theta_pdx,lhs_ind,axi,M_full,maxits,tags_pde,lib_list,sparsity_scale);

%% Display results

toggle_plot_basis_fcn = 1;
toggle_plot_sol =  1;
plotgap = 5;
toggle_plot_loss = 1;
toggle_plot_fft = 1;

print_loc = 1;
get_results;
display_results;
