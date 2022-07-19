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
close all;
clear all;

pde_num = 9;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','Sine_Gordon.mat','rxn_diff.mat','Nav_Stokes.mat','porous.mat','sod.mat'};
dr = ['datasets/',pde_names{pde_num}];
% dr = ['/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/',pde_names{pde_num}];

try
    U_obs = U_exact;
    xs_obs = xs;
catch
    load(dr,'U_exact','xs','lhs','true_nz_weights')
    U_obs = U_exact;
    xs_obs = xs;

    eq = 2;
    lhs = lhs(eq,:);
    true_nz_weights = true_nz_weights(eq);
end

dims = size(U_obs{1});
dim = length(dims);
n = length(U_obs);

%% Subsample data (if desired)

coarsen_data = repmat([0 1 1],dim,1);

%%% set row d of coarsen_data to [initial_frac inc final_frac] to subsample dth coordinate to 
%%% start at index initial_frac*L, where L is the number of points in dth coordinate
%%% end at index final_frac*L,
%%% skip every inc gridpoint.

coarsen_data(1,:) = [0 8 1];  % every 8th point in x
coarsen_data(2,:) = [0.5 2 1]; % start 50% of the way through time series, take every other point in time

[xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims);
dims = cellfun(@(x) length(x), xs_obs);

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
%%% phi_class = 1 for piecewise polynomial test function, 2 for Gaussian
phi_class = 1;          

%%% set convolution query point spacing:
s_x = 3;  %max(floor(length(xs_obs{1})/50),1);
s_t = 3;  %max(floor(length(xs_obs{end})/50),1);

%%% set reference test function parameters using spectrum of data:
tauhat = 0;      %%% if tauhat<=0, explicit vals for m_x,m_t,p_x,p_t used. 
tau = 10^-10;

%%% set reference test function parameters explicitly:
m_x = 15;
m_t = 15;
p_x = 5;
p_t = 5;

%%% toggle rescale state variables and spatiotemporal coordinates
toggle_scale = 2;

%---------------- model library
max_dx = 1;
max_dt = 1;
polys = [0:3];
trigs = [];
use_all_dt = 0;
use_cross_dx = 0;
custom_add = [];
custom_remove = {@(mat) mat(:,3)>1};
% lhs = [1 1 0 0 1];
% true_nz_weights = {};

%% set plots

toggle_plot_basis_fcn = 0;
toggle_plot_sol =  0;
plotgap = 0;
toggle_plot_loss = 1;
toggle_plot_fft = 0;

%% Build Linear System

%---------------- find test function hyperparams using Fourier spectrum of U
if tauhat > 0
    tauhat_inds = 1:n;
    [m_x,m_t,p_x,p_t,sig_est,corners] = findcorners(cellfun(@(x) x.^1, U_obs(tauhat_inds), 'uni',0),xs_obs,tau,tauhat,max_dx,max_dt,phi_class);
else
    m_x = min(m_x,floor((length(xs_obs{1})-1)/2));
    m_t = min(m_t,floor((length(xs_obs{end})-1)/2));
end
tols = [-p_x -p_t];

%---------------- build linear system
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

print_loc = 1;
get_results;
display_results;
