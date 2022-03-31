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

function [axi,tags_pde,lib_list,pdx_list,lhs_ind,Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds,scales,M_full,Theta_pdx] = wsindy_pde_fun(U_obs,xs,true_nz_weights,...
    lhs,max_dx,max_dt,polys,trigs,custom_add,custom_remove,use_all_dt,use_cross_dx,...
    toggle_scale,m_x,m_t,s_x,s_t,tols,phi_class)

    
n = length(U_obs);
dims = size(squeeze(U_obs{1}));
dim = length(dims);

tic,
[tags_pde,lib_list,pdx_list,lhs_ind,axi] = get_lib_tags(n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_dt,custom_remove,custom_add,true_nz_weights);

max_dx=max(reshape(lib_list(:,n+1:end-1),[],1));
max_dt=max(lib_list(:,end));
polys = unique(real(reshape(lib_list(:,1:n),[],1)));

[Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds] = get_testfcn_weights(dims,xs,max_dx,max_dt,m_x,m_t,s_x,s_t,tols,phi_class);
if length(toggle_scale) == 1
    if toggle_scale > 0
        scale_u = zeros(1,n);
        if toggle_scale == 1
            scale_u_fcn = @(v) norm(v(:)/norm(v(:),1)^(1/max(polys)),max(polys))^(max(polys)/(max(polys)-1));
        elseif toggle_scale == 2
            scale_u_fcn = @(v) min(max(norm(v(:)/norm(v(:),2)^(1/max(polys)),2*max(polys))^(max(polys)/max(max(polys)-1,1)),eps),1/eps);
        elseif toggle_scale == Inf
            scale_u_fcn = @(v) norm(v,inf)/(10^(1/max(polys)));
        end
        for k=1:n 
            scale_u(k) = scale_u_fcn(U_obs{k}(:));
        end
    else
        scale_u = [];
    end
else
    scale_u = toggle_scale;
end
if isempty(scale_u)
    scales = [];
    M_full = [];
else
    [scales,M_full] = get_scales(scale_u,p_x,m_x,dx,max_dx,p_t,m_t,dt,max_dt,lib_list,dim,phi_class,Cfs_x,Cfs_t);
end

Theta_pdx = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,m_x,m_t,dx,dt,sub_inds,dim,scales);
end
