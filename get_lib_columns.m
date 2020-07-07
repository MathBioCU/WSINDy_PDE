%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: compute G and b, listed together in Theta_pdx, 
%%%%%%%%%%%% using convolutions over separable test functions via convNDfft
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function Theta_pdx = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,supp_phi_x,supp_phi_t,dx,dt,sub_inds,dim)

Theta_pdx = [];
ind = 1;

while ind<size(lib_list,1)+1
    tags = lib_list(ind,1:n);
    if isreal(tags)
        fcn = U_obs{1}.^tags(1);
        for k=2:n
            fcn = fcn.*(U_obs{k}.^tags(k));
        end
    else
        ind_freq = find(imag(tags(1:n)));
        freq = sum(imag(tags(1:n)));        
        if freq<0
            fcn = sin(abs(freq)*U_obs{ind_freq});
        else
            fcn = cos(abs(freq)*U_obs{ind_freq});
        end
    end
    
    while all(lib_list(ind,1:n) == tags)
        test_conv_cell = {};
        for k=1:dim-1
            test_conv_cell{k} = Cfs_x(lib_list(ind,n+k)+1,:)' * (supp_phi_x*dx)^(-lib_list(ind,n+k)) * dx;
        end
        test_conv_cell{dim} = Cfs_t(lib_list(ind,end)+1,:)' * (supp_phi_t*dt)^(-lib_list(ind,end)) * dt;
        fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds);
        Theta_pdx(:,ind) = fcn_conv(:);
        ind = ind+1;
        if ind > size(lib_list,1)
            break
        end
    end
end
end