%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: get test function values on scaled reference grid,
%%%%%%%%%%%% generate indices for query points (sub_inds)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [Cfs_x,Cfs_t,dx,dt,p_x,p_t,sub_inds] = get_testfcn_weights(dims,xs,max_dx,max_dt,m_x,m_t,s_x,s_t,tols,phi_class)

dim = length(dims);

if or(2*m_x+1>min(dims(1:end-1)),2*m_t+1>min(dims(end)))
    disp('ERROR: Test function not compactly supported')
    return
end

sub_inds = cell(1,dim);
ss = [repmat(s_x,1,dim-1) s_t];
mm = [repmat(m_x,1,dim-1) m_t];

for j=1:dim
    N = dims(j);
    m = mm(j);
    s = ss(j);
    sub_inds{j} = 1:s:N-2*m;
%     divs = divisors(N-2*m-1);
%     [~,s] = min(abs(divs - s));    
%     sub_inds{j} = 1:divs(s):N-2*m;
end

dx = xs{1}(2)-xs{1}(1);
dt = xs{end}(2)-xs{end}(1);
[Cfs_x,p_x] = phi_int_weights(m_x,max_dx,tols(1),phi_class);
[Cfs_t,p_t] = phi_int_weights(m_t,max_dt,tols(2),phi_class);

end
