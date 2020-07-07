%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: print results of WSINDy_PDE to 'filename', with
%%%%%%%%%%%% filename = [] leading to results printed in command window
%%%%%%%%%%%%  
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function print_results(W,G,b,resid,dW,filename,dims,polys,trigs,max_dx,max_dt,lambda,gamma,lhs_ind,tags_pde,supp_phi_x,supp_phi_t,p_x,p_t,s_x,s_t,ET_wsindy)

[m,n] = size(W);

if isempty(filename)
    filename = 1; 
else
    filename = fopen(filename,'a');
end

for k=1:n
    tags_pde_rdx = tags_pde([1:lhs_ind(k)-1 lhs_ind(k)+1:end]);
    str_wsindy = print_pde(W(:,k),tags_pde_rdx,tags_pde{lhs_ind(k)});
    fprintf(filename,['Recovered PDE: ',str_wsindy]);
    fprintf(filename,'\nRel. Resid: ||b-G*W||_2 = %4e',norm(resid(:,k))/norm(b(:,k)));
    fprintf(filename,'\nMax Weight Error: max|W-W_{true}| = %6e\n', dW{1}(k));
end

fprintf(filename,'      \n');
fprintf(filename,'polys = ');
fprintf(filename,'%u ',polys);
fprintf(filename,'\ntrigs = ');
fprintf(filename,'%u ',trigs);
fprintf(filename,'\nMax derivs [t x] = ');
fprintf(filename,'%u ',[max_dt max_dx]);
fprintf(filename,'\n[m_x m_t] = ');
fprintf(filename,'%u ',[supp_phi_x supp_phi_t]);
fprintf(filename,'\n[s_x s_t] = ');
fprintf(filename,'%u ',[s_x s_t]);
fprintf(filename,'\n[p_x p_t] = ');
fprintf(filename,'%u ',[p_x p_t]) ;

fprintf(filename,'\n      \n');
fprintf(filename,'Size of dataset = ');
fprintf(filename,'%u ',dims);
fprintf(filename,'\nSize G = ');
fprintf(filename,'%u ',size(G{k}));
if gamma >0
    fprintf(filename,'\nlog_{10}(Cond G) = %4.2f',log10(cond([G{k};gamma*eye(m)])));
else
    fprintf(filename,'\nlog_{10}(Cond G) = %4.2f',log10(cond(G{k})));
end
fprintf(filename,'\n[lambda gamma] = ');
fprintf(filename,'%u ',[lambda gamma]);


fprintf(filename,'\n      \n');
fprintf(filename,'Elapsed time WSINDy = %4f \n',ET_wsindy);

if ~all(filename==1)
    fclose(filename);
end

end