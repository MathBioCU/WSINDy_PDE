warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
[m,n] = size(W);

str_wsindy = cell(n,1);
for k=1:n
    str_wsindy{k} = print_pde(W(:,k),tags_pde_G,tags_pde{lhs_ind(k)});
end
if ~isempty(axi)
    Tps = tpscore({W},axi);
else
    Tps=NaN;
end
dW = wnorm({W},axi,Inf);

try
    nu_learned = W(ismember(tags_pde_G,{'u^{1}_{lap}','u^{1}_{xx}'}));
catch
    nu_learned = 0;
end

if ~isequal(print_loc,0)
    if ~isequal(print_loc,1)
        print_loc = fopen(print_loc,'a');
    end

    for k=1:n
        fprintf(print_loc,['\nRecovered PDE: ',str_wsindy{k}]);
        fprintf(print_loc,'\nRelative Res: ||b-G*W||_2/||b||_2 = %.2e',norm(resid(:,k)));
        if ~isempty(dW)
            fprintf(print_loc,'\nMax Weight Error: max|W-W_{true}| = %.2e\n', dW(k));
        end
    end
    fprintf(print_loc,'TP Score = %1.2f\n', Tps);
    fprintf(print_loc,'      \n');
    fprintf(print_loc,'polys = ');
    fprintf(print_loc,'%u ',polys);
    fprintf(print_loc,'\ntrigs = ');
    fprintf(print_loc,'%u ',trigs);
    fprintf(print_loc,'\nMax derivs [t x] = ');
    fprintf(print_loc,'%u ',[max_dt max_dx]);
    fprintf(print_loc,'\n[m_x m_t] = ');
    fprintf(print_loc,'%u ',[m_x m_t]);
    fprintf(print_loc,'\n[s_x s_t] = ');
    fprintf(print_loc,'%u ',[diff(sub_inds{1}(1:min(length(sub_inds{1}),2))') diff(sub_inds{end}(1:min(length(sub_inds{end}),2))')]);
    fprintf(print_loc,'\n[p_x p_t] = ');
    pps=zeros(1,2);
    ps=[p_x p_t];
    for i=1:2
        if isequal(phi_class,1)
            pps(i)=ps(i);
        elseif isequal(phi_class,2)
            pps(i)=1/ps(i);
        else
            pps(i)=NaN;
        end
    end
    fprintf(print_loc,'%u ',pps);
    fprintf(print_loc,'\n scales = ');
    fprintf(print_loc,'%.2e ',scales) ;
    fprintf(print_loc,'\n      \n');
    fprintf(print_loc,'Size of dataset = ');
    fprintf(print_loc,'%u ',dims);
    fprintf(print_loc,'\nSize G = ');
    fprintf(print_loc,'%u ',size(G));
    if gamma >0
        fprintf(print_loc,'\nCond G = %.2e',cond([G;gamma*norm(G)*eye(m)]));
    else
        fprintf(print_loc,'\nCond G = %.2e',cond(G));
    end
    fprintf(print_loc,'\n[lambda_hat gamma] = ');
    fprintf(print_loc,'%.3e ',[lambda_hat gamma*norm(G)]);
    fprintf(print_loc,'\n[sigma_NR sigma] = ');
    fprintf(print_loc,'%.3e ',[sigma_NR sigma]);
    fprintf(print_loc,'\n      \n');
    fprintf(print_loc,'Elapsed time WSINDy = %4.4f \n',ET_wsindy);
    fprintf(print_loc,'STLS its = ');
    fprintf(print_loc,'%u ',its_all);
    fprintf(print_loc,'\n ');

    if ~all(print_loc==1)
        fclose(print_loc);
    end
end