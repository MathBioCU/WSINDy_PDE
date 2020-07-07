%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate 'tags' for terms to include in model,
%%%%%%%%%%%% indices for true model terms, true model weights 
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [pdx_list,tags_pde,lib_list,lhs_ind,axi] = get_lib_tags(n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_pt,custom_remove,custom_add,axi_tags)

if isempty(custom_remove)
    custom_remove = ones(1,n+dim)*Inf;
end

[~,pdx_list,lib_list] = build_fcn_lib_tags(n,dim,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_pt);
lib_list = [lib_list;custom_add];
[tags_pde,lib_list] = build_str_tags(lib_list,dim,n,custom_remove,use_all_pt,lhs);

lhs_ind=zeros(size(lhs,1));
try
    for k=1:size(lhs,1)
        lhs_ind(k) = find(ismember(lib_list,lhs(k,:),'rows'));
    end
catch
    disp('ERROR: LHS not computed')
    return
end
if ~isempty(axi_tags)
    num_eqs = size(lhs,1);
    try 
        axi = zeros(size(lib_list,1)-1,num_eqs);
        for k=1:num_eqs
            axi_temp = tags2axi(axi_tags{k},lib_list);
            axi(:,k) = axi_temp([1:lhs_ind(k)-1 lhs_ind(k)+1:end]);
        end
    catch
        u_input = input('True terms missing from library, proceed anyway?'); 
        if u_input~=1 
            return;
        else
            axi = [];
        end
    end
end
end
