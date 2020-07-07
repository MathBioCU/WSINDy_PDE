%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: generate true model weights 
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function axi = tags2axi(axi_tags,lib_list)
  m = size(lib_list,1);
  axi = zeros(m,1);
  [~,loc] = ismember(axi_tags(:,1:end-1),lib_list,'rows');
  axi(loc) = axi_tags(:,end);
end