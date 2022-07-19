function [xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims)

if ~isempty(coarsen_data)
    dim = length(size(U_obs{1}));
    n = length(U_obs);
    inds = cell(1,dim);
    for j=1:dim
        inds{j} = 1+floor(coarsen_data(j,1)*dims(j)):ceil(coarsen_data(j,2)):ceil(coarsen_data(j,3)*dims(j));
        xs_obs{j} = xs_obs{j}(inds{j});
    end
    for j=1:n
        U_obs{j} = U_obs{j}(inds{:});
    end
end


end