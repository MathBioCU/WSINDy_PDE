function [xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims)

if ~and(all(reshape(coarsen_data(:,1:2)==1,[],1)),all(coarsen_data(:,end)==dims'))
    nstates = length(U_obs);
    dim = length(size(U_obs{1}));
    inds = cell(1,dim);
    for j=1:dim
        N = length(xs_obs{j});
        inds{j} = coarsen_data(j,1):coarsen_data(j,2):coarsen_data(j,3);
        xs_obs{j} = xs_obs{j}(inds{j});
    end
    for j=1:nstates
        U_obs{j} = U_obs{j}(inds{:});
    end
end

end