function [U_obs,xs] = conv_movavg(U_obs,xs,w)

n = length(U_obs);
Ns = size(U_obs{1});
dim = length(Ns);
vec = cell(dim,1);

if ~isempty(w)
    if ~and(length(w)==1,w(1)==0)
        if length(w)==1
            w=repmat(w,dim,1);
        end
        for k=1:dim
            if w(k)~=0
                vec{k} = ones(1,2*w(k)+1)/(2*w(k)+1);
            end
        end
        for nn=1:n
            U_temp = U_obs{nn};
            vec_ffts = cell(dim,1);
            sub_inds = cell(dim,1);
            for k=1:dim
                if ~isempty(vec{k})
                    vec_ffts{k} = [zeros(1,Ns(k)-(2*w(k)+1)) vec{k}];
                    vec_ffts{k} = fft(vec_ffts{k});
                    sub_inds{k} = 1:Ns(k)-2*w(k)-1;
                end
            end
            U_obs{nn} = convNDfft(U_temp,vec_ffts,sub_inds,2);
        end
    end
end
for k=1:dim
    xs{k} = xs{k}(w(k)+1:size(U_obs{1},k)+w(k));
end

end