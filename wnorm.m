function dW = wnorm(Ws,axi,p)
    if ~isempty(axi)
        eq=size(axi,2);
        if isequal(class(Ws),'cell')
            M=length(Ws);
            dW = zeros(M,eq);
            for m=1:M
                W = Ws{m};
                dW(m,:) = vecnorm(W-axi,p,1)./vecnorm(axi,p,1);
            end
        elseif isequal(class(Ws),'double')
            dims = size(Ws);
            M=dims(end);
            dW = zeros(M,eq);
            for m=1:M
                if length(dims)==3
                    W = Ws(:,:,m);
                else
                    W = Ws(:,m);
                end
                dW(m,:) = vecnorm(W-axi,p,1)./vecnorm(axi,p,1);
            end
        end
    else
        dW=[];
    end
end
