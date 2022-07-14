nstate=1; %% choose solution component to plot
toggle_plot = 1; %% increase to speed up animation

figure(1); clf

Xcoordcell = get_Xcoords(sub_inds,xs_obs,m_x,m_t);

colormap(turbo(50));
zlims=[min(U_obs{min(nstate,end)}(:)) max(U_obs{min(nstate,end)}(:))];
t_ind = find(cellfun(@(x) ~isempty(x),Xcoordcell),1);

for j=1:toggle_plot:floor(length(xs_obs{end}))
    if dim==3
        imagesc(xs_obs{1},xs_obs{2},squeeze(U_obs{min(nstate,end)}(:,:,j))')
        set(gca,'YDir','normal')
        hold on
        if ~isempty(Xcoordcell{j})
            t_ind = j;
            viscircles(Xcoordcell{j},Xcoordcell{j}(:,1)*0+m_x*dx,'color','r','linewidth',1);
        else
            viscircles(Xcoordcell{t_ind},Xcoordcell{t_ind}(:,1)*0+m_x*dx,'color','r','linewidth',1);
        end
        hold off
        axis equal
        xlim(xs_obs{1}([1 end])); ylim(xs_obs{2}([1 end]));
        xlabel('$x$','interpreter','latex','fontsize',14)
        ylabel('$y$','interpreter','latex','fontsize',14)            
        title(num2str(xs_obs{3}(j)))
        drawnow
    elseif dim==2
        plot(xs_obs{1},U_obs{min(nstate,end)}(:,j)','.-')
        hold on
        if ~isempty(Xcoordcell{j})
            t_ind = j;
            for k=1:length(Xcoordcell{j})
                ind = find(Xcoordcell{j}(k)==xs_obs{1});
                plot(xs_obs{1}(1:2*m_x+1) + xs_obs{1}(max(mod(ind,length(xs_obs{1})),1))-xs_obs{1}(1),Cfs{1}{1}(1,:)/max(Cfs{1}{1}(1,:))*mean(abs(zlims))/2+zlims(1),'r-o');
            end
        else
            for k=1:length(Xcoordcell{t_ind})
                ind = find(Xcoordcell{t_ind}(k)==xs_obs{1});
                plot(xs_obs{1}(1:2*m_x+1) + xs_obs{1}(max(mod(ind,length(xs_obs{1})),1))-xs_obs{1}(1),Cfs{1}{1}(1,:)/max(Cfs{1}{1}(1,:))*mean(abs(zlims))/2+zlims(1),'r-o');
            end
        end
        hold off
        ylim(zlims)
        drawnow
    end
end

function X=get_Xcoords(sub_inds,xs_obs,m_x,m_t)

    dim = length(sub_inds);
    sub_dims = cellfun(@(x) length(x),sub_inds,'uni',0);
    num_tf = prod(cell2mat(sub_dims));
    
    subindgrid = cell(1,dim);
    [subindgrid{:}] = ndgrid(sub_inds{:});
    
    coords = zeros(num_tf,dim);
    for i=1:num_tf
        coord = cellfun(@(x) x(i),subindgrid);
        coords(i,:) = coord + [ones(1,dim-1)*m_x m_t];
    end
    
    times = unique(coords(:,end));
    Lt = length(times);
    X = cell(1,length(xs_obs{end}));
    for i=1:Lt
        tinds = find(coords(:,end)==times(i));
        for k=1:length(tinds)
            newp = [];
            for ll=1:dim-1
                newp = [newp xs_obs{ll}(max(mod(coords(tinds(k),ll),length(xs_obs{ll})),1))];
            end
            X{times(i)}=[X{times(i)};newp];
        end
    end
end
