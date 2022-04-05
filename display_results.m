figind=0;
% close all
set(0,'DefaultFigurePosition',[4000 800 580 406])

%% plot basis fcn

if toggle_plot_basis_fcn
    figind = figind +1;
    figure(figind); 
    pos=[];
    plot_basis_fcn(Cfs_x,Cfs_t,m_x,dx,m_t,dt,max_dx,max_dt,pdx_list,pos,scales(end-n:end));
end

%% plot data

if toggle_plot_sol>0
for jj=1:min(length(toggle_plot_sol),length(U_obs))
    figind = figind +1;
    figure(figind); 
    colormap(turbo(50))
    if dim==2
        surf(xs_obs{1},xs_obs{2},U_obs{min(toggle_plot_sol(jj),length(U_obs))}', 'EdgeColor','none')
%         view([0 90])       
        zlabel('$u$','interpreter','latex','fontsize',14)
        set(gca, 'TickLabelInterpreter','latex','fontsize',14)
        xlabel('$x$','interpreter','latex','fontsize',14)
        ylabel('$t$','interpreter','latex','fontsize',14)
        xlim([xs_obs{1}(1) xs_obs{1}(end)])
        ylim([xs_obs{2}(1) xs_obs{2}(end)])
        colorbar
        title(['U(',num2str(min(toggle_plot_sol(jj),length(U_obs))),')'])
    elseif dim==3
        for j=1:plotgap:floor(length(xs_obs{3}))
            for i=1:length(toggle_plot_sol)
                subplot(length(toggle_plot_sol),1,i)
                surf(1:length(xs_obs{1}),1:length(xs_obs{2}),squeeze(U_obs{min(toggle_plot_sol(jj),length(U_obs))}(:,:,j))', 'EdgeColor','none')
                xlabel('$x$','interpreter','latex','fontsize',14)
                ylabel('$y$','interpreter','latex','fontsize',14)            
                title(['U(',num2str(min(toggle_plot_sol(jj),length(U_obs))),')'])
                view([0 90])     
                title(num2str(xs_obs{3}(j)))
                colorbar
                caxis([min(reshape(U_obs{toggle_plot_sol(jj)},[],1)) max(reshape(U_obs{toggle_plot_sol(jj)},[],1))])
            end
            drawnow
        end
    end
end
end 

%% plot loss fcn (MSTLS)

if size(lossvals,2)>1 && toggle_plot_loss
    figind = figind +1;
    figure(figind); 
    loglog(lossvals(2,:),lossvals(1,:),'o-')
    xlabel('$\lambda$','interpreter','latex','fontsize',14)
    ylabel('$\mathcal{L}$','interpreter','latex','fontsize',14)
    set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    xtick =10.^linspace(log10(lossvals(2,1)),log10(lossvals(2,end)),4);
    xticks(xtick)
    xticklb = num2str(log10(xtick)');
    xticklabels(strcat('$10^{',xticklb(:,1:min(4,end)),'}$'))
end

%% plot data fft

if toggle_plot_fft>0
    Cfs = {{Cfs_x,Cfs_t}};
    figind = figind +1;
    figure(figind)
    coords=1:dim;
    for j=1:dim-1
        coords=[coords;circshift(coords(end,:),-1)];
    end
    for j=1:dim
        coordsj = coords(j,:);
        dd=coordsj(1);
        Ufft = abs(fft(permute(U_obs{1},coordsj)));
        Ufft = reshape(Ufft,size(Ufft,1),[]);
        Ufft = mean(Ufft(floor(end/2):end,:),2);
        L = length(Ufft)-1;
        ks = -L:L;
        Ufft = [Ufft; flipud(Ufft(1:end-1))]/max(Ufft);
        subplot(dim,1,j)
            semilogy(ks,Ufft)
            hold on
            if dd<dim
                Cfs_ffts = fft([zeros(1,length(xs_obs{dd})-2*m_x-1) Cfs{1}{1}(1,:)]);
            else
                Cfs_ffts = fft([zeros(1,length(xs_obs{dd})-2*m_t-1) Cfs{1}{2}(1,:)]);
            end
            Cfs_ffts=abs(Cfs_ffts(floor(end/2):end));
            Cfs_ffts=[Cfs_ffts fliplr(Cfs_ffts(1:end-1))];
            Cfs_ffts=Cfs_ffts/max(Cfs_ffts);
            semilogy(ks,Cfs_ffts)
            if exist('corners','var')
                k = corners{1}{dd}(2);
                semilogy([-k k],Ufft(L+1-k)*[1 1],'o','markersize',12)
            end    
            hold off     
            ylim([min(Ufft)*0.1 max(Ufft)])
            legend({'$\mathcal{F}(U_d)$','$\mathcal{F}(\phi_d)$','$k_d^*$'},'interpreter','latex','fontsize',14)
            title(['coord',num2str(dd)])
    end
    xlabel('k')
end
