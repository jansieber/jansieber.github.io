%% Plot fiber image and vertical line image
pl=@(x,y,s,cl)plot(x,y,s,'color',clr(cl,:),ldeco{:});
figure(2);clf;
hold on;
cltraj=(ones(2,3)+clr([2,4],:))/2;
plot(squeeze(Mg_pre(1,:,:)),squeeze(Mg_pre(2,:,:)),'-','color',cltraj(1,:),ldeco{:});
plot(squeeze(Mvert(1,:,:)),squeeze(Mvert(2,:,:)),'-','color',cltraj(2,:),ldeco{:});
pV=pl(zeros(1,ny)+vw0(1),vw_g_pre(2,:),'+-',4);
pgp=pl(vw_g_pre(1,:),vw_g_pre(2,:),'o-',2);
pVe=pl(squeeze(Mvert(1,end,:)),squeeze(Mvert(2,end,:)),'x',4);
pgpe=pl(squeeze(Mg_pre(1,end,:)),squeeze(Mg_pre(2,end,:)),'s',2);
pbase=plot(vw0(1),vw0(2),'ks','markerfacecolor','k',ldeco{:});
legend([pbase(1),pgp(1),pV(1),pgpe(1),pVe(1)],...
    {sprintf('(%s0,%s0)',coords{1},coords{2}),'stable fiber gpre',...
    'vertical line','M(t0,gpre)','M(t0,vertical line)'},...
    'location','NorthWest')
xlabel(coords{1});
ylabel(coords{2});
title('Trajectories starting from same stable fiber');
set(gca,'box','on');
hold off
