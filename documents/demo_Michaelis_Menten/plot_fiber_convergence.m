%% plot convergence of points inside fiber to each other
nrm=@(x,i)sqrt(sum(x.^2,i));
figure(3);clf
pgp=semilogy(tr,nrm(max(Mg_pre,[],3)-min(Mg_pre,[],3),1),'color',clr(2,:),ldeco{:});
hold on
pV=semilogy(tr,nrm(max(Mvert,[],3)-min(Mvert,[],3),1),'color',clr(4,:),ldeco{:});
legend([pV(1),pgp(1)],...
    {'radius of M(t;gpre)','radius of M(t;vertical line)'},...
    'location','East')
xlabel('time t');
ylabel('radius');
title('Radius of ball of phase space points over time');
set(gca,'box','on');
hold off
