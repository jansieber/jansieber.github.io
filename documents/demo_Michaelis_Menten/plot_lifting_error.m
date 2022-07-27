%% Plot error of implicit lifting scheme
figure(4);clf;
Phi_err=abs(Phi_tskip_fit(1:nt,:)-Phi_star_fit(ones(1,nt),:));
anchor_ind=1;
p_err=[];
lgtext={};
npl=0;
for j=3:-1:1
    npl=npl+1;
    p_err(npl)=semilogy(tskip,Phi_err(:,j),'o-','color',clr(j,:),ldeco{:});
    hold on
    lgtext{npl}=sprintf('d^%d/dx^%d   (Phi_tskip-Phi_star)',3-j,3-j);
end
npl=npl+1;
asymptote=log(Phi_err(anchor_ind,3))+(tskip-tskip(anchor_ind))*(-d_tr)-1;
p_err(npl)=semilogy(tskip,exp(asymptote),'k--',ldeco{:});
lgtext{npl}=sprintf('asymptote exp(-dtr*tskip)');
legend(p_err,lgtext,'Interpreter','none');
grid on
ylim([min(Phi_err(:)),max(Phi_err(:))]);
xlabel('tskip');
title('Difference to asymptotic expansion in epsilon');
