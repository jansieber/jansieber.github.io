%% Plot basic geometry
ny=100;
y_dev=linspace(-0.2,0.7,ny);
nfib=7;
xfibbase=linspace(min(xrange),max(xrange),nfib);
vwbd=[min(h_eps(1,:)),max(h_eps(1,:));min(h_eps(2,:)),L_level+0.1]; 
fibcol=[1,1,1]*0.7;
vtraj0=[];
figure(1);clf;hold on
for i=nfib:-1:1
    vwf=rotmat*fiber_graph_xy(xfibbase(i),y_dev);
    fsel= vwf(1,:)>=vwbd(1,1) & vwf(1,:)<=vwbd(1,2) &...
         vwf(2,:)>=vwbd(2,1) & vwf(2,:)<=vwbd(2,2);
    [~,imin]=min(abs(vwf(2,:)-L_level));
    vtraj0(i)=vwf(1,imin);
    ptmp=plot(vwf(1,fsel),vwf(2,fsel),'-','color',fibcol,ldeco{:});
    if ~isempty(ptmp)
      pfib=ptmp;
    end
end
pLdom=plot(vwbd(1,:),L_level+[0,0],'-','color',clr(2,:),ldeco{:});
% slow manifold
mfcol=clr(5,:);
ph=plot(h_eps(1,:),h_eps(2,:),'-','linewidth',3,'color',mfcol);
% Example trajectories
vtraj0=vtraj0(vtraj0>=vwbd(1,1) & vtraj0<=vwbd(1,2));
ntraj=length(vtraj0);
trajcol=clr(4,:);
for i=ntraj:-1:1
    [~,~,vwt]=M(2/sqrt(epsilon),[vtraj0(i);L_level]);
    indsel= vwt(1,:)<=vwbd(1,2);
    ptraj=plot(vwt(1,indsel),vwt(2,indsel),'-','color',trajcol,ldeco{:});
end
legend([ph,pLdom,ptraj(1),pfib(1)],...
    {'slow manifold','dom L','trajectories','stable fibers gpre'},...
    'location','NorthWest')
title('Geometry in phase space')
xlabel(coords{1});
ylabel(coords{2});
hold off
