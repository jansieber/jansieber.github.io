%% Convert gain as contructed by Brunovsky to true gain
% Then one can use |u=K[xtilde-x]|.
%
% Calling sequence
%
%   function [K,Delta]=gain_of_x(x,K0,x0,y0,delta,delta_change,rho)
%
%%
function [K,Delta]=gain_of_x(x,K0,x0,y0,delta,delta_change,rho)
%% Inputs
%
% * |x| point (n x nvec)
% * |K0| gain (nx1) (will be scaled by |1/delta|)
% * |x0| point (nx1) on periodic orbit where gain should start to be equal
% to |K0/delta|
% * |y0| approximate derivative of periodic orbit in |x0|
% * |delta| approximate length of time interval where gains are to be equal
% to |K0/delta|. This is approximate because time t is approximated by
% t(x)=y0^T(x-x0)/(y0^Ty0). )
% * |rho|: cutoff radius for gains (output K is zero if |x-x0|>rho) 
%% Outputs
% * |K| (n x nvec) true gains
% * |Delta| (n x nvec): |K| is |Delta*K0|.
%%
xshape=size(x);
[n,nvec]=size(x);
xr=reshape(x,[n,nvec]);
K0=K0(:);
on=ones(1,n);
onv=ones(1,nvec);
xdist=sqrt(sum((xr-x0(:,onv)).^2));
rg=xdist<rho;
yt=y0'/(y0'*y0);
t=yt*(xr-x0(:,onv));
Delta=zeros(1,nvec);
s_up=t<0&t>-delta_change;
Delta(s_up)=smooth_step((t(s_up)+delta_change)/delta_change)/delta;
s_on=t>=0&t<=delta;
Delta(s_on)=1/delta;
s_down=t>delta&t<delta+delta_change;
Delta(s_down)=smooth_step((delta+delta_change-t(s_down))/delta_change)/delta;
Delta(~rg)=0;
K=Delta(on,:).*K0(:,onv);
K=reshape(K,xshape);
end
%%
function y=smooth_step(x,delta)
if nargin<2
    delta=1;
end
nx=size(x);
x=x(:);
y=zeros(size(x));
xd=x/delta;
rg=xd>0&xd<1;
s=exp(-1./xd(rg));
s1=exp(-1./(1-xd(rg)));
y(rg)=s./(s+s1);
y(xd>=1)=1;
y=reshape(y,nx);
end
