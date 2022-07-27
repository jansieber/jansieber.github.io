%% Wrapper around user-defined function to be used in DDE-Biftool
%
%
%   function y=etdf_controlled_rhs(xx,p,oderhs,varargin)
%
%
% If |varargin| is empty then |y=oderhs(xx,p,0)| where |xx| was reshaped:
% |xx=reshape(xx,[n,nvec])|. Thus, oderhs can be used to perform
% bifurcation analysis on the ODE |x'=oderhs(x,p,0)|.
%
% If options arguments |'epsilon'| and |'K'| are given then |y=[y1;y2]|
% where |y1=f(x(t),p,K(x(t))[xtilde(t)-x(t)])| and
% |y2=xtilde(t)-(1-epsilon)xtilde(t-T)-epsilon x(t-T)|
%
%%
function y=etdf_controlled_rhs(xx,p,oderhs,varargin)
default={'K',[],'epsilon',[]};
options=dde_set_options(default,varargin);
n=size(xx,1);
nd=size(xx,2);
nvec=size(xx,3);
if isempty(options.K)
    x=reshape(xx,[n,nvec]);
    y=oderhs(x,p,zeros(1,nvec));
else
    assert(nd==2);
    nx=n/2;
    assert(mod(nx,1)==0);
    x=reshape(xx(1:nx,1,:),[nx,nvec]);
    x_tau=reshape(xx(1:nx,2,:),[nx,nvec]);
    xtilde=reshape(xx(nx+1:2*nx,1,:),[nx,nvec]);
    xtilde_tau=reshape(xx(nx+1:2*nx,2,:),[nx,nvec]);
    K=options.K(x);
    u=sum(K.*(xtilde-x));
    y_ode=oderhs(x,p,u);
    y_recursion=xtilde-((1-options.epsilon)*xtilde_tau+options.epsilon*x_tau);
    y=cat(1,y_ode,y_recursion);
end
end
