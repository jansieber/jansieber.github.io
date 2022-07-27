%% Solve nonlinear system with as many variables as equations
%%
function [x,converged,J]=ScSolve(f,x0,varargin)
%% optional arguments
defaults={'maxit',10,'tolerance',1e-12,'print',1,'damping',1,'df',[],'hjac',1e-4};
options=ScSetOptions(defaults,varargin,'pass_on');
if isempty(options.df)
    options.df=@(y)ScJacobian(f,y,options.hjac);
end
%% loop of Newton iteration
x=x0;
converged=0;
for i=1:options.maxit
    y=f(x);
    J=options.df(x);
    cor=-J\y;
    ynorm=norm(y,'inf');
    cornorm=norm(cor,'inf');
    x=x+options.damping*cor;
    if options.print>0
        fprintf('it=%d, |cor|=%g, |res|=%g\n',i,cornorm,ynorm);
    end
    if max(cornorm,ynorm)<options.tolerance
        converged=1;
        break
    end
end
end
