function g=SpecExpAssign(P,b,lambda)
%% assign spectrum to P*expm(b*g) to lambda
L=SpecAssign(P,-P*b,lambda);
b_orth=null(b(:)');
log1pLb=log(1+L*b);
mat=[b,b_orth];
rhs=[log1pLb,log1pLb/(L*b)*(L*b_orth)];
g=(mat'\(rhs'))';
end
function gains=SpecAssign(A,b,lambda)
%% find gains such that A-b*gains has eigenvalues lambdas
n=length(b);
%%
S=b(:,ones(1,n));
for i=2:n
    S(:,i)=A*S(:,i-1);
end
en=[zeros(1,n-1),1];
lrow=en/S;
%% calcluate lrow*P(A)
for i=1:n
    lrow=lrow*(A-lambda(i)*eye(n));
end
gains=real(lrow);
end
