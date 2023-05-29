function [ x,iter,resvec,tag] = DCG( A,b,x0,W,M,tol,maxiter)
%   DPCG Summary of this function goes here
%   Detailed explanation goes her

if isempty(W)
    P=@(x)x;
    Q=@(x)0;  
else
    WTAW=W'*A*W;
    WTAW_inv = inv(WTAW);
    WTAW_inv_Wt_A = WTAW_inv*W'*A;
    Q=@(x)W*(WTAW_inv*((W')*x));  
    P=@(x)x-W*(WTAW_inv_Wt_A*x);
end

if isempty(M)
    M=@(x)x;
end
if isa(M, 'numeric')
    M_mat=M;
    M=@(x)M_mat\x;
end
if isempty(x0)
    x0=0*b;
end
A = @(x)A*x;

r0=b-A(x0);
x=x0+Q(r0);
b_norm=norm(b);
res=norm(A(x)-b)/b_norm;
if res<tol || maxiter==0 || b_norm==0
    tag=0;
    resvec=res;
    iter=0;
    return
end
r=b-A(x);
z=M(r);
p=P(z);

gamma_old=dot(r,z);
tag=3;
resvec=zeros(maxiter+1,1);
resvec(1)=res;
for j=1:maxiter
    s=A(p);
    %s=P(A(p));
    alpha=gamma_old/dot(s,p);
    x=x+alpha*p;
    %res=norm(A(x)-b)/b_norm;
    r=r-alpha*s;
    res = norm(r)/b_norm;
    resvec(j+1)=res;
    if res<tol
        tag=1;
        break;
    end
    z=M(r);
    gamma_new=dot(r,z);
    beta=gamma_new/gamma_old;
    p=P(z)+beta*p;
    %p=z+beta*p;
    gamma_old=gamma_new;
end

resvec=resvec(1:j+1);
iter=j;
end