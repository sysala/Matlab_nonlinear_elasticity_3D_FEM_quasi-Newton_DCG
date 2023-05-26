function [ x,iter,resvec,tag, krylov_recycle] = ...
    DCG_adaptive( A,b,x0,M,tol,maxiter,krylov_recycle)
%DPCG Summary of this function goes here
%   Detailed explanation goes her
if isempty(krylov_recycle.W)
    P=@(x)x;
    Q=@(x)0;
else
    WTAW_inv=diag(1./krylov_recycle.krylov_space_vals_all);
    Q=@(x)krylov_recycle.W*(WTAW_inv*(x'*krylov_recycle.W)');
    P=@(x)x-krylov_recycle.W*(krylov_recycle.WTAW_inv_Wt_A*x);
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
W = {};
AW = {};
krylov_space_vals_all = {};
for j=1:maxiter
    s=A(p);
    %s=P(A(p));
    tmp_dot = dot(s,p);
    alpha=gamma_old/tmp_dot;
    x=x+alpha*p;
    
    W{end+1} = p;
    AW{end+1} = s;
    krylov_space_vals_all{end+1} = tmp_dot;
    
    r=r-alpha*s;
    res = norm(r)/b_norm;
    resvec(j+1)=res;
    if res<tol
        tag=1;
        break;
    end
    z=M(r);
    z=P(z);
    
    for jj=2:j
        kk = jj-1;
        beta2 = - dot(z,AW{end-kk})/krylov_space_vals_all{end-kk};
        z=z+beta2*W{end-kk};
    end
    
    gamma_new=dot(r,z);
    beta=gamma_new/gamma_old;
    
    p=z+beta*p;
    
    gamma_old=gamma_new;
end
krylov_recycle.max_iter = max(krylov_recycle.max_iter,j);
if size(krylov_recycle.W,2)<krylov_recycle.max_size && j>ceil(sqrt(krylov_recycle.max_iter))
%     fprintf("Expand space.  ")
    W = cell2mat(W);
    AW = cell2mat(AW);
    krylov_space_vals_all=cell2mat(krylov_space_vals_all);
    
    krylov_recycle.W = [krylov_recycle.W W];
    %krylov_recycle.AW{end+1} = AW;
    krylov_recycle.WTAW_inv_Wt_A=[krylov_recycle.WTAW_inv_Wt_A;diag(1./krylov_space_vals_all)*AW'];
    krylov_recycle.krylov_space_vals_all = [krylov_recycle.krylov_space_vals_all krylov_space_vals_all];
end
resvec=resvec(1:j+1);
iter=j;
end
