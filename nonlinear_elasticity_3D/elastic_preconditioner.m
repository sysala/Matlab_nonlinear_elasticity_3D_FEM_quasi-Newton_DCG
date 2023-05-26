function [y] = elastic_preconditioner(x)
%ELASTIC_PRECONDITIONER Summary of this function goes here
%   Detailed explanation goes here


global krylov_recycle A_elast L_elast

if isempty(krylov_recycle)
    regularizacni_parametr = 1e-4;
    %ichol_params = struct('michol','on');%
    ichol_params = struct('type', 'ict', 'droptol', 1e-3);
    L_elast = ichol(A_elast + regularizacni_parametr*diag(diag(A_elast)), ichol_params);
    krylov_recycle = struct();
    krylov_recycle.max_size = 200;
    krylov_recycle.max_iter = 0;
    krylov_recycle.WTAW_inv_Wt_A=[];
    krylov_recycle.krylov_space_vals_all =[];
    krylov_recycle.W = [];
end

tol = 1e-2;
max_iter = 100;

[y, iter, resvec, flag,krylov_recycle] = ...
    DCG_adaptive( A_elast, x, x*0,@(x)L_elast'\(L_elast\x),tol,max_iter,krylov_recycle);

% fprintf("%d \n",iter);

end

