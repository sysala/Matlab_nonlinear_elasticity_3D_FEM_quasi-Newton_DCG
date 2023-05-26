function [y] = prec_blkdiag(x,A_N1,A_N2,A_N3,idx1,idx2,idx3)
%PREC_BLKDIAG Summary of this function goes here
%   Detailed explanation goes here
y = x*0;

global P2 P3 P1 deflation_all

%% inner DCG solver
if isempty(deflation_all)
    deflation_all = cell(3,1);
    deflation_all{1} = [];
    deflation_all{2} = [];
    deflation_all{3} = [];
end

if isempty(P1)
    regularizacni_parametr = 1e-4;
    ichol_params = struct('type', 'ict', 'droptol', 4e-3);
    
    L1 = ichol(A_N1 + regularizacni_parametr*diag(diag(A_N1)), ichol_params);
    P1 = @(x)L1'\(L1\x);
    
    L2 = ichol(A_N2 + regularizacni_parametr*diag(diag(A_N2)), ichol_params);
    P2 = @(x)L2'\(L2\x);
    
    L3 = ichol(A_N3 + regularizacni_parametr*diag(diag(A_N3)), ichol_params);
    P3 = @(x)L3'\(L3\x);
end

tol = 2e-1;
max_iter = 100;

[y(idx1), iter1, resvec, flag] = ...
    DCG( A_N1, x(idx1), x(idx1)*0,deflation_all{1},P1,tol,max_iter);
if iter1 > 3
    [ W_orth] = my_orth_simple( deflation_all{1}, y(idx1));
    deflation_all{1} = [deflation_all{1} W_orth];
end

[y(idx2), iter2, resvec, flag] = ...
    DCG( A_N2, x(idx2), x(idx2)*0,deflation_all{2},P2,tol,max_iter);
if iter2 > 3
    [ W_orth] = my_orth_simple( deflation_all{2}, y(idx2));
    deflation_all{2} = [deflation_all{2} W_orth];
end

[y(idx3), iter3, resvec, flag] = ...
    DCG( A_N3, x(idx3), x(idx3)*0,deflation_all{3},P3,tol,max_iter);
if iter3 > 3
    [ W_orth] = my_orth_simple( deflation_all{3}, y(idx3));
    deflation_all{3} = [deflation_all{3} W_orth];
end
%fprintf("(%d,%d,%d)",iter1,iter2,iter3);

%% direct solver
% if isempty(P1)
%    
% 
%     P1 = @(x)A_N1\x;
%     P2 = @(x)A_N2\x;
%     P3 = @(x)A_N3\x;
%     
% end
% 
% y(idx1)=P1(x(idx1));
% y(idx2)=P2(x(idx2));
% y(idx3)=P3(x(idx3));

end

