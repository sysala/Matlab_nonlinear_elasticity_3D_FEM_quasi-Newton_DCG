function [shur_prec_handle] = diag_prec(A, Q)
    % =========================================================================
    %
    % Fuction that costructs a block diagonal preconditioner.
    % Each block corresponds to a one direction in the vector field.
    % Depends on global variable AGMG_present if the inner solver is AGMG or
    % Incomplete Cholesky.
    %
    % Input data:
    %  A     - matrix
    %  Q     - vector of indices of boundary nodes
    %
    % Output data:
    %  shur_prec_handle - function handle to the preconditioner
    %
    % =========================================================================
    %

    Q1 = Q;
    Q1(2:3, :) = 0;
    Q1 = Q1(Q);

    Q2 = Q;
    Q2([1 3], :) = 0;
    Q2 = Q2(Q);

    Q3 = Q;
    Q3([1 2], :) = 0;
    Q3 = Q3(Q);

    A1 = A(Q1, Q1);
    A2 = A(Q2, Q2);
    A3 = A(Q3, Q3);

    global AGMG_present

    if AGMG_present
        agmg(A1, [], 1, [], [], -1, [], 1, 8, 4, 1);
        agmg(A2, [], 1, [], [], -1, [], 1, 8, 4, 2);
        agmg(A3, [], 1, [], [], -1, [], 1, 8, 4, 3);

        A1_solver = @(x)agmg(A1, x, 1, 1e-6, 10, -1, [], 3, 8, 4, 1);
        A2_solver = @(x)agmg(A2, x, 1, 1e-6, 10, -1, [], 3, 8, 4, 2);
        A3_solver = @(x)agmg(A3, x, 1, 1e-6, 10, -1, [], 3, 8, 4, 3);
    else
        params = struct('type', 'ict', 'droptol', 1e-3);
        L1 = ichol(A1, params);
        L2 = ichol(A2, params);
        L3 = ichol(A3, params);

        A1_solver = @(x)L1' \ (L1 \ x);
        A2_solver = @(x)L2' \ (L2 \ x);
        A3_solver = @(x)L3' \ (L3 \ x);
    end

    shur_prec_handle = @(x) prec_apply(x, Q1, Q2, Q3, A1_solver, A2_solver, A3_solver);

end

function [res] = prec_apply(x, Q1, Q2, Q3, A1, A2, A3)
    b1 = x(Q1);
    b2 = x(Q2);
    b3 = x(Q3);

    y1 = A1(b1);
    y2 = A2(b2);
    y3 = A3(b3);
    res = x * 0;

    res(Q1) = y1;
    res(Q2) = y2;
    res(Q3) = y3;
end
