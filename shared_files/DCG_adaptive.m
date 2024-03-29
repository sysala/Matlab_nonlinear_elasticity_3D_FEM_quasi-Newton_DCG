function [x, iter, resvec, tag, krylov_recycle] = ...
        DCG_adaptive(A, b, x0, M, tol, maxiter, krylov_recycle)
    % =========================================================================
    %
    %  Conjugate gradient method with deflation using Krylov subspace recycling
    %  i.e. deflation space is constructed from the previous Krylov subspaces
    %  from solution with the same matrix A
    %
    %  input data:
    %    A     - the matrix of the system, size(A)=(n,n)
    %    b     - the right hand side vector, size(b)=(n,1)
    %    x0    - the initial guess for the solution, size(x0)=(n,1)
    %    M     - preconditioner (operator)
    %    tol   - the tolerance parameter
    %    maxiter - the maximum number of iterations
    %    krylov_recycle - structure with the following fields:
    %      max_size - maximum size of the deflation space
    %      max_iter - maximum number of iterations
    %      W        - the deflation space
    %      WTAW_inv_Wt_A - the matrix (W'*A*W)^(-1)*W'*A
    %      krylov_space_vals_all - the values of the matrix W'*A*W
    %
    %  output data:
    %    x     - the computed solution
    %    iter  - the number of iterations performed
    %    resvec - the history of residual norms
    %    tag   - the convergence tag
    %    krylov_recycle - updated structure
    %
    % ======================================================================
    %

    if isempty(krylov_recycle.W)
        P = @(x)x;
        Q = @(x)0;
    else
        WTAW_inv = diag(1 ./ krylov_recycle.krylov_space_vals_all);
        Q = @(x)krylov_recycle.W * (WTAW_inv * (x' * krylov_recycle.W)');
        P = @(x)x - krylov_recycle.W * (krylov_recycle.WTAW_inv_Wt_A * x);
    end

    if isempty(M)
        M = @(x)x;
    end

    if isa(M, 'numeric')
        M_mat = M;
        M = @(x)M_mat \ x;
    end

    if isempty(x0)
        x0 = 0 * b;
    end

    A = @(x)A * x;

    r0 = b - A(x0);
    x = x0 + Q(r0);
    b_norm = norm(b);
    res = norm(A(x) - b) / b_norm;

    if res < tol || maxiter == 0 || b_norm == 0
        tag = 0;
        resvec = res;
        iter = 0;
        return
    end

    r = b - A(x);
    z = M(r);
    p = P(z);

    gamma_old = dot(r, z);
    tag = 3;
    resvec = zeros(maxiter + 1, 1);
    resvec(1) = res;
    W = {};
    AW = {};
    krylov_space_vals_all = {};

    for j = 1:maxiter

        for jj = 2:j
            kk = jj - 2;
            beta2 =- dot(p, AW{end - kk}) / krylov_space_vals_all{end - kk};
            p = p + beta2 * W{end - kk};
        end

        s = A(p);
        %s=P(A(p));
        tmp_dot = dot(s, p);
        alpha = gamma_old / tmp_dot;
        x = x + alpha * p;

        W{end + 1} = p;
        AW{end + 1} = s;
        krylov_space_vals_all{end + 1} = tmp_dot;

        r = r - alpha * s;
        res = norm(r) / b_norm;
        resvec(j + 1) = res;

        if res < tol
            tag = 1;
            break;
        end

        z = M(r);
        z = P(z);

        gamma_new = dot(r, z);
        beta = gamma_new / gamma_old;

        p = z + beta * p;

        gamma_old = gamma_new;
    end

    krylov_recycle.max_iter = max(krylov_recycle.max_iter, j);

    if size(krylov_recycle.W, 2) < krylov_recycle.max_size && j > ceil(sqrt(krylov_recycle.max_iter))
        %     fprintf("Expand space.  ")
        W = cell2mat(W);
        AW = cell2mat(AW);
        krylov_space_vals_all = cell2mat(krylov_space_vals_all);

        krylov_recycle.W = [krylov_recycle.W W];
        %krylov_recycle.AW{end+1} = AW;
        krylov_recycle.WTAW_inv_Wt_A = [krylov_recycle.WTAW_inv_Wt_A; diag(1 ./ krylov_space_vals_all) * AW'];
        krylov_recycle.krylov_space_vals_all = [krylov_recycle.krylov_space_vals_all krylov_space_vals_all];
    end

    resvec = resvec(1:j + 1);
    iter = j;
end
