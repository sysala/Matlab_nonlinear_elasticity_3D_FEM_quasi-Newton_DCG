function [U, it, crit_hist] = newton(U_ini, WEIGHT, K_V, B, f, Q, shear, bulk, alpha, gamma, eps)

    %--------------------------------------------------------------------------
    % The Newton method for solution of the system
    %              find U:   F(U)=f
    % Linearized systems are solved by DCG with incomplete Cholesky
    % preconditioner.
    %
    % Input data:
    %   U_ini - initial choice of U
    %   WEIGHT - weight coefficients of integration points, size(WEIGHT)=(1,n_int)
    %   K_V    - volumetric part of the elastic stiffness matrix, size(K_V)=(3*n_n,3*n_n)
    %   B      - the strain-displacement matrix, size(B)=(6*n_int,3*n_n)
    %   f      - vector of external forces, size(f)=(3,n_n)
    %   Q      - logical array indicating the nodes where the homogeneous
    %            Dirichlet boundary condition is considered, size(Q)=(3,n_n)
    %   shear, bulk - elastic material parameters defined for each integration
    %                 point, size(shear)=size(bulk)=(1,n_int)
    %   alpha,gamma,eps - other material parameters
    %
    % Output data:
    %   U - approximation of the solution, size(U)=(3,n_n)
    %   it - number of Newton's iteration
    %   crit_hist - evolution of the stopping criterion
    %
    %--------------------------------------------------------------------------

    %
    % Auxiliary arrays and initialization
    %
    n_int = length(WEIGHT); % number of integration points
    n_n = size(U_ini, 2); % number of nodes
    dU = zeros(3, n_n); % Newton's increment (in displacement)
    F = zeros(3, n_n); % vector of internal forces
    E = zeros(6, n_int); % values of the strain tensor at integration points
    U = U_ini; % initial approximation of the solution

    % auxiliary arrays for the sparse stiffness matrix assembly
    AUX = reshape(1:6 * n_int, 6, n_int);
    iD = repmat(AUX, 6, 1);
    jD = kron(AUX, ones(6, 1));

    %
    it_max = 100; % maximal number of iterations
    crit_hist = zeros(1, it_max);
    deflation_basis = [];

    %
    % Newton's solver (the semismooth Newton method)
    %

    it = 0; % iteration number
    itcg = 0;

    while true

        it = it + 1;

        % constitutive operator and its derivative
        E(:) = B * U(:); % strain at integration points
        % solution of the constitutive problem
        [S, DS] = constitutive_problem(E, shear, bulk, alpha, gamma, eps);

        % vector of internal forces
        F(:) = B' * reshape(repmat(WEIGHT, 6, 1) .* S, 6 * n_int, 1);
        b_N = f(Q) - F(Q);

        % tangential stiffness matrix
        vD = repmat(WEIGHT, 36, 1) .* DS;
        D_p = sparse(iD(:), jD(:), vD(:), 6 * n_int, 6 * n_int);
        K_tangent = K_V + B' * D_p * B;
        A_N = K_tangent(Q, Q);
        A_N = (A_N + A_N') / 2;
        precond = diag_prec(A_N, Q);

        % stopping criterion
        criterion = norm(b_N) / norm(f(Q));
        crit_hist(it) = criterion;

        if criterion < 1e-10
            crit_hist = crit_hist(1:it);
            fprintf(' Standard Newton method converges ');
            fprintf(' number of iteration=%d  ', it);
            fprintf(' cumulative cg iterations=%d  ', itcg);
            fprintf(' stopping criterion=%e  ', criterion);
            fprintf('\n');
            break
        end

        % test on number of iteration
        if it == it_max
            fprintf('     Standard Newton method converges slowly: stopping criterion=%e  ', criterion)
            fprintf('\n');
            break
        end

        % deflated pcg solver with incomplete cholesky preconditioner
        tol = 1e-4; % precision of DCG solver
        max_iter = 1000; % DCG max iters
        [dU(Q), iter, ~, ~] = DCG(A_N, b_N, b_N * 0, deflation_basis, precond, tol, max_iter);
        itcg = itcg + iter;
        [W_orth] = my_orth_simple(deflation_basis, dU(Q)); % orthogonalization of deflation space
        deflation_basis = [deflation_basis W_orth]; % expanding deflation space
        % fprintf('%d ', iter);

        % update of the unknown vector
        U = U + dU;

    end % true

end % function
