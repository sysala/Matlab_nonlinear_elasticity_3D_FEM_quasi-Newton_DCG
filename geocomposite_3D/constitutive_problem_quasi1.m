function [S, DS, m, M] = constitutive_problem_quasi1(E, shear, bulk, alpha, gamma, eps)

    % =========================================================================
    %
    % The aim of this function is to construct constitutive operator and
    % auxiliary arrays for the quasi-Newton method 1 at integration points
    % 1,2,...,n_int.
    %
    % Input data:
    %  E       - current strain tensor, size(E)=(3,n_int)
    %  shear   - shear moduli at integration points, size(shear)=(1,n_int)
    %  bulk    - bulk moduli at integration points, size(bulk)=(1,n_int)
    %  alpha,gamma,eps - other material parameters
    %
    % Output data:
    %  S      - stress tensors at integration points, size(S)=(4,n_int)
    %  DS     - auxiliary array for the quasi-Newton stiffness matrix.
    %           (It represent the function \mu_m.) size(DS)=(1,n_int)
    %  M,m    - parameters defining the quasi-Newton step
    %
    % =========================================================================
    %

    %
    % Deviatoric and volumetric 6x6 matrices
    %
    IOTA = [1; 1; 1; 0; 0; 0];
    VOL = IOTA * IOTA';
    DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;

    %
    % Deviatoric strain and its norm
    %
    dev_E = DEV * E; % deviatoric part of E
    z = max(0, sum(E .* dev_E)); % scalar product of the deviatoric strain

    %
    % Branching on the elastic, plastic and smoothing parts
    %
    test = 2 * shear .* sqrt(z);
    ind2 = (test > gamma - eps) & (test < gamma + eps);
    ind3 = (test >= gamma + eps);

    %
    % The array mu
    %
    mu = shear;
    mu(ind2) = (1 - alpha(ind2)) .* shear(ind2) + ...
        alpha(ind2) .* (gamma(ind2) - (2 * shear(ind2) .* sqrt(z(ind2)) - gamma(ind2) - eps(ind2)) .^ 2 ./ (4 * eps(ind2))) ./ (2 * sqrt(z(ind2)));
    mu(ind3) = (1 - alpha(ind3)) .* shear(ind3) + alpha(ind3) .* gamma(ind3) ./ (2 * sqrt(z(ind3)));

    %
    % The array Dmu
    %
    Dmu = zeros(1, length(shear));
    Dmu(ind2) = -alpha(ind2) .* (gamma(ind2) - (2 * shear(ind2) .* sqrt(z(ind2)) - gamma(ind2) - eps(ind2)) .^ 2 ./ (4 * eps(ind2))) ./ (4 * sqrt(z(ind2))) + ...
        -alpha(ind2) .* shear(ind2) .* (2 * shear(ind2) .* sqrt(z(ind2)) - gamma(ind2) - eps(ind2)) ./ (4 * eps(ind2));
    Dmu(ind3) = -alpha(ind3) .* gamma(ind3) ./ (4 * sqrt(z(ind3)));

    %
    % The stress tensor
    %
    S = repmat(bulk, 6, 1) .* (VOL * E) + 2 * repmat(mu, 6, 1) .* dev_E;

    %
    % Auxiliary array for the quasi-Newton stiffness matrix
    %
    delta = 0.5; % a parameter defining the quasi-Newton matrix
    mu_delta = mu + 2 * delta * Dmu;
    DS = 2 * mu_delta;

    %
    % Parameters m and M defining the quasi-Newton step
    %
    mu_1 = mu + 2 * Dmu;
    m = min(mu_1 ./ mu_delta);
    M = max(mu ./ mu_delta);

end
