function [S, DS] = constitutive_problem(E, shear, bulk, alpha, gamma, eps)

    % =========================================================================
    %
    % The aim of this function is to construct constitutive and consistent
    % tangent operators at integration points 1,2,...,n_int for the semismooth
    % Newton method.
    %
    % Input data:
    %  E       - current strain tensor, size(E)=(6,n_int)
    %  shear   - shear moduli at integration points, size(shear)=(1,n_int)
    %  bulk    - bulk moduli at integration points, size(bulk)=(1,n_int)
    %  alpha,gamma,eps - other material parameters
    %
    % Output data:
    %  S      - stress tensors at integration points, size(S)=(6,n_int)
    %  DS     - consistent tangent operators at integr. points,
    %           size(DS)=(36,n_int)
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
    s = sqrt(z); % norm of the deviatoric strain

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
    mu(ind2) = (1 - alpha) * shear(ind2) + ...
        alpha * (gamma - (2 * shear(ind2) .* sqrt(z(ind2)) - gamma - eps) .^ 2 / (4 * eps)) ./ (2 * sqrt(z(ind2)));
    mu(ind3) = (1 - alpha) * shear(ind3) + alpha * gamma ./ (2 * sqrt(z(ind3)));

    %
    % The array Dmu
    %
    Dmu = zeros(1, length(shear));
    Dmu(ind2) = -alpha * (gamma - (2 * shear(ind2) .* sqrt(z(ind2)) - gamma - eps) .^ 2 / (4 * eps)) ./ (4 * sqrt(z(ind2))) + ...
        -alpha * shear(ind2) .* (2 * shear(ind2) .* sqrt(z(ind2)) - gamma - eps) / (4 * eps);
    Dmu(ind3) = -alpha * gamma ./ (4 * sqrt(z(ind3)));

    %
    % The stress tensor
    %

    S = repmat(bulk, 6, 1) .* (VOL * E) + 2 * repmat(mu, 6, 1) .* dev_E;

    %
    % The tangent stiffness matrix
    %

    DS = 2 * DEV(:) * mu;
    IND = (s > 0);
    N_hat = dev_E(:, IND) ./ repmat(s(IND), 6, 1);
    NN_hat = repmat(N_hat, 6, 1) .* kron(N_hat, ones(6, 1));
    DS(:, IND) = DS(:, IND) + 4 * repmat(Dmu(IND), 36, 1) .* NN_hat;

end
