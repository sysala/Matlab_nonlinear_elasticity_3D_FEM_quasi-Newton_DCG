% =========================================================================
%
%  This program compares Newton and quasi-Newton methods for a static
%  elasto-plastic model with the von Mises yield criterion in 3D.
%  It is considered the strip-footing benchmark.
%  One can set optionally 2 types of finite elements,
%  different levels of mesh density and many other parameters and
%  inner solvers. More details can be found in the paper
%  J. Karatson, S. Sysala. M. Beres:  Quasi-Newton variable preconditioning
%  for nonlinear elasticity systems in 3D, 2024.
%
% ======================================================================
%

% add shared files between problems (solvers, iterative solvers)
addpath('../shared_files')

% if there is AGMG solver present in path (function agmg(...))
% "Notay, Y. AGMG software and documentation. Available at http://agmg.eu"
% Note: academic licence is free upon an email request
global AGMG_present
AGMG_present = 0;

%
% Mesh data
%

elem_type = 'P1'; % available choices of finite elements: 'P1', 'P2'
% for P2 elements, it is necessary to change

% the regularization parameter within inner solvers
eps_parameter = 0.1;
alpha = 0.5;

% a positive integer defining mesh density
% values 4, 8 and 12 are used in the paper
density = 4;

%
% Elastic material parameters
%
young = 1e7; % Young's modulus
poisson = 0.48; % Poisson's ratio

% other parameters of the von Mises model
% hardening parameter
gamma = 100; % yield stress
eps = eps_parameter * gamma; % regularization parameter

%
% constant surface force representing strip-footing
%
surface_force = [0, 0, -450];

%
% Mesh generation
%

% geometrical parameters (choose only integers)
size_xy = 10; % size of the body in direction x and y
size_z = 10; % size of the body in z-direction

% the mesh generation depending prescribed finite elements
switch (elem_type)
    case 'P1'
        [COORD, ELEM, SURF, NEUMANN, Q] = mesh_P1(density, size_xy, size_z);
        %         fprintf('P1 elements: \n')
    case 'P2'
        [COORD, ELEM, SURF, NEUMANN, Q] = mesh_P2(density, size_xy, size_z);
        fprintf('P2 elements: \n')
    otherwise
        disp('bad choice of element type');
end

%
% Data from the reference element
%

% quadrature points and weights for volume and surface integration
[Xi, WF] = quadrature_volume(elem_type);
[Xi_s, WF_s] = quadrature_surface(elem_type);

% local basis functions and their derivatives for volume and surface
[HatP, DHatP1, DHatP2, DHatP3] = local_basis_volume(elem_type, Xi);
[HatP_s, DHatP1_s, DHatP2_s] = local_basis_surface(elem_type, Xi_s);

%
% Number of nodes, elements and integration points + print
%
n_n = size(COORD, 2); % number of nodes
n_unknown = length(COORD(Q)); % number of unknowns
n_e = size(ELEM, 2); % number of elements
n_q = length(WF); % number of quadratic points
n_int = n_e * n_q; % total number of integrations points
%
%   fprintf('number of nodes =%d ',n_n);
%   fprintf('\n');
fprintf('number of unknowns =%d ', n_unknown);
fprintf('\n');
%   fprintf('number of elements =%d ',n_e);
%   fprintf('\n');
%   fprintf('number of integration points =%d ',n_int);
%   fprintf('\n');

%
% Load vector representing the strip-footing
%
n_e_s = size(NEUMANN, 2); % number of surface elements
n_q_s = length(WF_s); % number of quadrature points on a surface element
n_int_s = n_e_s * n_q_s; % number of integration points on the surface
% (on the upper side of the body)
f_t_int = surface_force' * ones(1, n_int_s); % size(f_t_int)=(3,n_int_s)
f = vector_traction(NEUMANN, COORD, f_t_int, HatP_s, DHatP1_s, DHatP2_s, WF_s);

%
% Values of material parameters at integration points
%
shear = young / (2 * (1 + poisson)); % shear modulus
bulk = young / (3 * (1 - 2 * poisson)); % bulk modulus
shear = shear * ones(1, n_int);
bulk = bulk * ones(1, n_int);

%
% assembling of auxiliary arrays for Newton's and quasi-Newton's methods
%
[B, K_V, K_elast, WEIGHT] = auxiliary_matrices(ELEM, COORD, bulk, shear, ...
    DHatP1, DHatP2, DHatP3, WF);

%
% Newton's and quasi-Newton's solvers
%

% initialization displacement
U_it = zeros(3, n_n);

% standard Newton method
tic;
[U_N, it_N, crit_hist_N] = newton(U_it, WEIGHT, K_V, B, f, Q, shear, bulk, alpha, gamma, eps);
time_Newton = toc;
fprintf("     solver's runtime:  " + time_Newton + "-----\n");

% quasi-Newton method - preconditioner 1 (variable elastic stiffness matrix)
tic;
[U_qN1, it_qN1, crit_hist_qN1, omega_hist_qN1] = newton_quasi1(U_it, WEIGHT, K_V, K_elast, B, f, Q, shear, bulk, alpha, gamma, eps);
time_qNewton1 = toc;
fprintf("     solver's runtime:  " + time_qNewton1 + "-----\n");

% quasi-Newton method - preconditioner 2 (fixed elastic stiffness matrix)
tic;
[U_qN2, it_qN2, crit_hist_qN2, omega_hist_qN2] = newton_quasi2(U_it, WEIGHT, K_elast, B, f, Q, shear, bulk, alpha, gamma, eps);
time_qNewton2 = toc;
fprintf("     solver's runtime:  " + time_qNewton2 + "-----\n");

%
% Visualization of selected results
%

% mesh
if density < 5
    draw_mesh(COORD, SURF, elem_type)
end

% total displacements + deformed shape - newton
U_total = sqrt(U_N(1, :) .^ 2 + U_N(2, :) .^ 2 + U_N(3, :) .^ 2);
draw_quantity(COORD, SURF, 10000 * U_N, U_total, elem_type, size_xy, size_z)

% changes of the shear modulus
E = reshape(B * U_N(:), 6, []);
IOTA = [1; 1; 1; 0; 0; 0];
VOL = IOTA * IOTA';
DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;
dev_E = DEV * E; % deviatoric part of E
z = max(0, sum(E .* dev_E)); % scalar product of the deviatoric strain
test = 2 * shear .* sqrt(z);
ind2 = (test > gamma - eps) & (test < gamma + eps);
ind3 = (test >= gamma + eps);
mu = shear;
mu(ind2) = (1 - alpha) * shear(ind2) + ...
    alpha * (gamma - (2 * shear(ind2) .* sqrt(z(ind2)) - gamma - eps) .^ 2 / (4 * eps)) ./ (2 * sqrt(z(ind2)));
mu(ind3) = (1 - alpha) * shear(ind3) + alpha * gamma ./ (2 * sqrt(z(ind3)));
mu_node = transformation(mu, ELEM, WEIGHT);
draw_quantity(COORD, SURF, 0 * U_N, mu_node, elem_type, size_xy, size_z);
colorbar off; colorbar('location', 'eastoutside')

% convergence of the Newton-like solvers
figure_convergence(1:it_N, crit_hist_N, ...
    1:it_qN1, crit_hist_qN1, ...
    1:it_qN2, crit_hist_qN2)

% line search coefficients
figure_omega(1:it_qN1 - 1, omega_hist_qN1, ...
    1:it_qN2 - 1, omega_hist_qN2)
