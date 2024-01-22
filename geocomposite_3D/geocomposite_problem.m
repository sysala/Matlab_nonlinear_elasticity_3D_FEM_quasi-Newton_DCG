% =========================================================================
%
%  This program compares Newton and quasi-Newton methods for a static
%  elasto-plastic model with the von Mises yield criterion in 3D.
%  It is considered the geocomposite benchmark with two different sets
%  of data. One can set optionally 2 types of finite elements,
%  and many other parameters and inner solvers. For more details, we refer
%  to J. Karatson, S. Sysala. M. Beres:  Quasi-Newton variable
%  preconditioning for nonlinear elasticity systems in 3D, 2024.
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
% Parameters specifying the first and second materials
%
% first material (softer)
young1 = 300; % Young's modulus
poisson1 = 0.4; % Poisson's ratio
gamma1 = 10; % initial yield stress
alpha1 = 0.5; % hardening parameter
eps1 = 0.01 * gamma1; % regularization parameter
%
shear1 = young1 / (2 * (1 + poisson1)); % shear modulus
bulk1 = young1 / (3 * (1 - 2 * poisson1)); % bulk modulus

% second material (harder)
young2 = 2000; % Young's modulus
poisson2 = 0.3; % Poisson's ratio
gamma2 = 30; % initial yield stress
alpha2 = 0.5; % hardening parameter
eps2 = 0.01 * gamma2; % regularization parameter
%
shear2 = young2 / (2 * (1 + poisson2)); % shear modulus
bulk2 = young2 / (3 * (1 - 2 * poisson2)); % bulk modulus

%
% Constant uniaxial pressure (applied in y-direction)
%
surface_force = [0, 0, -20];

%
% Geometry, mesh data and creation of the mesh
%

% choice of finite elements
elem_type = 'P1'; % available choices of finite elements: 'P1', 'P2'
% for P2 elements, it is necessary to change
% the regularization parameter within inner solvers

% size of the investigated cube
size_xy = 13; % size of the body in direction x and y
size_z = 13; % size of the body in z-direction

% heterogeneity - there are two available input files:
%         either 'heter_mesh1.mat' or 'heter_mesh2.mat'.
C_heter = load('heter_mesh1.mat');
threshold = 135; % threshold parameter defining the heterogeneity
%
heterogeneity = repmat(C_heter.a(:)', 6, 1);
heter = heterogeneity(:);
Q_elem = heter <= threshold; % heterogeneity in each element
%
C_heter1 = C_heter.a(:, :, 1);
heterogeneity1 = repmat(C_heter1(:)', 2, 1);
heter1 = heterogeneity1(:);
Q_surf1 = heter1 <= threshold; % heterogeneity on the bottom face
%
C_heter2 = C_heter.a(:, :, end);
heterogeneity2 = repmat(C_heter2(:)', 2, 1);
heter2 = heterogeneity2(:);
Q_surf2 = heter2 <= threshold; % heterogeneity on the top face
%
C_heter3 = C_heter.a(:, 1, :);
heterogeneity3 = repmat(C_heter3(:)', 2, 1);
heter3 = heterogeneity3(:);
Q_surf3 = heter3 <= threshold; % heterogeneity on the front face
%
C_heter4 = C_heter.a(end, :, :);
heterogeneity4 = repmat(C_heter4(:)', 2, 1);
heter4 = heterogeneity4(:);
Q_surf4 = heter4 <= threshold; % heterogeneity on the right face
%
C_heter5 = C_heter.a(:, end, :);
heterogeneity5 = repmat(C_heter5(:)', 2, 1);
heter5 = heterogeneity5(:);
Q_surf5 = heter5 <= threshold; % heterogeneity on the back face
%
C_heter6 = C_heter.a(1, :, :);
heterogeneity6 = repmat(C_heter6(:)', 2, 1);
heter6 = heterogeneity6(:);
Q_surf6 = heter6 <= threshold; % heterogeneity on the left face
%
Q_surf = [Q_surf1; Q_surf2; Q_surf3; Q_surf4; Q_surf5; Q_surf6];

% mesh density in each direction
N_x = size(C_heter.a, 1);
N_z = size(C_heter.a, 3);

% the mesh generation depending prescribed finite elements
switch (elem_type)
    case 'P1'
        [COORD, ELEM, SURF, NEUMANN, Q] = mesh_P1(size_xy, size_z, N_x, N_z);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD, ELEM, SURF, NEUMANN, Q] = mesh_P2(size_xy, size_z, N_x, N_z);
        fprintf('P2 elements: \n')
    otherwise
        disp('bad choice of element type');
end

draw_mesh(COORD, SURF, Q_surf)

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
fprintf('number of nodes =%d ', n_n);
fprintf('\n');
fprintf('number of unknowns =%d ', n_unknown);
fprintf('\n');
fprintf('number of elements =%d ', n_e);
fprintf('\n');
fprintf('number of integration points =%d ', n_int);
fprintf('\n');

%
% Material parameters at integration points w.r.t. the heterogeneity
%
shear = zeros(1, n_int);
shear(Q_elem) = shear1;
shear(~Q_elem) = shear2;
%
bulk = zeros(1, n_int);
bulk(Q_elem) = bulk1;
bulk(~Q_elem) = bulk2;
%
gamma = zeros(1, n_int);
gamma(Q_elem) = gamma1;
gamma(~Q_elem) = gamma2;
%
alpha = zeros(1, n_int);
alpha(Q_elem) = alpha1;
alpha(~Q_elem) = alpha2;
%
eps = zeros(1, n_int);
eps(Q_elem) = eps1;
eps(~Q_elem) = eps2;

%
% Load vector representing the uniaxial compression
%
n_e_s = size(NEUMANN, 2); % number of surface elements
n_q_s = length(WF_s); % number of quadrature points on a surface element
n_int_s = n_e_s * n_q_s; % number of integration points on the surface
% (on the upper side of the body)
f_t_int = surface_force' * ones(1, n_int_s); % size(f_t_int)=(3,n_int_s)
f = vector_traction(NEUMANN, COORD, f_t_int, HatP_s, DHatP1_s, DHatP2_s, WF_s);

%
% Assembling of auxiliary arrays for Newton's and quasi-Newton's methods
%
[B, K_V, K_elast, WEIGHT] = auxiliary_matrices(ELEM, COORD, bulk, shear, ...
    DHatP1, DHatP2, DHatP3, WF);

%
% Newton's and quasi-Newton's solvers
%

% initialization displacement
U_it = zeros(3, n_n);
% standard Newton method, dcg + incomplete Cholesky
tic;
[U_N, it_N, crit_hist_N] = newton(U_it, WEIGHT, K_V, B, f, Q, shear, bulk, alpha, gamma, eps);
time_Newton = toc;
fprintf("     solver's runtime:  " + time_Newton + "-----\n");

% quasi-Newton method - preconditioner 1, dcg + incomplete Cholesky
tic;
[U_qN1, it_qN1, crit_hist_qN1, omega_hist_qN1] = newton_quasi1(U_it, WEIGHT, K_V, K_elast, B, f, Q, shear, bulk, alpha, gamma, eps);
time_qNewton1 = toc;
fprintf("     solver's runtime:  " + time_qNewton1 + "-----\n");

% quasi-Newton method - preconditioner 2, dcg + block-diagonal precond.
tic;
[U_qN2, it_qN2, crit_hist_qN2, omega_hist_qN2] = newton_quasi2(U_it, WEIGHT, K_elast, B, f, Q, shear, bulk, alpha, gamma, eps);
time_qNewton2 = toc;
fprintf("     solver's runtime:  " + time_qNewton2 + "-----\n");

%
% Visualization of selected results
%

% mesh
draw_mesh(COORD, SURF, Q_surf)

% total displacements + deformed shape - newton
U_total = sqrt(U_N(1, :) .^ 2 + U_N(2, :) .^ 2 + U_N(3, :) .^ 2);
draw_quantity(COORD, SURF, 10 * U_N, U_total, elem_type, size_xy, size_z)

% changes of the shear modulus
E = reshape(B * U_N(:), 6, []);
IOTA = [1; 1; 1; 0; 0; 0];
VOL = IOTA * IOTA';
DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;
dev_E = DEV * E; % deviatoric part of E
z = max(0, sum(E .* dev_E)); % scalar product of the deviatoric strain
test = 2 * shear .* sqrt(z);
ind1 = (test <= gamma - eps);
ind2 = (test > gamma - eps) & (test < gamma + eps);
ind3 = (test >= gamma + eps);
mu = shear;
mu(ind2) = (1 - alpha(ind2)) .* shear(ind2) + ...
    alpha(ind2) .* (gamma(ind2) - (2 * shear(ind2) .* sqrt(z(ind2)) - gamma(ind2) - eps(ind2)) .^ 2 ./ (4 * eps(ind2))) ./ (2 * sqrt(z(ind2)));
mu(ind3) = (1 - alpha(ind3)) .* shear(ind3) + alpha(ind3) .* gamma(ind3) ./ (2 * sqrt(z(ind3)));
mu_node = transformation(mu ./ shear, ELEM, WEIGHT);
draw_quantity(COORD, SURF, 0 * U_N, mu_node, elem_type, size_xy, size_z);
colorbar off; colorbar('location', 'eastoutside')

% convergence of the Newton-like solvers
figure_convergence(1:it_N, crit_hist_N, ...
    1:it_qN1, crit_hist_qN1, ...
    1:it_qN2, crit_hist_qN2)

% line search coefficients
figure_omega(1:it_qN1 - 1, omega_hist_qN1, ...
    1:it_qN2 - 1, omega_hist_qN2)
