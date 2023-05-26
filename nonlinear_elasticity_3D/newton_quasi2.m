function [U, it, crit_hist, omega_hist]=newton_quasi2(U_ini,WEIGHT,K_elast,B,f,Q,shear,bulk,model_type)

%--------------------------------------------------------------------------
% The quasi-Newton method with Preconditioner 5 for solution of the system
%              find U:   F(U)=f
% Linearized systems are solved by DCG with incomplete Cholesky
% preconditioner. It is used the fact that the stiffness matrix is
% constant.
%
% Input data:
%   U_ini - initial choice of U
%   WEIGHT- weight coefficients of integration points, size(WEIGHT)=(1,n_int)
%   K_elast - the elastic stiffness matrix, size(K_V)=(3*n_n,3*n_n)
%   B     - the strain-displacement matrix, size(B)=(6*n_int,3*n_n)
%   f     - vector of external forces, size(f)=(3,n_n)
%   Q     - logical array indicating the nodes where the homogeneous
%           Dirichlet boundary condition is considered, size(Q)=(3,n_n)
%   shear, bulk - elastic material parameters defined for each integration
%                 point, size(shear)=size(bulk)=(1,n_int)
%   model_type  - type of a nonlinear model
%
% Output data:
%   U - approximation of the solution, size(U)=(3,n_n)
%   it - number of Newton's iteration
%   crit_hist  - evolution of the stopping criterion
%   omega_hist - evolution of the damped coefficients
%
%--------------------------------------------------------------------------

%
% Auxiliary arrays and initialization
%

  n_int=length(WEIGHT);% number of integration points
  n_n=size(U_ini,2);   % number of nodes
  dU = zeros(3,n_n) ;  % Newton's increment (in displacement)
  F = zeros(3,n_n) ;   % vector of internal forces
  E = zeros(6,n_int);  % values of the strain tensor at integration points
  U=U_ini;             % initial approximation of the solution

  % elastic stiffness matrix (related to Precondtioner 2)
  A_N = K_elast(Q,Q);  
  A_N = (A_N + A_N')/2;

  % incomplete Cholesky for the elastic stiffness matrix
  reg_parameter = 1e-4;
  ichol_params = struct('type', 'ict', 'droptol', 1e-3);
  L = ichol(A_N + reg_parameter*diag(diag(A_N)), ichol_params);

  % arrays for recycling whole Krylov spaces
  krylov_recycle = struct();
  krylov_recycle.max_size = 1000;
  krylov_recycle.max_iter = 0;
  krylov_recycle.WTAW_inv_Wt_A=[];
  krylov_recycle.krylov_space_vals_all =[];
  krylov_recycle.W = [];
  
  %
  it_max=100;
  crit_hist=zeros(1,it_max);
  omega_hist=zeros(1,it_max);

%  
% Quasi-Newton's solver (Preconditioner 2)
%

  it=0;               % iteration number

  while true
    
    it=it+1;
    
    % solution of the constitutive problem
    E(:) = B*U(:) ;   % strain at integration points
    [S,m,M]=constitutive_problem_quasi2(E,shear,bulk,model_type);
                      % solution of the constitutive problem
    
    % vector of internal forces
    F(:) = B'*reshape(repmat(WEIGHT,6,1).*S, 6*n_int,1) ;
    b_N = f(Q)-F(Q);
    
    % stopping criterion
    criterion = norm(b_N)/norm(f(Q));
    crit_hist(it)=criterion;
    if  criterion < 1e-12
        crit_hist=crit_hist(1:it);
        omega_hist=omega_hist(1:it-1);
        fprintf(' Quasi-Newton method 2 converges ');
        fprintf(' number of iteration=%d  ',it);
        fprintf(' stopping criterion=%e  ',criterion);
        fprintf('\n');
        break
    end
    
    % test on number of iteration
    if  it == it_max
        omega_hist=omega_hist(1:it-1);
        fprintf('     Quasi-Newton method 2 converges slowly: stopping criterion=%e  ',criterion)
        fprintf('\n');
        break
    end
    
    % deflated pcg solver with incomplete cholesky preconditioner
    tol = 1e-1;
    max_iter = 100;   
    [dU(Q), iter, resvec, flag,krylov_recycle] = ...
        DCG_adaptive( A_N, b_N, b_N*0,@(x)L'\(L\x),tol,max_iter,krylov_recycle);
    if it == 1
       krylov_recycle.max_size = min(4*iter,krylov_recycle.max_size);
    end

    % update of the unknown vector
    omega=2/(M+m);
    omega_hist(it)=omega;
    U= U + omega*dU ;
    
  end % true

end % function