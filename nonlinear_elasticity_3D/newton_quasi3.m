function [U, it, crit_hist, omega_hist]=newton_quasi3(U_ini,WEIGHT,B,B_V,B_grad,f,Q,shear,bulk,model_type)

%--------------------------------------------------------------------------
% The quasi-Newton method with Preconditioner 4 for solution of the system
%              find U:   F(U)=f
% Linearized systems are solved by DCG with a block-diagonal
% preconditioner.
%
% Input data:
%   U_ini - initial choice of U
%   WEIGHT - weight coefficients of integration points, size(WEIGHT)=(1,n_int)
%   B - the gradient matrix, size(B)=(6*n_int,3*n_n)
%   B_V    - matrix representing volumetric strain, size(B_V)=(n_int,3*n_n)
%   B_grad - matrix representing gradient of U, size(B_grad)=(9*n_int,3*n_n)
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

  n_int=length(WEIGHT); % number of integration points
  n_n=size(U_ini,2);    % number of nodes
  dU = zeros(3,n_n) ;   % Newton's increment (in displacement)
  F = zeros(3,n_n) ;    % vector of internal forces
  E = zeros(6,n_int);   % values of the strain tensor at integration points
  U=U_ini;              % initial approximation of the solution
  DS_old=zeros(1,n_int);% array for the freezing technique

  % arrays with indices for block diagonal preconditioner
  Q1=Q(1,:); Q2=Q(2,:); Q3=Q(3,:);
  idx_all = reshape(1:(3*size(Q,2)),size(Q,2),3)'*0;
  idx_all(Q) = 1:(sum(Q(:)));
  idx1 = idx_all(1,Q1);
  idx2 = idx_all(2,Q2);
  idx3 = idx_all(3,Q3);

  % global variables for inner preconditioners
  global deflation_all P1

  %
  deflation_all=[];
  deflation_basis = [];
  
  %
  it_max=100;
  crit_hist=zeros(1,it_max);
  omega_hist=zeros(1,it_max);

%  
% Quasi-Newton's solver (Preconditioner 1)
%

  it=0;               % iteration number
  while true
    
    it=it+1;    
    
    % solution of the constitutive problem
    E(:) = B*U(:) ;   % strain at integration points
    [S,DS_D,m,M]=constitutive_problem_quasi3(E,shear,bulk,model_type);
                      % solution of the constitutive problem
    
    % stiffness matrix related to Preconditioner 3
    % if the matrix change is minimal we let the matrix from previous step
    if max(abs(DS_D-DS_old)./abs(DS_D))>1e-1
      D_V = sparse( 1:n_int,1:n_int,WEIGHT.*(bulk-DS_D/3));
      D_D = sparse( 1:9*n_int,1:9*n_int,kron(WEIGHT.*DS_D,[1 1 1 1 1 1 1 1 1])) ;
      K_tangent = B_V'*D_V*B_V + B_grad'*D_D*B_grad;   
      A_N = K_tangent(Q,Q);
      A_N = (A_N + A_N')/2;
      % diagonal preconditioner set-up     
      A_N1=A_N(idx1,idx1);
      A_N2=A_N(idx2,idx2);
      A_N3=A_N(idx3,idx3);
      precond = @(x)prec_blkdiag(x,A_N1,A_N2,A_N3,idx1,idx2,idx3);
      P1 =[]; % resseting cholesky preconditioners for each block diag system
      DS_old=DS_D;
    end
       
    % vector of internal forces
    F(:) = B'*reshape(repmat(WEIGHT,6,1).*S, 6*n_int,1) ;
    b_N = f(Q)-F(Q);
    
    % stopping criterion
    criterion = norm(b_N)/norm(f(Q));
    crit_hist(it)=criterion;
    if  criterion < 1e-12
        crit_hist=crit_hist(1:it);
        omega_hist=omega_hist(1:it-1);
        fprintf(' Quasi-Newton method 3 converges ');
        fprintf(' number of iteration=%d  ',it);
        fprintf(' stopping criterion=%e  ',criterion);
        fprintf('\n');
        break
    end
    
    % test on number of iteration
    if  it == it_max
        omega_hist=omega_hist(1:it-1);
        fprintf('\n     Quasi-Newton method 3 converges slowly: stopping criterion=%e  \n',criterion)
        break
    end        
    
    % solver with block-diagonal preconditioner
    tol = 1e-1;
    max_iter = 100; 
    [dU(Q), iter, resvec, flag] = DCG( A_N, b_N, b_N*0,deflation_basis,precond,tol,max_iter);
    if iter>3
        [ W_orth] = my_orth_simple( deflation_basis, dU(Q));
        deflation_basis = [deflation_basis W_orth];
    end
    
    % update of the unknown vector
    omega=2/(M+m);
    omega_hist(it)=omega;
    U= U + omega*dU ;

  end % true

end % function