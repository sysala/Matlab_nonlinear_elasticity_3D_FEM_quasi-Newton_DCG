function [U, it, crit_hist, omega_hist]=newton_quasi1(U_ini,WEIGHT,K_V,B,f,Q,shear,bulk,model_type)
                               
%--------------------------------------------------------------------------
% The quasi-Newton method with Preconditioner 1 for solution of the system 
%              find U:   F(U)=f
% Linearized systems are solved by DCG with incomplete Cholesky
% preconditioner.
%
% Input data:
%   U_ini - initial choice of U
%   WEIGHT- weight coefficients of integration points, size(WEIGHT)=(1,n_int)
%   K_V   - volumetric part of the elastic stiffness matrix, size(K_V)=(3*n_n,3*n_n)
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
  DS_old=zeros(1,n_int);
  
  % auxiliary arrays for the sparse stiffness matrix assembly
  IOTA=[1;1;1;0;0;0];  
  VOL=IOTA*IOTA'; 
  DEV=diag([1,1,1,1/2,1/2,1/2])-VOL/3; 
  AUX=reshape(1:6*n_int,6,n_int);
  iD=repmat(AUX,6,1); 
  jD=kron(AUX,ones(6,1));
  r=[1 2 3 7 8 9 13 14 15 22 29 36];
  iD=iD(r,:); jD=jD(r,:);
  
  %
  it_max=100;
  crit_hist=zeros(1,it_max);
  omega_hist=zeros(1,it_max);
  deflation_basis = [];
 
%  
% Quasi-Newton's solver (Preconditioner 1)
%
      
  it=0;               % iteration number
   
  while true         
      
     it=it+1;  

     % solution of the constitutive problem
     E(:) = B*U(:) ;      % strain at integration points
     [S,DS_D,m,M]=constitutive_problem_quasi1(E,shear,bulk,model_type);
                          % solution of the constitutive problem
     
     % stiffness matrix related to Preconditioner 1
     % if the matrix change is minimal we let the matrix from previous step
     if max(abs(DS_D-DS_old)./abs(DS_D))>1e-1
        vD=DEV(:)*DS_D;
        vD = repmat(WEIGHT,36,1).*vD ;      
        vD=vD(r,:);
        D_p = sparse( iD(:),jD(:),vD(:), 6*n_int,6*n_int ) ;   
        K_tangent = K_V + B'*D_p*B;  
        A_N = K_tangent(Q,Q);
        A_N = (A_N + A_N')/2;
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
        fprintf(' Quasi-Newton method 1 converges ');
        fprintf(' number of iteration=%d  ',it); 
        fprintf(' stopping criterion=%e  ',criterion); 
        fprintf('\n'); 
        break
     end      

     % test on number of iteration
     if  it == it_max
        omega_hist=omega_hist(1:it-1); 
        fprintf('     Quasi-Newton method 1 converges slowly: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        break
     end     
        
     % deflated pcg solver with incomplete cholesky preconditioner   
     tol = 1e-1;
     max_iter = 1000;
     reg_parameter = 1e-4;
     ichol_params = struct('type', 'ict', 'droptol', 4e-3);
     L = ichol(A_N + reg_parameter*diag(diag(A_N)), ichol_params);
     [dU(Q), iter, resvec, flag] = DCG( A_N, b_N, b_N*0,deflation_basis,@(x)L'\(L\x),tol,max_iter);
     [ W_orth] = my_orth_simple( deflation_basis, dU(Q));
     deflation_basis = [deflation_basis W_orth];

     % update of the unknown vector
     omega=2/(M+m);
     omega_hist(it)=omega;
     U= U + omega*dU ;
     
  end % true     
  
end % function