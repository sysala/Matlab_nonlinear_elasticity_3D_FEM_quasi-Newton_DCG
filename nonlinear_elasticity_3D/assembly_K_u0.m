function [K_u0,Dmu]=assembly_K_u0(U0,B,WEIGHT,shear,model_type)

% =========================================================================
%
% The aim of this function is to assemble an auxiliary stiffness matrix 
% generated using the vector u0.
%
% Input data:
%  U0      - given displacement vector, size(U)=(3,n_n)
%   B      - the strain-displacement matrix, size(B)=(6*n_int,3*n_n)
%  shear   - shear moduli at integration points, size(shear)=(1,n_int)
%  model_type - type of the nonlinear constitutive model
%
% Output data:
%  K_u0 - auxiliary stiffness matrix, size(K_u0)=(3*n_n,3*n_n)
%  Dmu  - derivative of the function mu
%
% =========================================================================
%

% Strain tensor related to U0
  n_int=length(shear);
  E = zeros(6,n_int);
  E(:) = B*U0(:) ;   % strain at integration points

%
% Deviatoric and volumetric 6x6 matrices
% 
  IOTA=[1;1;1;0;0;0];  
  VOL=IOTA*IOTA'; 
  DEV=diag([1,1,1,1/2,1/2,1/2])-VOL/3; 
 
%
% Deviatoric strain and its norm
%
  dev_E=DEV*E;                 % deviatoric part of E
  z=max(0,sum(E.*dev_E));      % scalar product of the deviatoric strain
  s=sqrt(z);                   % norm of the deviatoric strain
  
%
% Choice of the nonlinear constitutive model 
%

 switch(model_type)
    case 'M1'
        [mu, Dmu]=model1(z,shear);
    case 'M2'
        [mu, Dmu]=model2(z,shear);
    case 'M3'
        [mu, Dmu]=model3(z,shear);    
    otherwise
        disp('bad choice of model type');
  end           

%
% Constitutive matrix representing U0, size(DS)=(36,n_int)
%

  DS=zeros(36,n_int);
  IND=(s>0);
  N_hat=dev_E(:,IND)./repmat(s(IND),6,1);
  NN_hat=repmat(N_hat,6,1).*kron(N_hat,ones(6,1));
  DS(:,IND)=DS(:,IND)+4*repmat(Dmu(IND),36,1).*NN_hat; 

%
% Assembling of the matrix K_u0
%
  AUX=reshape(1:6*n_int,6,n_int);
  iD=repmat(AUX,6,1); 
  jD=kron(AUX,ones(6,1));
  vD = repmat(WEIGHT,36,1).*DS ;      
  D_p = sparse( iD(:),jD(:),vD(:), 6*n_int,6*n_int ) ;   
  K_u0 = B'*D_p*B;   

 end
