
function [S,DS,m,M]=constitutive_problem_quasi1(E,shear,bulk,model_type)

% =========================================================================
%
% The aim of this function is to construct constitutive operator and
% auxiliary arrays for the quasi-Newton method at integration points 
% 1,2,...,n_int. These auxiliary arrays are related to Preconditioner 1 
% and have similar structure for Models 1-3.
% 
% Input data:
%  E       - current strain tensor, size(E)=(3,n_int)
%  shear   - shear moduli at integration points, size(shear)=(1,n_int)
%  bulk    - bulk moduli at integration points, size(bulk)=(1,n_int)
%  model_type - type of the nonlinear constitutive model
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
  IOTA=[1;1;1;0;0;0];  
  VOL=IOTA*IOTA'; 
  DEV=diag([1,1,1,1/2,1/2,1/2])-VOL/3; 
 
%
% Deviatoric strain and its norm
%
  dev_E=DEV*E;                 % deviatoric part of E
  z=max(0,sum(E.*dev_E));      % scalar product of the deviatoric strain
  
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
% The stress tensor
%
  S=repmat(bulk,6,1).*(VOL*E)+2*repmat(mu,6,1).*dev_E;   

%
% Auxiliary array for the quasi-Newton stiffness matrix
%
  delta=0.5; % a parameter defining the quasi-Newton matrix
  mu_delta=mu+2*delta*Dmu;
  DS=2*mu_delta;

%
% Parameters m and M defining the quasi-Newton step
%
  mu_1=mu+2*Dmu;
  mu_max=max(mu,mu_1);
  mu_min=min(mu,mu_1);
  m=min(mu_min./mu_delta);
  M=max(mu_max./mu_delta);
 
 end
