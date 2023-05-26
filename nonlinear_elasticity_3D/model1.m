function [mu, Dmu]=model1(z,shear)

% =========================================================================
%
% Model 1 on non-linear elasticity - this model reduces the shear modulus
%
% Input data:
%  z       - given values at each integration point; sqrt(z) represents the
%            norm of the deviatoric strain, size(z)=(1,n_int)
%  shear   - shear modulus at integration points, size(shear)=(1,n_int)
%
% Output data:
%  mu     - modified shear modulus at integr. points, size(mu)=(1,n_int);
%           mu=mu(z) - values of a scalar function representing Model 3
%  Dmu    - size(mu)=(1,n_int); Dmu=z*mu'(z)
%
% =========================================================================
%

% other parameters of Model 1
shear0=0.1*shear;   % reduced shear modulus
eps=10000;          % parameter defining slope of the function mu(z)

% the array mu
mu=shear0+(shear-shear0)./(1+eps.*sqrt(z));

% the array Dmu
Dmu=-eps.*sqrt(z).*(shear-shear0)./(2*(1+eps.*sqrt(z)).^2);

end % function
     