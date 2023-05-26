function [mu, Dmu]=model2(z,shear)

% =========================================================================
%
% Model 2 on non-linear elasticity - this model enlarges the shear modulus
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

% other parameters of Model 2
shear1=10*shear;   % enlarged shear modulus
eps=10000;         % parameter defining slope of the function mu(z)

% the array mu
mu=shear1-(shear1-shear)./(1+eps.*sqrt(z));

% the array Dmu
Dmu=eps.*sqrt(z).*(shear1-shear)./(2*(1+eps.*sqrt(z)).^2);

end % function
     