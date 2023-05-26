function [mu, Dmu]=model3(z,shear)

% =========================================================================
%
% Model 3 on non-linear elasticity - this model represent a smooth version
% of the elastoplastic oprerator with isotropic hardening; 
%
% Input data:
%  z       - given values at each integration point; sqrt(z) represents the
%            norm of the deviatoric strain, size(z)=(1,n_int)
%  shear   - shear moduli at integration points, size(shear)=(1,n_int)
%
% Output data:
%  mu     - modified shear moduli at integr. points, size(mu)=(1,n_int);
%           mu=mu(z) - values of a scalar function representing Model 3
%  Dmu    - size(mu)=(1,n_int); Dmu=z*mu'(z)
%
% =========================================================================
%

% other parameters of Model 3
alpha=0.5;            % hardening parameter
gamma=100;            % yield stress
eps=0.01*gamma;       % regularization parameter

% branching within the definition of Model 3
test=2*shear.*sqrt(z);
ind2=(test>gamma-eps)&(test<gamma+eps);
ind3=(test>=gamma+eps);

% the array mu
mu=shear;
mu(ind2)=(1-alpha)*shear(ind2)+...
    alpha*(gamma-(2*shear(ind2).*sqrt(z(ind2))-gamma-eps).^2/(4*eps))./(2*sqrt(z(ind2)));
mu(ind3)=(1-alpha)*shear(ind3)+alpha*gamma./(2*sqrt(z(ind3)));

% the array Dmu
Dmu=zeros(1,length(shear));
Dmu(ind2)=-alpha*(gamma-(2*shear(ind2).*sqrt(z(ind2))-gamma-eps).^2/(4*eps))./(4*sqrt(z(ind2)))+...
    -alpha*shear(ind2).*(2*shear(ind2).*sqrt(z(ind2))-gamma-eps)/(4*eps);
Dmu(ind3)=-alpha*gamma./(4*sqrt(z(ind3)));

end % function
     