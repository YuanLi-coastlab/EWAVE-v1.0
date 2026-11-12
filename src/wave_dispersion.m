%  wave theory dispersion relationship with linear shear current
% 
%   NOTE: there is no provision for a stopping current:
%         be carefull with large U 
%
%
% Call:  [k,error] = lin_disp(omega,g,h,U,tol,linear_shear_current)
%
%       k:     Wave number    (1/L)
%       omega: Wave frequency (1/T)
%       g:     Gravitational acceleration (defines units)
%       h:     Water depth
%       U:     Eulerian current (U>0: assisting, U<0: opposing)
%       error: error in disp relation (gk tanh(kh)-(omega-kU)^2) / k
%       tol:   max allowable error
%       linear_shear_current: 0=vertically uniform current; 1=linear shear curent
%       if no U is given, it is taken as zero
%       if no tol is given, default is 10*eps (~2e-15 on a ieee machine)
%                                      1e-14 in mex file version


function [k,err] = wave_dispersion(omega,g,h,U,tol,linear_shear_current)

maxiter = 20;
if nargin < 5;tol = 10*eps;end
if nargin < 4;U = 0; end
if (nargin < 3)
  disp('error in lin_disp: must have 3 or 4 or 5 arguments')
  error(' ')
end 

pi = 3.14159;
if linear_shear_current == 0
    U0 = U;
elseif linear_shear_current == 1
    U0 = 2*U;
end

k = power (tanh(power(h,3/4)),-2/3); c = omega/k;
err = -999; i = 1;
while (abs(err) > tol)&(i<maxiter)
c0 = c;
k = omega/c;
kh = k*h;
c = sqrt(g/k*tanh(kh)-(linear_shear_current>0)*U0/kh*(c-U0)*tanh(kh))+U0;
err = (c-c0)/c0;
i = i+1;
end
k = omega/c;
k = real(k);

if i>=maxiter
  disp('lin_disp has not converged, check tol and stopping current')
%   error(' ')
end
if k<0
  disp('lin_disp has computed negative k, check for stopping current')
  error(' ')
end  












