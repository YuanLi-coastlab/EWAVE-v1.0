%  calculate the bottom friction-induced wave energy dissipation rate
 % Jonsson 1966  Baldock 2012  Renier 2004  TG1983...
function [df,ubx,usx] = wave_friction_dissipation_functions (Hsw,Tsw,dd50,kh,angle)
grav=9.81; pi=3.14159; rho=1000;
ks=2.5*dd50; z0=ks/33+1e-19;
AA=0.5*Hsw/sinh(kh)+1e-19; 
ub=pi*Hsw/Tsw/sinh(kh);
us=pi*Hsw/Tsw/sinh(kh)*cosh(kh);
fricf=1.39*power(AA/z0,-0.52);
if z0==0; fricf=0; end
df = 0.5/sqrt(pi)*rho*fricf*power(ub,3);
ubx=ub*cos(angle);
usx=us*cos(angle);