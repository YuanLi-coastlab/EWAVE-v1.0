%  calculate the offshore wave parameters 

function [H0,L0,s0] = offshore_wave_functions (Hsw,Tsw,fwx,shape,L,s0_method)
grav=9.81; pi=3.14159; rho=1000;
if strcmpi(s0_method,'deepwater') == 1
    L0=grav*Tsw^2/2/pi; 
    c0=L0/Tsw; 
    cg0=c0/2; 
    H0 = sqrt(fwx/shape/rho/grav/cg0); 
    s0=H0/L0; 
elseif strcmpi(s0_method,'offshore') == 1
    L0=L; 
    c0=L0/Tsw; 
    H0 = Hsw; 
    s0=H0/L0;
end