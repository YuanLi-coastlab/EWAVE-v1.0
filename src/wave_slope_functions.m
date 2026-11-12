%  the wave front slope functions (Zhang et al. 2017)

function [beta] = wave_slope_functions(beta_formula,Ursell,HL)
if strcmpi(beta_formula,'Z2017') == 1
    beta=(Ursell<=1)*atan(2*HL)+(Ursell>1)*atan((2+0.6*(log10(Ursell))^3)*HL); 
else 
    beta=str2double(beta_formula);  
end







