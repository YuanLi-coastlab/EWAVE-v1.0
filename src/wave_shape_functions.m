%  the wave shape factor functions (Zhang et al. 2018)

function [shape] = wave_shape_functions(B_formula,Ursell)
if strcmpi(B_formula,'Z2018') == 1
    shape =min(0.125, 0.15-0.03*log10(Ursell)); 
    shape =max(0.04,shape);
else 
    shape =str2double(B_formula);  
end







