%  the gamma functions
% gamma_formula:
% BJ1978  TG1983  BS1985  N1990
% R2003: Ruessink et al. (2003)
% A2008_JB2007  A2008_TG1983  A2008_B1998  A2008_BJ1978
% Sal2015: Salmon et al. (2015); 
% Lin2017: Lin and Sheng (2017);
% ZL2020: Zhang and Li et al. (2020)
% Mad1976: Madsen (1976); 
% RS2000: Rattanapitikon and Shibayama (2000);
% TM2002: Tajima and Madsen (2002); 
% T2001: Ting (2001);
% SH1985: Sallenger and Holman (1985); 
% SH1989: Sallenger and Howd (1989);
% Lip1996: Lippmann (1996); 
% G2010: Goda (2010)
% NE1987

function [gamma] = gamma_functions(gamma_formula,s0,H0,kh,p,h)

if strcmpi(gamma_formula,'ZL2020') == 1   
    kh_range = [0.3 1.2];  
    s0_range = [0.005 0.05];
    gamma_range = [0 1.2]; % 0 1   0.718            
    s0=max(s0,s0_range(1)); s0=min(s0,s0_range(2)); 
    kh=max(kh,kh_range(1)); kh=min(kh,kh_range(2));
    gamma = (237*s0^2-34.81*s0+1.46)*exp((1.96*log(38.64*s0))*kh);


elseif strcmpi(gamma_formula,'R2003') == 1;
    gamma=0.76*(kh)+0.29;
    
elseif strcmpi(gamma_formula,'Sal2015') == 1
    if kh<=0.75
      gamma=0.54+7.59*min(abs(p),0.1);
    elseif kh>0.75
      gamma1=0.54+7.59*min(abs(p),0.1);
      gamma2=-8.06+8.09*1.33*kh;    
      gamma=gamma1./tanh(gamma1./gamma2);
    end
    if p==0
      gamma=0.95;
    end
                  
elseif strcmpi(gamma_formula,'Lin2017') == 1
    if kh<=0.75
      gamma=0.54+0.47*exp(-0.018/min(abs(p),0.1));
    elseif kh>0.75
      gamma1=0.54+0.47*exp(-0.018/min(abs(p),0.1));
      gamma2=-8.06+8.09*1.33*kh;    
      gamma=gamma1./tanh(gamma1./gamma2);
    end
    
elseif strcmpi(gamma_formula,'N1990') == 1 %%%%%%%JB2007%%%%%%%%%
   gamma=0.39+0.56*tanh(33*s0); 
   
elseif strcmpi(gamma_formula,'BS1985') == 1%%%%%%%AB2007%%%%%%
   gamma=0.5+0.4*tanh(33*s0);
   
elseif strcmpi(gamma_formula,'A2008_JB2007') == 1 %%%%%%%JB2007_tuned%%%%%%
   gamma=0.11+0.55*tanh(1*H0);

elseif strcmpi(gamma_formula,'TG1983') == 1 %%%%%%%TG1983%%%%%%
    gamma=0.42;
elseif strcmpi(gamma_formula,'A2008_TG1983') == 1 %%%%%%%TG1983_tuned%%%%%%
    gamma=0.18+0.4*tanh(0.9*H0); 
elseif strcmpi(gamma_formula,'SH1985') == 1
    gamma=0.3+3.2*p;
elseif strcmpi(gamma_formula,'SH1989') == 1 
    gamma=0.24+2.7*p;
elseif strcmpi(gamma_formula,'Lip1996') == 1 
    gamma=0.23+1.42*p;
elseif strcmpi(gamma_formula,'A2008_B1998') == 1 %%%%%%%B1998_tuned%%%%%%
    gamma=0.24+0.45*tanh(H0); 
elseif strcmpi(gamma_formula,'A2008_BJ1978') == 1 %%%%%%%BJ1978_tuned%%%%%%
    gamma=0.3+0.45*tanh(0.9*H0);
elseif strcmpi(gamma_formula,'BJ1978') == 1 %%%%%%%BJ1978_swan%%%%%%
    gamma=0.73;
elseif strcmpi(gamma_formula,'Mad1976') == 1
     gamma = 0.72*(1+6.4*p);
elseif strcmpi(gamma_formula,'RS2000') == 1
     gamma = 0.91+5.01*p-11.21*p^2;
elseif strcmpi(gamma_formula,'TM2002') == 1
     gamma = 0.3+4*p;
elseif strcmpi(gamma_formula,'T2001') == 1
     gamma = 0.17+1.53*kh;
elseif strcmpi(gamma_formula,'G2010') == 1
     A=0.17;  L0=H0/s0; hl0=h/L0;
    gamma = (A/hl0)*(1-exp(-1.5*pi*hl0*(1+11*p^(4/3))));
elseif strcmpi(gamma_formula,'NE1987') == 1
    gamma = 0.55+0.88*exp(-0.012/p);
else 
    gamma=str2double(gamma_formula);
end








