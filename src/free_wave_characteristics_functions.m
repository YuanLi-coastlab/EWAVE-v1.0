%  calculate the local wave parameters 

function [k,c,L,n,kh,HL,Hh,Ursell,shape,beta,ew,cg,fwx,sxx,qwx] ...
= free_wave_characteristics_functions (Hsw,Tsw,angle,depth,mwl,Uc,...
      c_formula,B_formula,beta_formula)
    grav = 9.81; rho = 1000;
    hh = depth+mwl+strcmpi(c_formula,'nonlinear')*(0.5*Hsw); 

    UU = 0; omega=2*pi/Tsw; 
    [k,~] = wave_dispersion (omega,grav,hh,UU,1e-3,0);
    c=2*pi/Tsw/k;
    L=c*Tsw;
    n=0.5*( 1+ (2*omega*(depth+mwl)/c) / sinh(2*omega*(depth+mwl)/c)  );
    kh = k*(depth+mwl);
    HL=Hsw/L;
    Hh=Hsw/(depth+mwl);
    Ursell=Hsw*L^2/(depth+mwl)^3;

    [shape] = wave_shape_functions(B_formula,Ursell);
    [beta] = wave_slope_functions(beta_formula,Ursell,HL);

    ew=shape*rho*grav*Hsw^2;
    cg=c*n;
    fwx=ew*(cg*cos(angle)+Uc);
    sxx=ew*((1+cos(angle)^2)*n-0.5);
    qwx=ew/rho/c*cos(angle);   