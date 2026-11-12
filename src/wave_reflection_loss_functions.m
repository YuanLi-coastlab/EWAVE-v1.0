%  (1) calculate the local energy loss of incoming wave due to partially/fully   
%       reflection from the underwater sloping bottom.
%  (2) calculate the reflection coefficient at the shoreline boundary using 
%       the formulas of coastal structure reflection.
%  Kr_formula = 'LL2005' 'M1983' 'L1932' for bottom slope-induced refelction
%               'ZV2008' 'M2017' 'V2005' 'B1974' 'YT2011' for swash-related reflection

%% 
function     [D_reflected,Kr] = wave_reflection_loss_functions (direction,reflection_source,setup_Kr_formula,Fi,H,Tt,c,depth,mwl, x,i0,ncal,dx,p,c1)

    H  = H * sqrt(2) +1e-19;   %
    i1 = i0+direction*1;      i1 = max (i1,1);      i1 = min (i1,ncal);
    Fi0 = Fi (i0);                 H0 = H(i0);            T = Tt (i0);
    h = depth+mwl;          h0 = h(i0);             h1 = h(i1);
    c0 = c(i0);                    L = c.*Tt;                L00 = 9.81*T^2/2/pi;
                          
    correction = setup_Kr_formula(3);  
    if strcmpi(correction,'original') == 1
        v0 = sqrt(9.81*h0);        v1 = sqrt(9.81*h1);        L0 = sqrt(9.81*h0)*T;     L1 = sqrt(9.81*h1)*T;    v   = power(9.81*h,0.5);
    elseif strcmpi(correction,'corrected') == 1
        v0 = c0;                 v1 = c1;                 L0 = c0*T;                      L1 = c1*T;                     v   = c;
    end
   

    if strcmpi(reflection_source,'bottom') == 1            
 
        Kr_formula = setup_Kr_formula(1);   L_fac = str2double (setup_Kr_formula(4));   L_fac = max ( L_fac, dx/L0+1e-19);

        ivec = find(direction*x>=direction*x(i0)-L_fac*L0  &  direction*x<direction*x(i0)); 
        ivec(ivec>ncal)=[];  if isempty(ivec)==1; ivec=1; end
        Kr_0 = zeros(length(x),1); Kr_1 = zeros(length(x),1); D_r = zeros(length(x),1); 
        for i = ivec(1):ivec(end)
            distance0 = abs (x(i) - x(i0)) + 1e-19;   distance1 = abs (x(i) - x(i1)) + 1e-19;
            dhdx0      = (h(i) - h0) / (x(i) - x(i0));    tanb0    = abs (-dhdx0);       tanb0 = max(-dhdx0*direction,0);
            dhdx1      = (h(i) - h1) / (x(i) - x(i1));    tanb1    = abs (-dhdx1);       tanb1 = max(-dhdx1*direction,0);
            
% theoretical formulas   
            if strcmpi(Kr_formula,'L1932') == 1
                Kr_0(i) = Kr_L1932 (v(i),v0);     
                Kr_1(i) = Kr_L1932 (v(i),v1);     
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2) ) / dx * Fi0/(1-Kr_0(i)^2) ); 
            elseif strcmpi(Kr_formula,'LL2005') == 1
                Kr_0(i) = Kr_LL2005 (v(i)/v0,L0,distance0);  
                Kr_1(i) = Kr_LL2005 (v(i)/v1,L1,distance1);
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi0/(1-Kr_0(i)^2) ); 
            elseif strcmpi(Kr_formula,'M1983') == 1
                Kr_0(i) = Kr_M1983 (v(i)^2,v0^2,L(i),L0,distance0);
                Kr_1(i) = Kr_M1983 (v(i)^2,v1^2,L(i),L1,distance1);
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi0/(1-Kr_0(i)^2) ); 
% structure formulas            
            elseif strcmpi(Kr_formula,'YT2011') == 1
                Kr_0(i) = Kr_YT2011 (-h0/H(i));
                Kr_1(i) = Kr_YT2011 (-h1/H(i));
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi(i) );  
            elseif strcmpi(Kr_formula,'ZV2008') == 1
                Kr_0(i) = Kr_ZV2008 (tanb0/sqrt(H(i)/L00),-h0/H(i));   
                Kr_1(i) = Kr_ZV2008 (tanb1/sqrt(H(i)/L00),-h1/H(i));  
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi(i) ); 
            elseif strcmpi(Kr_formula,'M2017') == 1
                Kr_0(i) = Kr_M2017 (tanb0/sqrt(H(i)/L00),-h0/H(i));
                Kr_1(i) = Kr_M2017 (tanb1/sqrt(H(i)/L00),-h1/H(i));
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi(i) );  
            elseif strcmpi(Kr_formula,'V2005') == 1
                Kr_0(i) = Kr_V2005 (tanb0/sqrt(H(i)/L00),-h0/H(i));
                Kr_1(i) = Kr_V2005 (tanb1/sqrt(H(i)/L00),-h1/H(i));
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi(i) ); 
            elseif strcmpi(Kr_formula,'B1974') == 1
                Kr_0(i) = Kr_B1974 (tanb0/sqrt(H(i)/L00),-h0/H(i));
                Kr_1(i) = Kr_B1974 (tanb1/sqrt(H(i)/L00),-h1/H(i));
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi(i) ); 
            else
                Kr_0(i) = str2double(Kr_formula);  
                Kr_1(i) = str2double(Kr_formula);  
                D_r(i)   = max (0 , ( power(Kr_1(i),2)  - power(Kr_0(i),2)  ) / dx * Fi(i) );  
            end  
        end

        imax = find (D_r == max (D_r));
        Kr = Kr_1 (imax(1));      D_reflected = D_r (imax(1));
        
        
        %% test the original formulations %%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        ip=max(i0-1,1);  pi0 = max(p(i0)*direction,0) ;   pip = max(p(ip)*direction,0) ;
        if strcmpi(Kr_formula,'L1932_old') == 1
                Kr = Kr_L1932 (v0,v1);     
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);    
            elseif strcmpi(Kr_formula,'LL2005_old') == 1
                Kr = Kr_LL2005 (v0/v1,L0,dx);  
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);     
            elseif strcmpi(Kr_formula,'M1983_old') == 1
                Kr = Kr_M1983 (v0^2,v1^2,L0,L1,dx);
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);     
% structure formulas            
            elseif strcmpi(Kr_formula,'YT2011_old') == 1
                Kr_0 = Kr_YT2011 (-h0/H(ip));
                Kr_1 = Kr_YT2011 (-h1/H(i0));  Kr = Kr_1;
                D_reflected   = max (0 , ( power(Kr_1,2)  - power(Kr_0,2)  ) / dx * Fi0 );                  
            elseif strcmpi(Kr_formula,'ZV2008_old') == 1
                Kr_0 = Kr_ZV2008 (pip/sqrt(H(ip)/L00),-h0/H(ip));
                Kr_1 = Kr_ZV2008 (pi0/sqrt(H(i0)/L00),-h1/H(i0));   Kr = Kr_1;               
                D_reflected   = max (0 , ( power(Kr_1,2)  - power(Kr_0,2)  ) / dx * Fi0 );  
            elseif strcmpi(Kr_formula,'M2017_old') == 1
                Kr_0 = Kr_M2017 (pip/sqrt(H(ip)/L00),-h0/H(ip));
                Kr_1 = Kr_M2017 (pi0/sqrt(H(i0)/L00),-h1/H(i0));   Kr = Kr_1;
                D_reflected   = max (0 , ( power(Kr_1,2)  - power(Kr_0,2)  ) / dx * Fi0 );  
            elseif strcmpi(Kr_formula,'V2005_old') == 1
                Kr_0 = Kr_V2005 (pip/sqrt(H(ip)/L00),-h0/H(ip));
                Kr_1 = Kr_V2005 (pi0/sqrt(H(i0)/L00),-h1/H(i0));   Kr = Kr_1;
                D_reflected   = max (0 , ( power(Kr_1,2)  - power(Kr_0,2)  ) / dx * Fi0 ); 
            elseif strcmpi(Kr_formula,'B1974_old') == 1
                Kr_0 = Kr_B1974 (pip/sqrt(H(ip)/L00),-h0/H(ip));
                Kr_1 = Kr_B1974 (pi0/sqrt(H(i0)/L00),-h1/H(i0));   Kr = Kr_1;
                D_reflected   = max (0 , ( power(Kr_1,2)  - power(Kr_0,2)  ) / dx * Fi0 );   
         end  
         % test the original formulations %%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
     
    elseif strcmpi(reflection_source,'swash') == 1   
            Kr_formula = setup_Kr_formula(2); 
            h1 = -2*H0;        pi0 = max(p(i0)*direction,0) ;
% structure formulas            
            if strcmpi(Kr_formula,'YT2011') == 1
                Kr = Kr_YT2011 (-h1/H(i0)); 
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0); 
            elseif strcmpi(Kr_formula,'ZV2008') == 1
                Kr = Kr_ZV2008 (pi0/sqrt(H(i0)/L00),-h1/H(i0));
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);  
            elseif strcmpi(Kr_formula,'M2017') == 1
                Kr = Kr_M2017 (pi0/sqrt(H(i0)/L00),-h1/H(i0));
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);  
            elseif strcmpi(Kr_formula,'V2005') == 1
                Kr = Kr_V2005 (pi0/sqrt(H(i0)/L00),-h1/H(i0));
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);  
            elseif strcmpi(Kr_formula,'B1974') == 1
                Kr = Kr_B1974 (pi0/sqrt(H(i0)/L00),-h1/H(i0));
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);  
            else
                Kr = str2double(Kr_formula);  
                D_reflected   = max ( power(Kr,2)/dx * Fi0  , 0);   
            end  
          
           
    end
    
    
end 

    %% bottom slope-induced refelction 
    function [Kr] = Kr_Z2020 (v0,v1,distance,Kr0)
        dcdx = abs (v1-v0) / distance;
        dKrdx = dcdx * (1-Kr0^2) / (2*v0);
        Kr = Kr0 + distance * dKrdx;   
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
    
  
    %% bottom slope-induced refelction 
    function [Kr] = Kr_L1932 (v0,v1)
        Kr = abs (  (v0-v1)/(v0+v1)  );
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
        
        
    %% bottom slope-induced refelction   
    function [Kr] = Kr_M1983 (v0,v1,L0,L1,distance)  
        vm = 0.5*(v0+v1); Lm =0.5*(L0+L1) * sqrt(v0+v1)/sqrt(0.5)/(sqrt(v0)+sqrt(v1)); km = 2*pi/Lm;
        Kr = abs (  (v0-v1)/(4*vm) * sin(km*distance)/(km*distance)  );
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
    
    %% bottom slope-induced refelction   an infinite step behind a linear slope
    function [Kr] = Kr_LL2005 (s01,L0,distance)  
    if s01~=1
        if s01 < 1; s01=power(s01,-1); end
        M1 = distance/L0; 
        alpha0 = 4*pi*M1*(s01^2)/(s01^2-1);     
        alpha1 = 4*pi*M1*(s01^1)/(s01^2-1);
        delta0 = besselj(0,alpha0)^2  + besselj(1,alpha0)^2;
        delta1 = besselj(0,alpha1)^2  + besselj(1,alpha1)^2;
        miu0   = bessely(0,alpha0)^2 + bessely(1,alpha0)^2;
        miu1   = bessely(0,alpha1)^2 + bessely(1,alpha1)^2;
        v0       = besselj(0,alpha0)*bessely(0,alpha0) + besselj(1,alpha0)*bessely(1,alpha0);
        v1       = besselj(0,alpha1)*bessely(0,alpha1) + besselj(1,alpha1)*bessely(1,alpha1);
        % J0  and J1 are the Bessel function of the    first    kind of order 0 and 1;
        % Y0 and Y1 are the Bessel function of the second kind of order 0 and 1;
        Kr =   sqrt( delta0*miu1 + delta1*miu0 - 2*v0*v1 - 8*s01/pi^2/alpha0^2 ) ...
                / sqrt( delta0*miu1 + delta1*miu0 - 2*v0*v1 + 8*s01/pi^2/alpha0^2 );  
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    else
        Kr = 0;
    end
    end
        
    
    %% bottom slope/swash-related reflection 
    function [Kr] = Kr_YT2011 (RH)  
        Kr = 0.53 * exp (+0.85*RH);  
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
 
        
    %%  swash-related reflection    
    function [Kr] = Kr_ZV2008 (ir,RH)   
         a = 0.16; b = 1.43;   % the coefficients for smooth sloping bed
         if RH > 0.5
            Kr = tanh (a*ir^b) ;
         elseif RH <= 0.5
            Kr =tanh (a*ir^b) * (0.67 + 0.37*RH);
         end
         Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
         
         
    %%  swash-related reflection    
    function [Kr] = Kr_M2017 (ir,RH)  
        Kr = 0.3*ir^2/(9.96+ir^2) * (1.271 + 0.177*RH); 
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
        
        
    %% swash-related reflection   
    function [Kr] = Kr_V2005 (ir,RH)  
        if RH >= 0.5
            Kr = 0.14 * ir^0.73;
        elseif RH < 0.5
            Kr = 0.14 * ir^0.73 * (0.9 + 0.2*RH);
        end
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end
        
        
    %%  swash-related reflection   
    function [Kr] = Kr_B1974 (ir,RH)   
        if RH >= 0.5
            Kr = 0.1 * ir^2;
        elseif RH < 0.5 
            Kr = 0;
        end
        Kr = real(Kr);  Kr = min (Kr, 1); Kr = max (Kr, 0); 
    end  