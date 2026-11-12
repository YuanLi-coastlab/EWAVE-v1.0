function [output] = EWAVE_functions_2020d (setup, be, bc, index)  

%% constants
pi=3.14159; grav=9.81; rho=1000;                                        %  
nowtime = bc(index.t);                                                            % Current simulation time.
swl=bc(index.swl);                                                                   % Still-water level.
depth=swl-be(:,2);                                                                   % Still-water depth.
nx=length(be);                                                                         % Number of computational grid points.
dxmean=(be(end,1)-be(1,1))/(nx-1);                                       % Mean grid spacing; only uniform cross-shore grids are supported.
p=bed_slope_functions(be,dxmean);                                       % Local bed slope.
dd50=zeros(nx,1);   dd50(:)=str2double(setup.d50);              % Bed roughness / median grain diameter d50; spatially uniform.
rhor=zeros(nx,1);    rhor(:)=str2double(setup.rhor);                % Surface-roller density; spatially uniform.
Cf=zeros(nx,1);       Cf(:)=str2double(setup.Cf);                        % Current friction coefficient; spatially uniform.
mwl=zeros(nx,1);    mwl(1)=bc(index.mwl);                             % Mean water-level displacement (wave setup/setdown).

%% Short-wave variables
% Short-wave period.
Tsw=zeros(nx,1);           Tsw(:)=bc(index.Tsw);

% Short-wave direction.
angle_sw=zeros(nx,1);   angle_sw(1)=bc(index.angle_sw);

% Incoming short waves.
Hsw_in=zeros(nx,1);              Hsw_in(1)=bc(index.Hsw_in);
k_sw_in=zeros(nx,1);              L_sw_in=zeros(nx,1);        kh_sw_in=zeros(nx,1);
c_sw_in=zeros(nx,1);              n_sw_in=zeros(nx,1);        cg_sw_in=zeros(nx,1); 
Hh_sw_in=zeros(nx,1);           HL_sw_in=zeros(nx,1);      Ursell_sw_in=zeros(nx,1);
shape_sw_in=zeros(nx,1);      ew_sw_in=zeros(nx,1);      gamma_sw_in=zeros(nx,1);
ubx_sw_in=zeros(nx,1);         usx_sw_in=zeros(nx,1);      Kr_sw_inout=zeros(nx,1); 
fwx_sw_in=zeros(nx,1);          sxx_sw_in=zeros(nx,1);      qwx_sw_in=zeros(nx,1);

% Energy source/sink terms.
Dsw_in_br=zeros(nx,1);            Dsw_in_fc=zeros(nx,1);       
Dsw_in_swout=zeros(nx,1); 
Dsw_in_LFin=zeros(nx,1);         % WME exchange between incoming short waves and low-frequency motions (IG waves/setup).

% Outgoing short waves.
Hsw_out=zeros(nx,1);              Hsw_out(1)=bc(index.Hsw_out);
k_sw_out=zeros(nx,1);              L_sw_out=zeros(nx,1);        kh_sw_out=zeros(nx,1);
c_sw_out=zeros(nx,1);              n_sw_out=zeros(nx,1);        cg_sw_out=zeros(nx,1); 
Hh_sw_out=zeros(nx,1);           HL_sw_out=zeros(nx,1);      Ursell_sw_out=zeros(nx,1);
shape_sw_out=zeros(nx,1);      ew_sw_out=zeros(nx,1);      gamma_sw_out=zeros(nx,1);
ubx_sw_out=zeros(nx,1);         usx_sw_out=zeros(nx,1);      Kr_sw_outin=zeros(nx,1);
fwx_sw_out=zeros(nx,1);          sxx_sw_out=zeros(nx,1);      qwx_sw_out=zeros(nx,1);
               
% Energy source/sink terms.
Dsw_out_br=zeros(nx,1);      Dsw_out_fc=zeros(nx,1); 
Dsw_out_swin=zeros(nx,1); 

% Total short-wave height.
Hsw_tt=zeros(nx,1);      Hsw_tt(1)=bc(index.Hsw_tt);

%% IG-wave variables
% IG-wave period.
Tig=zeros(nx,1);             Tig(:)=bc(index.Tig);

% IG-wave direction.
angle_ig=zeros(nx,1);    angle_ig(1)=bc(index.angle_ig);

% Incoming IG waves.
Hig_in=zeros(nx,1);               Hig_in(1)=bc(index.Hig_in);
k_ig_in=zeros(nx,1);              L_ig_in=zeros(nx,1);        kh_ig_in=zeros(nx,1);
c_ig_in=zeros(nx,1);              n_ig_in=zeros(nx,1);        cg_ig_in=zeros(nx,1);
Hh_ig_in=zeros(nx,1);           HL_ig_in=zeros(nx,1);      Ursell_ig_in=zeros(nx,1);
shape_ig_in=zeros(nx,1);      ew_ig_in=zeros(nx,1);       gamma_ig_in=zeros(nx,1);
ubx_ig_in=zeros(nx,1);         usx_ig_in=zeros(nx,1);      Kr_ig_inout=zeros(nx,1); 
fwx_ig_in=zeros(nx,1);         sxx_ig_in=zeros(nx,1);      qwx_ig_in=zeros(nx,1);
phi_ig_in=zeros(nx,1);

% Energy source/sink terms.
Dig_in_br=zeros(nx,1);      Dig_in_fc=zeros(nx,1); 
Dig_in_igout=zeros(nx,1); 

% Outgoing IG waves.
Hig_out=zeros(nx,1);               Hig_out(1)=bc(index.Hig_out);
k_ig_out=zeros(nx,1);              L_ig_out=zeros(nx,1);        kh_ig_out=zeros(nx,1);
c_ig_out=zeros(nx,1);              n_ig_out=zeros(nx,1);        cg_ig_out=zeros(nx,1);
Hh_ig_out=zeros(nx,1);           HL_ig_out=zeros(nx,1);      Ursell_ig_out=zeros(nx,1);
shape_ig_out=zeros(nx,1);      ew_ig_out=zeros(nx,1);       gamma_ig_out=zeros(nx,1);
ubx_ig_out=zeros(nx,1);         usx_ig_out=zeros(nx,1);      Kr_ig_outin=zeros(nx,1); 
fwx_ig_out=zeros(nx,1);         sxx_ig_out=zeros(nx,1);      qwx_ig_out=zeros(nx,1);
phi_ig_out=zeros(nx,1);
                
% Energy source/sink terms.
Dig_out_br=zeros(nx,1);      Dig_out_fc=zeros(nx,1); 
Dig_out_igin=zeros(nx,1); 

% Total IG-wave height.
Hig_tt=zeros(nx,1);      Hig_tt(1)=bc(index.Hig_tt);

%% Surface-roller variables
% Incoming short-wave roller.
rl_sw_in=zeros(nx,1);    rl_sw_in(1)=bc(index.Lr_sw);
er_sw_in=zeros(nx,1);    ra_sw_in=zeros(nx,1);   rh_sw_in=zeros(nx,1);    beta_sw_in=zeros(nx,1);
frx_sw_in=zeros(nx,1);    rxx_sw_in=zeros(nx,1);    qrx_sw_in=zeros(nx,1); 
         
% Energy source/sink terms.
Drsw_in_turb=zeros(nx,1);         
Drsw_in_LFin=zeros(nx,1);         % WME exchange between the incoming short-wave roller and low-frequency motions.

% Incoming IG-wave roller.
rl_ig_in=zeros(nx,1);      rl_ig_in(1)=bc(index.Lr_ig);
er_ig_in=zeros(nx,1);     ra_ig_in=zeros(nx,1);   rh_ig_in=zeros(nx,1);    beta_ig_in=zeros(nx,1);
frx_ig_in=zeros(nx,1);    rxx_ig_in=zeros(nx,1);    qrx_ig_in=zeros(nx,1); 
         
% Energy source/sink terms.
Drig_in_turb=zeros(nx,1);  Dmwl = zeros(nx,1);   

%% Wave-induced currents
% Currents during shoreward propagation.
Usw_in = zeros(nx,1);    % mean current  due to the mass flux of incoming short waves  (should be positive)
Ur_in    = zeros(nx,1);    % mean current  due to the mass flux of incoming rollers   (should be positive)
Uig_in    = zeros(nx,1);   % mean current  due to the mass flux of incoming long waves   (should be positive)
Uc_in     = zeros(nx,1);    % undertow current  due to the mass fluxes of incoming wave and rollers (should be negative)

% Currents during seaward propagation.
Usw_out = zeros(nx,1);     % mean current  due to the mass flux of outgoing short waves  (should be negative)
Ur_out    = zeros(nx,1);     % mean current  due to the mass flux of outgoing rollers  (should be negative)
Uig_out    = zeros(nx,1);   % mean current  due to the mass flux of outgoing long waves  (should be negative)
Uc_out    = zeros(nx,1);    % undertow current  due to the mass fluxes of outgoing wave and rollers (should be positive)

% Total currents.
UL = zeros(nx,1);    % total Laguralr current
US = zeros(nx,1);    % total Stokes drift
UE = zeros(nx,1);    % total undertow current (Euler) ( = Uc_in + Uc_out)
Uy = zeros(nx,1);     %  alongshore current velocity

%% Iteration settings (multiple passes only when reflection is enabled).
if strcmpi(setup.reflected_wave,'yes') == 1
     N_iteration = str2double(setup.N_iteration);    Errors = zeros(N_iteration,1);  
     mwl_in = zeros(nx,1);  mwl_out = zeros(nx,1);
     mwl_final = zeros(nx,N_iteration);
else
     N_iteration = 1;      Errors = zeros(N_iteration,1);  
     mwl_in = zeros(nx,1);  mwl_out = zeros(nx,1);
end

%% Main cross-shore calculation: alternating shoreward/seaward passes until setup convergence.
for iter = 1:N_iteration   

%% Incoming waves
%% Shoreward sweep
for ii=1:nx

dx=be(ii,1)-be(max(ii-1,1),1)+1e-19;  % At the first cell, a small positive dx avoids division by zero.
    

% Incoming short waves.
U = 0+  Uig_in(max(ii-1,1))  +Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  Uig_out(max(ii-1,1))  +  Uc_out(max(ii-1,1)); 
U = U *(strcmpi(setup.current_advection,'yes') == 1);
[k_sw_in(ii),c_sw_in(ii),L_sw_in(ii),n_sw_in(ii),kh_sw_in(ii),HL_sw_in(ii),Hh_sw_in(ii),Ursell_sw_in(ii),shape_sw_in(ii),beta_sw_in(ii),...
     ew_sw_in(ii),cg_sw_in(ii),fwx_sw_in(ii),sxx_sw_in(ii),qwx_sw_in(ii)] ...
= free_wave_characteristics_functions (Hsw_in(max(ii-1,1)),Tsw(max(ii-1,1)),angle_sw(max(ii-1,1)),depth(ii),mwl(max(ii-1,1)),U,...
                                                      setup.c_formula,setup.B_formula,setup.beta_formula);

if ii==1; [H0,~,s0] = offshore_wave_functions (Hsw_in(ii),Tsw(ii),fwx_sw_in(ii),0.125,L_sw_in(ii),setup.s0_method); end   % Evaluate offshore-equivalent wave properties.
if strcmpi(setup.igbc_LHS62(1),'yes') == 1 && bc(index.Hig_tt)==0 && bc(index.Hig_in)==0 && ii==1 % LHS62 incoming IG-wave boundary condition.
%     top = sxx_sw_in(ii);
    top = 0.5*ew_sw_in(ii);  % Modified numerator to avoid unrealistically large offshore IG-wave boundary heights.
    bc(index.Hig_in) = top/rho/(grav*(depth(ii)+mwl(ii))-cg_sw_in(ii)^2)  * str2double(setup.igbc_LHS62(2));  
    Hig_in(1)=bc(index.Hig_in);
end

fwx_sw_in(ii) = fwx_sw_in(max(ii-1,1)) ...    % For shoreward propagation, use the previous-cell source terms; seaward-to-shoreward reflection is local.
                    + dx*    Dsw_out_swin(ii) ...
                    -  dx*(   Dsw_in_br(max(ii-1,1))  +  Dsw_in_fc(max(ii-1,1))  +   Dsw_in_swout(max(ii-1,1))   ...
  + 1*Dsw_in_LFin(max(ii-1,1)) );   % The incoming short-wave WME energy exchange is subtracted here.

if fwx_sw_in(ii)<0   % If energy losses exceed available flux, limit the flux to zero and adjust breaking dissipation.
    Dsw_in_br(max(ii-1,1)) = Dsw_in_br(max(ii-1,1)) + fwx_sw_in(ii)/dx;   
    fwx_sw_in(ii)=0;    
end

angle_sw(ii)=asin(sin(angle_sw(1))/c_sw_in(1)*c_sw_in(ii));
ew_sw_in(ii)=fwx_sw_in(ii)/(cg_sw_in(ii)*cos(angle_sw(ii)) + U);
Hsw_in(ii)=sqrt(ew_sw_in(ii)/shape_sw_in(ii)/rho/grav);
    
U =  0+  Uig_in(max(ii-1,1))  +Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  Uig_out(max(ii-1,1))  +  Uc_out(max(ii-1,1)); 
U = U *(strcmpi(setup.current_advection,'yes') == 1);
[k_sw_in(ii),c_sw_in(ii),L_sw_in(ii),n_sw_in(ii),kh_sw_in(ii),HL_sw_in(ii),Hh_sw_in(ii),Ursell_sw_in(ii),shape_sw_in(ii),beta_sw_in(ii),...
              ew_sw_in(ii),cg_sw_in(ii),~,sxx_sw_in(ii),qwx_sw_in(ii)] ...
= free_wave_characteristics_functions  (Hsw_in(ii),Tsw(ii),angle_sw(ii),depth(ii),mwl(max(ii-1,1)),U,...
                                                      setup.c_formula,setup.B_formula,setup.beta_formula);

                                                  
% Total short-wave height (variance-based combination).
% NOTE: this simple expression neglects the phase-correlation between incident and reflected waves !  
% This error may be large under quasi-standing (resonant) conditions (e.g., in front of a strongly reflective structure)
Hsw_tt(ii)=sqrt(Hsw_in(ii)^2+Hsw_out(ii)^2);   
                                                 
                                                  
% Incoming short-wave roller.
U =  0+  Uig_in(max(ii-1,1))  +Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  Uig_out(max(ii-1,1))  +  Uc_out(max(ii-1,1)); 
U = U *(strcmpi(setup.current_advection,'yes') == 1);
frx_sw_in(ii)=max((frx_sw_in(max(ii-1,1))+Dsw_in_br(max(ii-1,1))*dx-Drsw_in_turb(max(ii-1,1))*dx)...
                    - Drsw_in_LFin(max(ii-1,1))*dx,   0 ) ...
                    *(strcmpi(setup.surface_roller,'yes') == 1);  % Surface-roller switch.
er_sw_in(ii)= frx_sw_in(ii)/(c_sw_in(ii)*cos(angle_sw(ii)) + U);  
ra_sw_in(ii)= 2*er_sw_in(ii)*Tsw(ii)/rhor(ii)/c_sw_in(ii);% Roller cross-sectional area.
rl_sw_in(ii)= sqrt(ra_sw_in(ii)/0.11)*cos(beta_sw_in(ii));% Roller length.
rh_sw_in(ii)= sqrt(ra_sw_in(ii)/0.11)*sin(beta_sw_in(ii));% Roller height.
rxx_sw_in(ii)=2*er_sw_in(ii)*cos(angle_sw(ii))^2;% Roller radiation stress.
qrx_sw_in(ii)=2*er_sw_in(ii)/rhor(ii)/c_sw_in(ii)*cos(angle_sw(ii));


% Incoming IG waves.
if strcmpi(setup.ig_wave,'yes') == 1
    
    U = Usw_in(max(ii-1,1))  +0+  Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  Uig_out(max(ii-1,1))  +  Uc_out(max(ii-1,1)); 
    U = U *(strcmpi(setup.current_advection,'yes') == 1);
    [k_ig_in(ii),c_ig_in(ii),L_ig_in(ii),n_ig_in(ii),kh_ig_in(ii),HL_ig_in(ii),Hh_ig_in(ii),Ursell_ig_in(ii),shape_ig_in(ii),beta_ig_in(ii),...
        ew_ig_in(ii),cg_ig_in(ii),fwx_ig_in(ii),sxx_ig_in(ii),qwx_ig_in(ii)] ...
    = bound_wave_characteristics_functions (Hig_in(max(ii-1,1)),Tig(max(ii-1,1)),Hsw_in(ii),Tsw(ii),angle_ig(max(ii-1,1)),depth(ii),mwl(max(ii-1,1)),U,...
          setup.c_formula,setup.B_formula,setup.beta_formula);
    
    if ii==1; [H0_ig,~,s0_ig] = offshore_wave_functions (Hig_in(ii),Tig(ii),fwx_ig_in(ii),0.125,L_ig_in(ii),setup.s0_method); end  % Evaluate offshore-equivalent IG-wave properties.
      
    fwx_ig_in(ii) = fwx_ig_in(max(ii-1,1)) ... 
                    + dx*(   Dig_out_igin(ii)    )...
                    + dx*( Dsw_in_LFin(max(ii-1,1)) * str2double(setup.SW_IG_transfer) ...
     + Drsw_in_LFin(max(ii-1,1)) * str2double(setup.roller_IG_transfer) ) ... % WME energy contributions from incoming short waves and their rollers.
                    -  dx*(   Dig_in_br(max(ii-1,1))  +  Dig_in_fc(max(ii-1,1))  +   Dig_in_igout(max(ii-1,1))    );                   
    if fwx_ig_in(ii)<0     % If energy losses exceed available IG flux, limit incoming IG flux to zero.
        fwx_ig_in(ii)=0;   
    end
   
    angle_ig(ii)=asin(sin(angle_ig(1))/c_ig_in(1)*c_ig_in(ii));
    ew_ig_in(ii)  = fwx_ig_in(ii)/(cg_ig_in(ii)*cos(angle_ig(ii)) + U); 
    Hig_in(ii) = real(sqrt(ew_ig_in(ii)/shape_ig_in(ii)/rho/grav));
    
    U = Usw_in(max(ii-1,1))  +0+  Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  Uig_out(max(ii-1,1))  +  Uc_out(max(ii-1,1)); 
    U = U *(strcmpi(setup.current_advection,'yes') == 1);
    [k_ig_in(ii),c_ig_in(ii),L_ig_in(ii),n_ig_in(ii),kh_ig_in(ii),HL_ig_in(ii),Hh_ig_in(ii),Ursell_ig_in(ii),shape_ig_in(ii),beta_ig_in(ii),...
        ew_ig_in(ii),cg_ig_in(ii),~,sxx_ig_in(ii),qwx_ig_in(ii)] ...
    = bound_wave_characteristics_functions (Hig_in(ii),Tig(ii),Hsw_in(ii),Tsw(ii),angle_ig(ii),depth(ii),mwl(max(ii-1,1)),U,...
          setup.c_formula,setup.B_formula,setup.beta_formula);
      

% Total IG-wave height (variance-based combination).
% NOTE: this simple expression neglects the phase-correlation between incident and reflected waves !  
% This error may be large under quasi-standing (resonant) conditions (e.g., in front of a strongly reflective structure)
Hig_tt(ii)=sqrt(Hig_in(ii)^2+Hig_out(ii)^2);

      
% Incoming IG-wave roller.
U = Usw_in(max(ii-1,1))  +0+  Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  Uig_out(max(ii-1,1))  +  Uc_out(max(ii-1,1)); 
U = U *(strcmpi(setup.current_advection,'yes') == 1);
frx_ig_in(ii)=max((frx_ig_in(max(ii-1,1))+Dig_in_br(max(ii-1,1))*dx-Drig_in_turb(max(ii-1,1))*dx),0)...
                         *(strcmpi(setup.surface_roller,'yes') == 1);  % Surface-roller switch.
er_ig_in(ii)= frx_ig_in(ii)/(cg_ig_in(ii)*cos(angle_ig(ii)) + U);
ra_ig_in(ii)= 2*er_ig_in(ii)*Tig(ii)/rhor(ii)/c_ig_in(ii);% Roller cross-sectional area.
rl_ig_in(ii)= sqrt(ra_ig_in(ii)/0.11)*cos(beta_ig_in(ii));% Roller length.
rh_ig_in(ii)= sqrt(ra_ig_in(ii)/0.11)*sin(beta_ig_in(ii));% Roller height.
rxx_ig_in(ii)=2*er_ig_in(ii)*cos(angle_ig(ii))^2;% Roller radiation stress.
qrx_ig_in(ii)=2*er_ig_in(ii)/rhor(ii)/c_ig_in(ii)*cos(angle_ig(ii)); 

end


% Currents (cross-shore, depth-integrated transport representation).
Usw_in(ii) = qwx_sw_in(ii)  /  (depth(ii)+mwl(max(ii-1,1)));
Ur_in(ii) = (qrx_sw_in(ii) + qrx_ig_in(ii))  /  (depth(ii)+mwl(max(ii-1,1)));
Uig_in(ii) = qwx_ig_in(ii)  /  (depth(ii)+mwl(max(ii-1,1)));
Uc_in(ii) = - (Usw_in(ii) + Ur_in(ii) + Uig_in(ii));
US(ii) = -Uc_in(ii) - Uc_out(ii);
UE(ii) = Uc_in(ii) + Uc_out(ii);
UL(ii) = US(ii) + UE(ii);


% Mean water-level variation.
term1=(sxx_sw_in(ii)-sxx_sw_in(max(ii-1,1)))/dx;
term2=(rxx_sw_in(ii)-rxx_sw_in(max(ii-1,1)))/dx;
term3=(sxx_ig_in(ii)-sxx_ig_in(max(ii-1,1)))/dx;
term4=(rxx_ig_in(ii)-rxx_ig_in(max(ii-1,1)))/dx;
term5=(sxx_sw_out(ii)-sxx_sw_out(max(ii-1,1)))/dx;
term6=(sxx_ig_out(ii)-sxx_ig_out(max(ii-1,1)))/dx;
term7=rho*Cf(ii)*UE(ii)*sqrt((1.16*ubx_sw_in(ii)/sqrt(2))^2+(UE(ii))^2);   % bottom shear stress due to current Ruessink 2001
term8=rho*(depth(ii)+mwl(max(ii-1,1)))*UL(ii)*(UL(ii)-UL(max(ii-1,1)))/dx;  % horizontal flow advection term
mwl(ii)=mwl(max(ii-1,1)) - (term1+term2+term3+term4+term5+term6+term7+term8)*dx/rho/grav/(depth(ii)+mwl(max(ii-1,1)));
mwl(ii)=mwl(ii)*(strcmpi(setup.mwl_calculation,'yes') == 1);
mwl_in (ii) = mwl(ii);


% Currents updated with the local mean water-level variation.
Usw_in(ii) = qwx_sw_in(ii)  /  (depth(ii)+mwl(ii));
Ur_in(ii) = (qrx_sw_in(ii) + qrx_ig_in(ii))  /  (depth(ii)+mwl(ii));
Uig_in(ii) = qwx_ig_in(ii)  /  (depth(ii)+mwl(ii));
Uc_in(ii) = - (Usw_in(ii) + Ur_in(ii) + Uig_in(ii));     
US(ii) = -Uc_in(ii) - Uc_out(ii);
UE(ii) = Uc_in(ii) + Uc_out(ii);
UL(ii) = US(ii) + UE(ii);


% Incoming short-wave energy dissipation (used in the next cell).
gamma_sw_in(ii) = real(gamma_functions (setup.gamma_formula,s0,H0,kh_sw_in(ii),p(ii),depth(ii)+mwl(ii)));  
Dsw_in_br(ii) = wave_breaking_dissipation_functions (Hsw_tt(ii),Tsw(ii),depth(ii),mwl(ii),...    % Total short-wave height includes the reflected-wave contribution to breaking.
                gamma_sw_in(ii),kh_sw_in(ii),k_sw_in(ii),setup.wave_type,setup.Db_formula);
            
[Dsw_in_fc(ii),ubx_sw_in(ii),usx_sw_in(ii)] = wave_friction_dissipation_functions (Hsw_in(ii),Tsw(ii),dd50(ii),kh_sw_in(ii),angle_sw(ii));

 if strcmpi(setup.reflected_wave,'yes') == 1
     [~,c1,~,~,~,~,~,~,~,~,~,~,~,~,~] ...
= free_wave_characteristics_functions  (Hsw_in(min(ii+1,nx)),Tsw(min(ii+1,nx)),angle_sw(min(ii+1,nx)),depth(min(ii+1,nx)),mwl(min(ii+1,nx)),U,...
                                                      setup.c_formula,setup.B_formula,setup.beta_formula); 
     [Dsw_in_swout(ii),Kr_sw_inout(ii)] =  wave_reflection_loss_functions ...
          (1,'bottom',setup.Kr_formula,fwx_sw_in,Hsw_in,Tsw,c_sw_in,depth,mwl, be(:,1),ii,nx,dxmean,p,c1);  % Shoreward primary wave: +1; seaward primary wave: -1.

      Dsw_in_swout(ii) = min (Dsw_in_swout(ii), fwx_sw_in(ii)/dx);
 end



% Incoming IG-wave energy losses: breaking, friction, reflection and other retained terms.
if strcmpi(setup.ig_wave,'yes') == 1
    gamma_ig_in(ii) = gamma_functions (setup.gamma_formula,s0_ig,H0_ig,kh_ig_in(ii),p(ii),depth(ii)+mwl(ii));  
    Dig_in_br(ii) = wave_breaking_dissipation_functions (Hig_tt(ii),Tig(ii),depth(ii),mwl(ii),...   
                    gamma_ig_in(ii),kh_ig_in(ii),k_ig_in(ii),setup.wave_type,setup.Db_formula);
             
    [Dig_in_fc(ii),ubx_ig_in(ii),usx_ig_in(ii)] = wave_friction_dissipation_functions (Hig_in(ii),Tig(ii),dd50(ii),kh_ig_in(ii),angle_ig(ii));  

    
    if strcmpi(setup.reflected_wave,'yes') == 1   
        [~,c1,~,~,~,~,~,~,~,~,~,~,~,~,~] ...
    = bound_wave_characteristics_functions (Hig_in(min(ii+1,nx)),Tig(min(ii+1,nx)),Hsw_in(min(ii+1,nx)),Tsw(min(ii+1,nx)),angle_ig(min(ii+1,nx)),depth(min(ii+1,nx)),mwl(min(ii+1,nx)),U,...
          setup.c_formula,setup.B_formula,setup.beta_formula);
       [Dig_in_igout(ii),Kr_ig_inout(ii)] = wave_reflection_loss_functions ...
           (1,'bottom',setup.Kr_formula,fwx_ig_in,Hig_in,Tig,c_ig_in,depth,mwl, be(:,1),ii,nx,dxmean,p,c1);  % Shoreward primary wave: +1; seaward primary wave: -1.
       Dig_in_igout(ii) = min (Dig_in_igout(ii), fwx_ig_in(ii)/dx);
    end                               
end               

                
% Short-wave roller dissipation.
Drsw_in_turb(ii)= 2*grav*er_sw_in(ii)*sin(beta_sw_in(ii))/c_sw_in(ii);     


% IG-wave roller dissipation.
if strcmpi(setup.ig_wave,'yes') == 1
    Drig_in_turb(ii)= 2*grav*er_ig_in(ii)*sin(beta_ig_in(ii))/c_ig_in(ii);          
end


% WME transfer from incoming short waves and rollers to IG waves and wave setup.
Drsw_in_LFin (ii) = c_sw_in(ii)^2/2/Tsw(ii)*(rhor(ii)*ra_sw_in(ii)-rhor(max(ii-1,1))*ra_sw_in(max(ii-1,1)))/dx;

% Setup-related potential-energy gradient; the full flux divergence is not evaluated.
i0 = max(ii-1,1);

if ii == 1
    Dmwl(ii) = 0;
else
    eta_i = mwl(ii);
    eta_0 = mwl(i0);

    % Total depth is used only for roller transport velocity
    H_i = depth(ii) + eta_i;

    % Short-wave and short-wave roller transport velocity
    U_mwl_i = Usw_in(ii) + Usw_out(ii) ...
            + qrx_sw_in(ii)/H_i;

    % Spatial gradient of setup-related potential energy
    dE_mwl_dx = 0.5*rho*grav ...
              *(eta_i^2 - eta_0^2)/dx;

    % Revised Dmwl
    Dmwl(ii) = U_mwl_i*dE_mwl_dx;
end


Dsw_in_LFin(ii) = ...
    (Dmwl(ii) ...
    - (1-str2double(setup.roller_IG_transfer)) ...
    *Drsw_in_LFin(ii)) ...
    /(1-str2double(setup.SW_IG_transfer));

Dsw_in_LFin(ii)  =  0.7 * Dsw_in_LFin(ii)  + 0.3 * Dsw_in_LFin(max(ii-1,1));
Drsw_in_LFin(ii) =  0.7 * Drsw_in_LFin(ii) + 0.3 * Drsw_in_LFin(max(ii-1,1));

% Detect the last active shoreward grid cell.
if (k_sw_in(ii)*(depth(min(ii+1,nx))+mwl(ii))<=str2double(setup.kh_cut) || depth(min(ii+1,nx))+mwl(ii)<=0.001)
    ncal = ii;     break;
else
    ncal = nx;  
end

end

%% Wave reflection
if strcmpi(setup.reflected_wave,'yes') == 1
%% Shoreline boundary: plane beach face; full runup or a non-overtopped structure.
[Dsw_in_swout(ncal),Kr_sw_inout(ncal)]   =  wave_reflection_loss_functions ... 
    (1,'swash',setup.Kr_formula,fwx_sw_in,Hsw_in,Tsw,c_sw_in,depth,mwl, be(:,1),ncal,ncal,dxmean,p,0);  % Shoreward primary wave: +1; seaward primary wave: -1.
Hsw_out(ncal) = Kr_sw_inout(ncal) * Hsw_in(ncal);    % Prescribe the reflected short-wave height at the shoreline boundary.
Dsw_in_br (ncal)   = (fwx_sw_in(ncal)*(1-Kr_sw_inout(ncal)^2))/dx;    % Residual incoming short-wave energy loss at the swash boundary.
Drsw_in_turb (ncal) = (frx_sw_in(ncal) + Dsw_in_br (ncal)*dx)/dx;  % Dissipate the remaining surface-roller energy at the swash boundary.
Dsw_in_fc (ncal)   = 0;             % Neglect bottom-friction loss in the swash zone.
Dsw_in_LFin (ncal) = 0;        % Disable short-wave WME exchange in the swash zone.
Drsw_in_LFin (ncal) = 0;        % Disable roller WME exchange in the swash zone.

[Dig_in_igout(ncal),Kr_ig_inout(ncal)]   =  wave_reflection_loss_functions ... 
    (1,'swash',setup.Kr_formula,fwx_ig_in,Hig_in,Tig,c_ig_in,depth,mwl, be(:,1),ncal,ncal,dxmean,p,0);  % Shoreward primary wave: +1; seaward primary wave: -1.
Hig_out(ncal) = Kr_ig_inout(ncal) * Hig_in(ncal);
Dig_in_br (ncal)   = (fwx_ig_in(ncal)*(1-Kr_ig_inout(ncal)^2))/dx;   
Drig_in_turb (ncal) = (frx_ig_in(ncal) + Dig_in_br (ncal)*dx)/dx;
Dig_in_fc (ncal)   = 0;             
    
% Update shoreline setup, including the contribution from reflected waves (subject to relaxation).
dx = (depth(ncal)+mwl(ncal))/(p(ncal)+1e-19)+1e-19;
sxx_sw_out(ncal)=(shape_sw_in(ncal)*rho*grav*Hsw_out(ncal)^2)*((1+cos(angle_sw(ncal))^2)*n_sw_in(ncal)-0.5);
sxx_ig_out(ncal)=(shape_ig_in(ncal)*rho*grav*Hig_out(ncal)^2)*((1+cos(angle_ig(ncal))^2)*n_ig_in(ncal)-0.5);
term1=(0-sxx_sw_in(ncal))/dx;
term2=(0-rxx_sw_in(ncal))/dx;
term3=(0-sxx_ig_in(ncal))/dx;
term4=(0-rxx_ig_in(ncal))/dx;
term5=(0-sxx_sw_out(ncal))/dx;
term6=(0-sxx_ig_out(ncal))/dx;
term7=rho*Cf(ncal)*UE(ncal)*sqrt((1.16*ubx_sw_in(ncal)/sqrt(2))^2+(UE(ncal))^2);   
term8=rho*(depth(ncal)+mwl(ncal))*UL(ncal)*(0-UL(ncal))/dx;  
mwl(ncal) = mwl(ncal)  -  (term1+term2+term3+term4+term5+term6+term7+term8)*dx/rho/grav/(depth(ncal)+mwl(ncal));
mwl(ncal)=mwl(ncal)*(strcmpi(setup.mwl_calculation,'yes') == 1);
mwl_in (ncal) = mwl(ncal);

%% Seaward sweep
for ii=ncal:-1:1
    
dx=be(min(ii+1,ncal),1)-be(ii,1)+1e-19;


% Outgoing short waves.
U = Usw_in(min(ii+1,ncal))+  Uig_in(min(ii+1,ncal))  +Uc_in(min(ii+1,ncal))  +  0  +  Uig_out(min(ii+1,ncal))  +  Uc_out(min(ii+1,ncal)); 
U = -U *(strcmpi(setup.current_advection,'yes') == 1);  % Reverse the current-advection sign for seaward propagation.
[k_sw_out(ii),c_sw_out(ii),L_sw_out(ii),n_sw_out(ii),kh_sw_out(ii),HL_sw_out(ii),Hh_sw_out(ii),Ursell_sw_out(ii),shape_sw_out(ii),~,...
     ew_sw_out(ii),cg_sw_out(ii),fwx_sw_out(ii),sxx_sw_out(ii),qwx_sw_out(ii)] ...
= free_wave_characteristics_functions (Hsw_out(min(ii+1,ncal)),Tsw(min(ii+1,ncal)),angle_sw(ii),depth(ii),mwl(ii),U,...
                                                      setup.c_formula,setup.B_formula,setup.beta_formula);

fwx_sw_out(ii) = fwx_sw_out(min(ii+1,ncal)) ...   % For seaward propagation, use the next-cell source terms; shoreward-to-seaward reflection is local.
                    + dx*    Dsw_in_swout(ii) ...
                    -  dx*(   Dsw_out_fc(min(ii+1,ncal))  +   Dsw_out_swin(min(ii+1,ncal)) );   

if fwx_sw_out(ii)<0  % If losses exceed the available flux, limit it to zero and adjust bottom-friction dissipation.
    Dsw_out_fc(min(ii+1,ncal)) = Dsw_out_fc(min(ii+1,ncal)) + fwx_sw_out(ii)/dx;  
    fwx_sw_out(ii)=0;    
end

ew_sw_out(ii)=fwx_sw_out(ii)/(cg_sw_out(ii)*cos(angle_sw(ii)) + U);
Hsw_out(ii)=sqrt(ew_sw_out(ii)/shape_sw_out(ii)/rho/grav);
    
U = Usw_in(min(ii+1,ncal))+  Uig_in(min(ii+1,ncal))  +Uc_in(min(ii+1,ncal))  +  0  +  Uig_out(min(ii+1,ncal))  +  Uc_out(min(ii+1,ncal)); 
U = -U *(strcmpi(setup.current_advection,'yes') == 1);  % Reverse the current-advection sign for seaward propagation.
[k_sw_out(ii),c_sw_out(ii),L_sw_out(ii),n_sw_out(ii),kh_sw_out(ii),HL_sw_out(ii),Hh_sw_out(ii),Ursell_sw_out(ii),shape_sw_out(ii),~,...
     ew_sw_out(ii),cg_sw_out(ii),~,sxx_sw_out(ii),qwx_sw_out(ii)] ...
= free_wave_characteristics_functions (Hsw_out(ii),Tsw(ii),angle_sw(ii),depth(ii),mwl(ii),U,...
                                                      setup.c_formula,setup.B_formula,setup.beta_formula);
                                                                            
    
% Total short-wave height (variance-based combination).
% NOTE: this simple expression neglects the phase-correlation between incident and reflected waves !  
% This error may be large under quasi-standing (resonant) conditions (e.g., in front of a strongly reflective structure)
Hsw_tt(ii)=sqrt(Hsw_in(ii)^2+Hsw_out(ii)^2); 
    

% Outgoing IG waves.
if strcmpi(setup.ig_wave,'yes') == 1
    U = Usw_in(max(ii-1,1))  +Uig_in(min(ii+1,ncal)) +  Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  0  +  Uc_out(max(ii-1,1));
    U = -U *(strcmpi(setup.current_advection,'yes') == 1);  % Reverse the current-advection sign for seaward propagation.
    [k_ig_out(ii),c_ig_out(ii),L_ig_out(ii),n_ig_out(ii),kh_ig_out(ii),HL_ig_out(ii),Hh_ig_out(ii),Ursell_ig_out(ii),shape_ig_out(ii),~,...
        ew_ig_out(ii),cg_ig_out(ii),fwx_ig_out(ii),sxx_ig_out(ii),qwx_ig_out(ii)] ...
    = free_wave_characteristics_functions (Hig_out(min(ii+1,ncal)),Tig(min(ii+1,ncal)),angle_ig(ii),depth(ii),mwl(ii),U,...
          setup.c_formula,setup.B_formula,setup.beta_formula);
      
    fwx_ig_out(ii) = fwx_ig_out(min(ii+1,ncal)) ... 
                    + dx*(   Dig_in_igout(ii)    )...
                    -  dx*(   Dig_out_fc(min(ii+1,ncal))  +   Dig_out_igin(min(ii+1,ncal))    );
   
     if fwx_ig_out(ii)<0  % If losses exceed the available flux, limit it to zero and adjust bottom-friction dissipation.
         Dig_out_fc(min(ii+1,ncal)) = Dig_out_fc(min(ii+1,ncal)) + fwx_ig_out(ii)/dx;
        fwx_ig_out(ii)=0;    
     end
    
    ew_ig_out(ii)  = fwx_ig_out(ii)/(cg_ig_out(ii)*cos(angle_ig(ii)) + U); 
    Hig_out(ii) = real(sqrt(ew_ig_out(ii)/shape_ig_out(ii)/rho/grav));
    
    U = Usw_in(max(ii-1,1))  +Uig_in(min(ii+1,ncal)) +  Uc_in(max(ii-1,1))  +  Usw_out(max(ii-1,1))  +  0  +  Uc_out(max(ii-1,1));
    U = -U *(strcmpi(setup.current_advection,'yes') == 1);  % Reverse the current-advection sign for seaward propagation.
    [k_ig_out(ii),c_ig_out(ii),L_ig_out(ii),n_ig_out(ii),kh_ig_out(ii),HL_ig_out(ii),Hh_ig_out(ii),Ursell_ig_out(ii),shape_ig_out(ii),~,...
        ew_ig_out(ii),cg_ig_out(ii),~,sxx_ig_out(ii),qwx_ig_out(ii)] ...
    = free_wave_characteristics_functions (Hig_out(ii),Tig(ii),angle_ig(ii),depth(ii),mwl(ii),U,...
          setup.c_formula,setup.B_formula,setup.beta_formula);     
end


% Total IG-wave height (variance-based combination).
% NOTE: this simple expression neglects the phase-correlation between incident and reflected waves !  
% This error may be large under quasi-standing (resonant) conditions (e.g., in front of a strongly reflective structure)
Hig_tt(ii)=sqrt(Hig_in(ii)^2+Hig_out(ii)^2);


% Currents updated in the seaward sweep.
Usw_out(ii) = - qwx_sw_out(ii)  /  (depth(ii)+mwl(ii));
Ur_out(ii) = 0;
Uig_out(ii) =- qwx_ig_out(ii)  /  (depth(ii)+mwl(ii));
Uc_out(ii) = - (Usw_out(ii) + Ur_out(ii) + Uig_out(ii));
US(ii) = -Uc_in(ii) - Uc_out(ii);
UE(ii) = Uc_in(ii) + Uc_out(ii);
UL(ii) = US(ii) + UE(ii);


% Mean water-level variation.
term1=(sxx_sw_in(min(ii+1,ncal))-sxx_sw_in(ii))/dx;
term2=(rxx_sw_in(min(ii+1,ncal))-rxx_sw_in(ii))/dx;
term3=(sxx_ig_in(min(ii+1,ncal))-sxx_ig_in(ii))/dx;
term4=(rxx_ig_in(min(ii+1,ncal))-rxx_ig_in(ii))/dx;
term5=(sxx_sw_out(min(ii+1,ncal))-sxx_sw_out(ii))/dx;
term6=(sxx_ig_out(min(ii+1,ncal))-sxx_ig_out(ii))/dx;
term7=rho*Cf(ii)*UE(ii)*sqrt((1.16*ubx_sw_in(ii)/sqrt(2))^2+(UE(ii))^2);   
term8=rho*(depth(ii)+mwl(ii))*UL(ii)*(UL(min(ii+1,ncal))-UL(ii))/dx;
mwl(ii)=mwl(min(ii+1,ncal))  +  (term1+term2+term3+term4+term5+term6+term7+term8)*dx/rho/grav/(depth(ii)+mwl(ii));
mwl(ii)=mwl(ii)*(strcmpi(setup.mwl_calculation,'yes') == 1);
mwl_out (ii) = mwl(ii);


% Relax the updated mean water level.
mwl(ii) = (1-str2double(setup.relax))*mwl_out(ii) +str2double(setup.relax)*mwl_in(ii);
mwl(ii)=mwl(ii)*(strcmpi(setup.mwl_calculation,'yes') == 1);

% Currents updated after mean-water-level relaxation.
Usw_out(ii) = - qwx_sw_out(ii)  /  (depth(ii)+mwl(ii));
Ur_out(ii) = 0;
Uig_out(ii) =- qwx_ig_out(ii)  /  (depth(ii)+mwl(ii));
Uc_out(ii) = - (Usw_out(ii) + Ur_out(ii) + Uig_out(ii));
US(ii) = -Uc_in(ii) - Uc_out(ii);
UE(ii) = Uc_in(ii) + Uc_out(ii);
UL(ii) = US(ii) + UE(ii);


% Outgoing short-wave energy dissipation.
[Dsw_out_fc(ii),ubx_sw_out(ii),usx_sw_out(ii)] = wave_friction_dissipation_functions (Hsw_out(ii),Tsw(ii),dd50(ii),kh_sw_out(ii),angle_sw(ii));

 if strcmpi(setup.reflected_wave,'yes') == 1     
     [~,c1,~,~,~,~,~,~,~,~,~,~,~,~,~] ...
= free_wave_characteristics_functions (Hsw_out(max(ii-1,1)),Tsw(max(ii-1,1)),angle_sw(max(ii-1,1)),depth(max(ii-1,1)),mwl(max(ii-1,1)),U,...
                                                      setup.c_formula,setup.B_formula,setup.beta_formula);
     [Dsw_out_swin(ii),Kr_sw_outin(ii)] =  wave_reflection_loss_functions ...
          (-1,'bottom',setup.Kr_formula,fwx_sw_out,Hsw_out,Tsw,c_sw_out,depth,mwl, be(:,1),ii,ncal,dxmean,p,c1);  % Shoreward primary wave: +1; seaward primary wave: -1.

      Dsw_out_swin(ii) = min (Dsw_out_swin(ii), fwx_sw_out(ii)/dx);
 end
 
    
% Outgoing IG-wave energy dissipation.
if strcmpi(setup.ig_wave,'yes') == 1
    [Dig_out_fc(ii),ubx_ig_out(ii),usx_ig_out(ii)] = wave_friction_dissipation_functions (Hig_out(ii),Tig(ii),dd50(ii),kh_ig_out(ii),angle_ig(ii));  % or current friction dissipation?         
%   [Dig_out_fc(ii),ubx_ig_out(ii),usx_ig_out(ii)] = wave_friction_dissipation_functions (Hig_out(ii),Tig(ii),Cf(ii),kh_ig_out(ii),angle_ig(ii));  % or current friction dissipation?
%    Dig_out_fc(ii) = rho*Cf(ii)* power( ubx_ig_out(ii) / sqrt(2), 2) * abs(ubx_ig_out(ii) / sqrt(2));   %  corrected expression for long-wave dissipation due to bottom friction
    
    if strcmpi(setup.reflected_wave,'yes') == 1
        [~,c1,~,~,~,~,~,~,~,~,~,~,~,~,~] ...
    = free_wave_characteristics_functions (Hig_out(max(ii-1,1)),Tig(max(ii-1,1)),angle_ig(max(ii-1,1)),depth(max(ii-1,1)),mwl(max(ii-1,1)),U,...
          setup.c_formula,setup.B_formula,setup.beta_formula);     
        [Dig_out_igin(ii),Kr_ig_outin(ii)] = wave_reflection_loss_functions ...
           (-1,'bottom',setup.Kr_formula,fwx_ig_out,Hig_out,Tig,c_ig_out,depth,mwl, be(:,1),ii,ncal,dxmean,p,c1);  % Shoreward primary wave: +1; seaward primary wave: -1.

       Dig_out_igin(ii) = min (Dig_out_igin(ii), fwx_ig_out(ii)/dx);
    end                                     
end         


end
end

% Apply a mass-conservation/reference-level correction to mean water level.
if strcmpi(setup.mwl_correction,'yes') == 1
    mwl(1:ncal) = mwl(1:ncal) - mean(mwl(1:ncal));   
else
    mwl(1:ncal) = mwl(1:ncal) - (mwl(1)-bc(index.mwl)) - 1e-19;   
end

% Check iterative convergence.
mwl_final(1:ncal,iter) = mwl(1:ncal);
if iter >1
    Errors (iter-1) = mean(abs((mwl_final(1:ncal,iter)-mwl_final(1:ncal,iter-1))./(0.5*mwl_final(1:ncal,iter)+0.5*mwl_final(1:ncal,iter-1)))); 
    if (Errors (iter-1)<=str2double(setup.converg))
        break;
    end
end

% Reset the incoming-wave boundary heights.
if (bc(index.Hsw_tt)~=0);   Hsw_in(1) = Hsw_in(1) * bc(index.Hsw_tt) /Hsw_tt(1) *0.7 +  Hsw_in(1)*0.3;   end
if (bc(index.Hig_tt)~=0);    Hig_in(1) = Hig_in(1) * bc(index.Hig_tt) /Hig_tt(1) *0.7 +  Hig_in(1)*0.3;       end
if (bc(index.Hsw_tt)==0 && bc(index.Hsw_in)~=0);   Hsw_in(1) = bc(index.Hsw_in);   end
if (bc(index.Hig_tt)==0 && bc(index.Hig_in)~=0);    Hig_in(1) = bc(index.Hig_in);      end

end

%% Model output (additional variables may be added)
output.nowtime      =       nowtime          ;
output.ncal         =       ncal             ;
output.Errors         =       Errors             ;
output.swl          =       swl              ;
output.mwl          =       mwl              ;
output.mwl_in          =       mwl_in              ;
output.mwl_out          =       mwl_out              ;
output.be           =       be               ;
output.depth        =       depth            ;
output.p            =       p                ;
output.dd50         =       dd50             ;
output.rhor         =       rhor             ;
output.Cf           =       Cf               ;
output.s0           =       s0               ;

output.UE           =       UE          ; 
output.US           =       US          ; 
output.UL           =       UL          ; 
output.Uy                 =       Uy                ; 

output.Hsw_tt             =       Hsw_tt            ; 
output.Tsw                =       Tsw               ; 
output.angle_sw           =       angle_sw          ; 
output.Hsw_in             =       Hsw_in            ; 
output.k_sw_in            =       k_sw_in           ; 
output.L_sw_in            =       L_sw_in           ; 
output.kh_sw_in           =       kh_sw_in          ; 
output.c_sw_in            =       c_sw_in           ; 
output.n_sw_in            =       n_sw_in           ; 
output.cg_sw_in           =       cg_sw_in          ; 
output.Hh_sw_in           =       Hh_sw_in          ; 
output.HL_sw_in           =       HL_sw_in          ; 
output.Ursell_sw_in       =       Ursell_sw_in      ; 
output.shape_sw_in        =       shape_sw_in       ; 
output.ew_sw_in           =       ew_sw_in          ; 
output.gamma_sw_in        =       gamma_sw_in       ; 
output.ubx_sw_in          =       ubx_sw_in         ; 
output.usx_sw_in          =       usx_sw_in         ; 
output.fwx_sw_in          =       fwx_sw_in         ; 
output.sxx_sw_in          =       sxx_sw_in         ; 
output.qwx_sw_in          =       qwx_sw_in         ; 
output.Kr_sw_inout      =    Kr_sw_inout ;
output.Kr_sw_outin      =    Kr_sw_outin ;
output.Hsw_out            =       Hsw_out           ; 
output.k_sw_out           =       k_sw_out          ; 
output.L_sw_out           =       L_sw_out          ; 
output.kh_sw_out          =       kh_sw_out         ; 
output.c_sw_out           =       c_sw_out          ; 
output.n_sw_out           =       n_sw_out          ; 
output.cg_sw_out          =       cg_sw_out         ; 
output.Hh_sw_out          =       Hh_sw_out         ; 
output.HL_sw_out          =       HL_sw_out         ; 
output.Ursell_sw_out      =       Ursell_sw_out     ; 
output.shape_sw_out       =       shape_sw_out      ; 
output.ew_sw_out          =       ew_sw_out         ; 
output.gamma_sw_out       =       gamma_sw_out      ; 
output.ubx_sw_out         =       ubx_sw_out        ; 
output.usx_sw_out         =       usx_sw_out        ; 
output.fwx_sw_out         =       fwx_sw_out        ; 
output.sxx_sw_out         =       sxx_sw_out        ; 
output.qwx_sw_out         =       qwx_sw_out        ; 

output.Hig_tt             =       Hig_tt            ; 
output.Tig                =       Tig               ; 
output.angle_ig           =       angle_ig          ; 
output.Hig_in             =       Hig_in            ; 
output.k_ig_in            =       k_ig_in           ; 
output.L_ig_in            =       L_ig_in           ; 
output.kh_ig_in           =       kh_ig_in          ; 
output.c_ig_in            =       c_ig_in           ; 
output.n_ig_in            =       n_ig_in           ; 
output.cg_ig_in           =       cg_ig_in          ; 
output.Hh_ig_in           =       Hh_ig_in          ; 
output.HL_ig_in           =       HL_ig_in          ; 
output.Ursell_ig_in       =       Ursell_ig_in      ; 
output.shape_ig_in        =       shape_ig_in       ; 
output.ew_ig_in           =       ew_ig_in          ; 
output.gamma_ig_in        =       gamma_ig_in       ; 
output.ubx_ig_in          =       ubx_ig_in         ; 
output.usx_ig_in          =       usx_ig_in         ; 
output.fwx_ig_in          =       fwx_ig_in         ; 
output.sxx_ig_in          =       sxx_ig_in         ; 
output.qwx_ig_in          =       qwx_ig_in         ; 
output.phi_ig_in         =       phi_ig_in        ;
output.Kr_ig_inout      =    Kr_ig_inout ;
output.Kr_ig_outin      =    Kr_ig_outin ;
output.Hig_out            =       Hig_out           ; 
output.k_ig_out           =       k_ig_out          ; 
output.L_ig_out           =       L_ig_out          ; 
output.kh_ig_out          =       kh_ig_out         ; 
output.c_ig_out           =       c_ig_out          ; 
output.n_ig_out           =       n_ig_out          ; 
output.cg_ig_out          =       cg_ig_out         ; 
output.Hh_ig_out          =       Hh_ig_out         ; 
output.HL_ig_out          =       HL_ig_out         ; 
output.Ursell_ig_out      =       Ursell_ig_out     ; 
output.shape_ig_out       =       shape_ig_out      ; 
output.ew_ig_out          =       ew_ig_out         ; 
output.gamma_ig_out       =       gamma_ig_out      ; 
output.ubx_ig_out         =       ubx_ig_out        ; 
output.usx_ig_out         =       usx_ig_out        ; 
output.fwx_ig_out         =       fwx_ig_out        ; 
output.sxx_ig_out         =       sxx_ig_out        ; 
output.qwx_ig_out         =       qwx_ig_out        ; 
output.phi_ig_out         =       phi_ig_out        ;

output.rl_sw_in           =       rl_sw_in          ; 
output.er_sw_in           =       er_sw_in          ; 
output.ra_sw_in           =       ra_sw_in          ; 
output.rh_sw_in           =       rh_sw_in          ; 
output.beta_sw_in         =       beta_sw_in        ; 
output.frx_sw_in          =       frx_sw_in         ; 
output.rxx_sw_in          =       rxx_sw_in         ; 
output.qrx_sw_in          =       qrx_sw_in         ; 

output.rl_ig_in           =       rl_ig_in          ; 
output.er_ig_in           =       er_ig_in          ; 
output.ra_ig_in           =       ra_ig_in          ; 
output.rh_ig_in           =       rh_ig_in          ; 
output.beta_ig_in         =       beta_ig_in        ; 
output.frx_ig_in          =       frx_ig_in         ; 
output.rxx_ig_in          =       rxx_ig_in         ; 
output.qrx_ig_in          =       qrx_ig_in         ; 

% Current components
output.Usw_in   = Usw_in;
output.Ur_in    = Ur_in;
output.Uig_in   = Uig_in;
output.Uc_in    = Uc_in;

output.Usw_out  = Usw_out;
output.Ur_out   = Ur_out;
output.Uig_out  = Uig_out;
output.Uc_out   = Uc_out;

output.US       = US;
output.UE       = UE;
output.UL       = UL;


output.Dsw_in_br          =       Dsw_in_br         ; 
output.Dsw_in_fc          =       Dsw_in_fc         ; 
output.Dsw_in_swout       =       Dsw_in_swout      ; 
output.Dsw_in_LFin        =       Dsw_in_LFin       ; 
output.Dsw_out_br         =       Dsw_out_br        ; 
output.Dsw_out_fc         =       Dsw_out_fc        ; 
output.Dsw_out_swin       =       Dsw_out_swin      ; 
output.Dig_in_br          =       Dig_in_br         ; 
output.Dig_in_fc          =       Dig_in_fc         ; 
output.Dig_in_igout       =       Dig_in_igout      ; 
output.Dig_out_br         =       Dig_out_br        ; 
output.Dig_out_fc         =       Dig_out_fc        ; 
output.Dig_out_igin       =       Dig_out_igin      ; 
output.Drsw_in_turb       =       Drsw_in_turb      ; 
output.Drig_in_turb       =       Drig_in_turb      ; 
output.Drsw_in_LFin        =       Drsw_in_LFin       ;
output.Dmwl = Dmwl;