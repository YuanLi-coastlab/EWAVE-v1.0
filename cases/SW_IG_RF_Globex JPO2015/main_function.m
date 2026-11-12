clearvars;
close all;

%% Locate the case directory and EWAVE source files
case_dir = fileparts(mfilename('fullpath'));
if isempty(case_dir)
    error('EWAVE:ScriptLocationUnknown','Run the saved main_function.m script.');
end

% The project root contains src/ alongside case/ or cases/.
project_dir = fileparts(case_dir);
while exist(fullfile(project_dir,'src','EWAVE_functions_2020d.m'),'file') ~= 2
    parent_dir = fileparts(project_dir);
    if isempty(parent_dir) || strcmp(parent_dir,project_dir)
        error('EWAVE:SourceNotFound', ...
            'Could not find src/EWAVE_functions_2020d.m above "%s".',case_dir);
    end
    project_dir = parent_dir;
end
src_dir = fullfile(project_dir,'src');
addpath(genpath(src_dir),'-begin');

%% Model inputs
casename = 'A2';
% Bathymetry points may be unevenly spaced. The first point is the offshore boundary.
input.be_file  = fullfile(case_dir,['be_' casename '.inp']);  % Bathymetry
input.bc_file  = fullfile(case_dir,['bc_' casename '.inp']);  % Boundary conditions
input.obs_file = fullfile(case_dir,['obs_' casename '.inp']); % Observations
input.index.t=1;  input.index.x=2;  input.index.swl=3;  input.index.mwl=4;  input.index.spectype=5;
input.index.Hsw_tt =  6;   input.index.Hsw_in =  7;   input.index.Hsw_out =  8;  input.index.Tsw = 9;  input.index.angle_sw = 10;
input.index.Hig_tt  =11;   input.index.Hig_in  = 12;  input.index.Hig_out   =13;  input.index.Tig  =14; input.index.angle_ig =15;   
input.index.Lr_sw =16;   input.index.Lr_ig =17;    input.index.V_wind =18;

%% Model settings
setup.wave_type = 'random';       % Wave type ('random' or 'regular').
setup.B_formula = '0.125';        % Wave energy shape factor ('0.125' or 'Z2018').
setup.c_formula = 'nonlinear';   % Wave celerity formula ('linear' or 'nonlinear').
setup.Db_formula = 'B1998';       % Breaking dissipation (BJ1978, TG1983, B1998, JB2007, BTG2020).
setup.gamma_formula = 'zl2020';   % Breaker index parameterization.
% Mad1976 SH1985 NE1987 SH1989 N1990 Lip1996 RS2000 T2001 TM2002 G2010
% BJ1978 TG1983 BS1985 R2003 Sal2015 Lin2017 zl2020
% A2008_JB2007 A2008_TG1983 A2008_B1998 A2008_BJ1978
setup.reflected_wave = 'yes';     % Include wave reflection ('yes' or 'no').
% Reflection options are listed in order.
% 1. Bed reflection formula ('LL2005', 'M1983', or 'L1932').
% 2. Swash reflection formula ('ZV2008', 'M2017', 'V2005', 'B1974', or 'YT2011').
% 3. Use the 'original' or 'corrected' formulation.
% 4. Search distance for the effective slope toe, in wavelengths.
setup.Kr_formula = {'M1983','ZV2008','corrected','1'};
setup.ig_wave = 'yes';            % Calculate infragravity (IG) waves ('yes' or 'no').
setup.mwl_calculation = 'yes';    % Calculate wave setup and setdown ('yes' or 'no').
setup.surface_roller = 'yes';     % Include surface rollers ('yes' or 'no').
setup.beta_formula = 'Z2017';     % Roller-front slope ('0.1' or 'Z2017').
% Fractions of wave-induced water-mass exchange (WME) allocated to IG waves.
setup.roller_IG_transfer = '0.7'; % Roller contribution (0 to 1).
setup.SW_IG_transfer = '0.3';     % Short-wave contribution (0 to 1).

setup.rhor = '1000';              % Roller density [kg/m^3].
setup.s0_method = 'offshore';     % Reference for wave steepness ('deepwater' or 'offshore').
setup.d50 = '0.0';                % Bed roughness scale for wave friction [m].
setup.Cf = '0.0';                 % Mean-current bottom friction coefficient for bed shear stress.
setup.Tig = '9';                  % IG-to-short-wave peak-period ratio used if the IG period is missing.
setup.dx = '0.1';                 % Computational grid spacing [m].
setup.kh_cut = '0.1';             % Minimum relative water depth (kh) at the shoreward grid limit.
setup.relax = '0';                % Relaxation factor for mean water level iterations.
setup.converg = '0.001';          % Relative water level convergence tolerance when reflection is enabled.
setup.N_iteration = '20';         % Maximum number of iterations.
setup.mwl_correction = 'no';      % Correct mean water level based on mass conservation.
setup.current_advection = 'no';   % Include the effect of currents on wave energy flux ('yes' or 'no').
setup.igbc_LHS62 = {'yes' , '1'}; % Use the LHS62 solution for offshore incident IG height if unspecified.


%% Read and preprocess inputs
input.be   = importdata (input.be_file);    input.be   =  input.be.data;
input.bc   = importdata (input.bc_file);    input.bc   =  input.bc.data;
input.obs = importdata (input.obs_file);  input.obs  = input.obs.data;
% Interpolate bathymetry onto the uniform computational grid.
x  = input.be(1,1):str2double(setup.dx):input.be(end,1); 
zb = interp1(input.be(:,1),input.be(:,2),x);
x=x';zb=zb'; 
input.be = []; input.be(:,1)=x; input.be(:,2)=zb; 
% Convert wave heights to RMS values (spectype 1 for RMS and 2 for significant).
ivec = [input.index.Hsw_tt   input.index.Hsw_in   input.index.Hsw_out ...
            input.index.Hig_tt    input.index.Hig_in    input.index.Hig_out ];
input.bc(:,ivec)   = input.bc(:,ivec)./    [sqrt(input.bc(:,input.index.spectype)) sqrt(input.bc(:,input.index.spectype)) ...
                                                            sqrt(input.bc(:,input.index.spectype)) sqrt(input.bc(:,input.index.spectype)) ...
                                                            sqrt(input.bc(:,input.index.spectype)) sqrt(input.bc(:,input.index.spectype))];
input.obs(:,ivec) = input.obs(:,ivec)./  [sqrt(input.obs(:,input.index.spectype)) sqrt(input.obs(:,input.index.spectype)) ...
                                                            sqrt(input.obs(:,input.index.spectype)) sqrt(input.obs(:,input.index.spectype)) ...
                                                            sqrt(input.obs(:,input.index.spectype)) sqrt(input.obs(:,input.index.spectype))];
input.bc(:,input.index.spectype)=1;
input.obs(:,input.index.spectype)=1;
clear ivec;
% Read the time records and count the gauges per record.
time=input.bc(:,1);
ngauge=size(input.obs,1)/length(time); % This count should be an integer.
% Fill missing boundary condition values.
for it=1:length(time)
    if isnan(input.bc(it,input.index.Tig));  input.bc(it,input.index.Tig) = str2double(setup.Tig) * input.bc(it,input.index.Tsw); end
    if isnan(input.bc(it,input.index.Hsw_in));   input.bc(it,input.index.Hsw_in) = input.bc(it,input.index.Hsw_tt);   end
    if isnan(input.bc(it,input.index.Hig_in));   input.bc(it,input.index.Hig_in) = input.bc(it,input.index.Hig_tt);   end
    input.bc(it, isnan(input.bc(it,:)))=0; 
end

%% Run EWAVE
for it=1:length(time)
    [output(it)] = EWAVE_functions_2020d (setup, input.be, input.bc(it,:),input.index);
end





%% SINGLE-PROFILE VALIDATION: RMSE and Relative bias
% RMSE and relative bias as reported by Leijnse et al. (2024), Table 1.
% One instant, cross-shore profile. No minimum water-depth threshold.
% Relative bias (%) = 100*mean(M-O)/mean(O).

if numel(time) ~= 1
    error('This validation block requires exactly one time record.');
end

it = 1;
itvec = 1:ngauge;
xgauge = input.obs(itvec,input.index.x);
xgauge = xgauge(:);
ncal = output(it).ncal;
if ncal < 2
    error('Need at least two valid model nodes for interpolation.');
end
xmodel = input.be(1:ncal,1);
xmodel = xmodel(:);

% 8 columns: observed/modeled total SW, total IG, incoming IG, outgoing IG.
stat_vars = {'Hsw_tt','Hig_tt','Hig_in','Hig_out'};
column_names = {'Hsw_obs','Hsw_model', ...
                'Hig_obs','Hig_model', ...
                'Hig_in_obs','Hig_in_model', ...
                'Hig_out_obs','Hig_out_model'};
comparison_data = nan(ngauge,8);

RMSE = nan(4,1);             % Root-mean-square error [m]
RelativeBias = nan(4,1);     % Signed relative bias [%]
N = zeros(4,1);              % Number of valid observation-model pairs


for iv = 1:numel(stat_vars)

    varname = stat_vars{iv};

    % Observations
    O = input.obs(itvec,input.index.(varname));
    O = O(:);

    % Model results
    Hmodel = output(it).(varname);
    Hmodel = Hmodel(1:ncal);
    Hmodel = Hmodel(:);

    % Sort model coordinates for interpolation
    [x_sorted,idx] = sort(xmodel);
    H_sorted = Hmodel(idx);

    % Check duplicate model coordinates
    if any(diff(x_sorted) == 0)
        error('Duplicate model coordinates detected.');
    end

    % Interpolate model results to observation locations
    M = interp1(x_sorted,H_sorted,xgauge,'linear',NaN);
    M = M(:);

    % Store original paired values
    comparison_data(:,2*iv-1) = O;
    comparison_data(:,2*iv)   = M;

    % Check physically invalid model predictions
    n_negative = sum(isfinite(M) & M < 0);

    if n_negative > 0
        warning('%s: %d negative modeled wave heights.', ...
                varname,n_negative);
    end

    % Valid observation-model pairs
    valid = isfinite(O) & isfinite(M) & ...
            O >= 0 & M >= 0;

    O_valid = O(valid);
    M_valid = M(valid);

    N(iv) = numel(O_valid);

    if N(iv) == 0
        continue
    end

    % Error statistics
    err = M_valid - O_valid;

    RMSE(iv) = sqrt(mean(err.^2));

    if sum(O_valid) > 0
        RelativeBias(iv) = ...
            100 * sum(err) / sum(O_valid);
    end

end


validation_summary = table(stat_vars(:),N,RMSE,RelativeBias, ...
    'VariableNames',{'Variable','N','RMSE_m','RelativeBias_percent'});

disp(' ');
disp(['SINGLE-PROFILE VALIDATION: ' casename]);
disp('Metrics: RMSE (m), Relative bias (%)');
disp(validation_summary);

% Same paired-data MAT output used by the previous validation scripts.
row_gauge_x_m = xgauge;
row_time_id = repmat(time(it),ngauge,1);
mat_file = fullfile(case_dir,['validation_' casename '_paired.mat']);
save(mat_file,'comparison_data','column_names','row_gauge_x_m', ...
    'row_time_id','validation_summary','stat_vars','-v7');
fprintf('Saved %s: %d gauge rows x 8 columns.\n', ...
    mat_file,size(comparison_data,1));



%% Plot model results
for it=1:length(time)
    itvec=(it-1)*ngauge+1:it*ngauge;
    swl=0;
    ncal=output(it).ncal;
    xgauge=input.obs(itvec,input.index.x);
    
    figure(it);
    
    set(gcf,'position',[100,100,1100,500])
        
        subplot(3,3,[1 4]); hold on; 
            xlim([-inf,inf]); %ylim([-0.5,4]);
            xlabel('{\itx} (m)','FontSize',8,'Fontname','Helvetica');
            ylabel('{\itz} (m)','Fontsize',8,'Fontname','Helvetica');
            set(gca,'FontSize',8,'Fontname','Helvetica');box('on');
                
            t1=title ({[casename,...
                ':   \ittime \rm= ', num2str(time(it)) ,...
                '    \itwavetype \rm= ', setup.wave_type,...
                ':   \itDb \rm= ', setup.Db_formula ,...
                '    \gamma = ', setup.gamma_formula,...
                '    \itB \rm= ', setup.B_formula,...
                '    \itc \rm= ', setup.c_formula];...
                 ['    \beta = ', setup.beta_formula,...
                '    \itroller \rm= ', setup.surface_roller...
               '    \itIGwave \rm= ', setup.ig_wave,...
                '    \itreflection \rm= ', setup.reflected_wave,...
                '    \itKr \rm= ',...
                ]},'Fontsize',8,'Fontname','Helvetica','color','b'); 
         
            plot(xgauge,input.obs(itvec,input.index.Hsw_tt)+swl,...
                'o','color',[0 .45 .74],'Markersize',5,'LineWidth',0.5);
         f1=plot(x(1:ncal),output(it).Hsw_tt(1:ncal)+swl,'-','color',[0 .45 .74],'LineWidth',1);

            plot(xgauge,input.obs(itvec,input.index.Hig_tt)+swl,...
                'o','color',[0.87 .49 0],'Markersize',5,'LineWidth',0.5);
         f2=plot(x(1:ncal),output(it).Hig_tt(1:ncal)+swl,'-','color',[0.87 .49 0],'LineWidth',1);

            plot(xgauge,input.obs(itvec,input.index.mwl)+swl,...
                'o','color','k','Markersize',5,'LineWidth',0.5);
         f3=plot(x(1:ncal),output(it).mwl(1:ncal)+swl,'-','color',[0 .65 0],'LineWidth',1);
         
         plot(xgauge,input.obs(itvec,input.index.angle_sw)+swl,...
                'o','color','b','Markersize',5,'LineWidth',0.5);
         f4=plot(x(1:ncal),output(it).angle_sw(1:ncal)+swl,'-','color','b','LineWidth',1);
         
         f5=plot(x(1:ncal),output(it).UE(1:ncal)+swl,'--','color',[0 0 0],'LineWidth',1);
          
            % plot(input.be(:,1),input.be(:,2),'-','color',[.5 .5 .5],'LineWidth',0.5);
            % plot([x(1) x(ncal)], [swl swl],'--','color',[.5 .5 .5],'LineWidth',0.5);
            
            legend([f1,f2,f3,f4,f5],'Hsw_T','Hig_T','MWL','\beta','Undertow','location','best');
            legend('boxoff');

        subplot(3,3,7); hold on; 
       xlim([-inf,inf]); %ylim([min(zb)-0.1,max(zb)+0.5]); 
             swl0=input.bc(1,3); fu=20;
            xlabel('{\itx} (m)','FontSize',8,'Fontname','Helvetica');
            ylabel('{\itz} (m)','Fontsize',8,'Fontname','Helvetica');
            set(gca,'FontSize',8,'Fontname','Helvetica');box('on');    
            plot(input.be(:,1),input.be(:,2),'-','color',[.5 .5 .5],'LineWidth',0.5);
            plot([x(1) x(ncal)], [swl0 swl0],'--','color',[.5 .5 .5],'LineWidth',0.5);
% load obs_U
% plot(obs.Uloc(1)+fu*obs.U1(:,1)*sqrt(9.81*interp1(x,output.depth+output.mwl,obs.Uloc(1))),interp1(input.be(:,1),input.be(:,2),obs.Uloc(1))+ interp1(x,output.depth+output.mwl,obs.Uloc(1))*obs.U1(:,2),'o','color','k','Markersize',3,'LineWidth',0.5);
% plot(obs.Uloc(2)+fu*obs.U2(:,1)*sqrt(9.81*interp1(x,output.depth+output.mwl,obs.Uloc(2))),interp1(input.be(:,1),input.be(:,2),obs.Uloc(2))+ interp1(x,output.depth+output.mwl,obs.Uloc(2))*obs.U2(:,2),'o','color','k','Markersize',3,'LineWidth',0.5);
% plot(obs.Uloc(3)+fu*obs.U3(:,1)*sqrt(9.81*interp1(x,output.depth+output.mwl,obs.Uloc(3))),interp1(input.be(:,1),input.be(:,2),obs.Uloc(3))+ interp1(x,output.depth+output.mwl,obs.Uloc(3))*obs.U3(:,2),'o','color','k','Markersize',3,'LineWidth',0.5);
% plot(obs.Uloc(4)+fu*obs.U4(:,1)*sqrt(9.81*interp1(x,output.depth+output.mwl,obs.Uloc(4))),interp1(input.be(:,1),input.be(:,2),obs.Uloc(4))+ interp1(x,output.depth+output.mwl,obs.Uloc(4))*obs.U4(:,2),'o','color','k','Markersize',3,'LineWidth',0.5);
% plot([obs.Uloc(1) obs.Uloc(1)],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(1)) swl0],'--','color',[.5 .5 .5],'LineWidth',0.5);
% plot([obs.Uloc(2) obs.Uloc(2)],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(2)) swl0],'--','color',[.5 .5 .5],'LineWidth',0.5);
% plot([obs.Uloc(3) obs.Uloc(3)],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(3)) swl0],'--','color',[.5 .5 .5],'LineWidth',0.5);
% plot([obs.Uloc(4) obs.Uloc(4)],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(4)) swl0],'--','color',[.5 .5 .5],'LineWidth',0.5);
% 
% plot([obs.Uloc(1)+fu*output.UE(where_nearest(x,obs.Uloc(1))) obs.Uloc(1)+fu*output.UE(where_nearest(x,obs.Uloc(1)))],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(1)) swl0-0.5*output.Hsw_in(where_nearest(x,obs.Uloc(1)))],'b-','LineWidth',0.5);
% plot([obs.Uloc(2)+fu*output.UE(where_nearest(x,obs.Uloc(2))) obs.Uloc(2)+fu*output.UE(where_nearest(x,obs.Uloc(2)))],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(2)) swl0-0.5*output.Hsw_in(where_nearest(x,obs.Uloc(2)))],'b-','LineWidth',0.5);
% plot([obs.Uloc(3)+fu*output.UE(where_nearest(x,obs.Uloc(3))) obs.Uloc(3)+fu*output.UE(where_nearest(x,obs.Uloc(3)))],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(3)) swl0-0.5*output.Hsw_in(where_nearest(x,obs.Uloc(3)))],'b-','LineWidth',0.5);
% plot([obs.Uloc(4)+fu*output.UE(where_nearest(x,obs.Uloc(4))) obs.Uloc(4)+fu*output.UE(where_nearest(x,obs.Uloc(4)))],[interp1(input.be(:,1),input.be(:,2),obs.Uloc(4)) swl0-0.5*output.Hsw_in(where_nearest(x,obs.Uloc(4)))],'b-','LineWidth',0.5);


        subplot(2,3,2); hold on;
%             xlim([-1000,0]); ylim([-0.5,4]);
            xlabel('{\itx} (m)','FontSize',8,'Fontname','Helvetica');
            ylabel('{Hsw} (m)','Fontsize',8,'Fontname','Helvetica');
            set(gca,'FontSize',8,'Fontname','Helvetica');box('on');  

            f1=plot(xgauge(1:ngauge),input.obs(itvec(1:ngauge),input.index.Hsw_in)+swl,...
                '>','color',[0 .65 0],'Markersize',5,'LineWidth',0.5,'Markerfacecolor','none');
            plot(x(1:ncal),output(it).Hsw_in(1:ncal)+swl,'-','color',[0 .45 .74],'LineWidth',1);

            f2=plot(xgauge(1:ngauge),input.obs(itvec(1:ngauge),input.index.Hsw_out)+swl,...
                '<','color',[1 .15 .74],'Markersize',5,'LineWidth',0.5,'Markerfacecolor','none');
            plot(x(1:ncal),output(it).Hsw_out(1:ncal)+swl,'-','color',[0 .45 .74],'LineWidth',1);

            f3=plot(x(1:ncal),output(it).Hsw_tt(1:ncal)+swl,'-.','color',[.5 .5 .5],'LineWidth',0.5);


%             plot(input.be(:,1),input.be(:,2),'-','color',[.5 .5 .5],'LineWidth',0.5);
%             plot([x(1) x(ncal)], [swl swl],'--','color',[.5 .5 .5],'LineWidth',0.5);

            legend([f1 f2 f3],'Hsw_{in}','Hsw_{out}', 'Hsw_{T}','location','best');
            legend('boxoff');


        subplot(2,3,5); hold on;
%             xlim([-1000,0]); ylim([-0.1,0.5]);
            xlabel('{\itx} (m)','FontSize',8,'Fontname','Helvetica');
            ylabel('{Hig} (m)','Fontsize',8,'Fontname','Helvetica');
            set(gca,'FontSize',8,'Fontname','Helvetica');box('on');  

            f1=plot(xgauge,input.obs(itvec,input.index.Hig_in)+swl,...
                '>','color',[0 .65 0],'Markersize',5,'LineWidth',0.5,'Markerfacecolor','none');
            plot(x(1:ncal),output(it).Hig_in(1:ncal)+swl,'-','color',[0 .65 0],'LineWidth',1);

            f2=plot(xgauge,input.obs(itvec,input.index.Hig_out)+swl,...
                '<','color',[1 .15 .74],'Markersize',5,'LineWidth',0.5,'Markerfacecolor','none');
            plot(x(1:ncal),output(it).Hig_out(1:ncal)+swl,'-','color',[1 .15 .74],'LineWidth',1);

            f3=plot(x(1:ncal),output(it).Hig_tt(1:ncal)+swl,'-.','color',[.5 .5 .5],'LineWidth',0.5);

%             plot(input.be(:,1),input.be(:,2),'-','color',[.5 .5 .5],'LineWidth',0.5);
%             plot([x(1) x(ncal)], [swl swl],'--','color',[.5 .5 .5],'LineWidth',0.5);

            legend([f1 f2 f3],'Hig_{in}','Hig_{out}','Hig_{T}','location','best');
            legend('boxoff');


         subplot(2,3,3); hold on;
%              xlim([-1000,0]); ylim([-0.5,4]);
            xlabel('{\itx} (m)','FontSize',8,'Fontname','Helvetica');
            ylabel('{velocity} (m/s)','Fontsize',8,'Fontname','Helvetica');
            set(gca,'FontSize',8,'Fontname','Helvetica');box('on');  

            f1=plot(x(1:ncal),output(it).ubx_sw_in(1:ncal)+swl,'-','color',[0 .45 .74],'LineWidth',1);
            f2=plot(x(1:ncal),-output(it).ubx_sw_out(1:ncal)+swl,'-','color',[0.87 .49 0],'LineWidth',1);
            f3=plot(x(1:ncal),output(it).ubx_ig_in(1:ncal)+swl,'--','color',[0 .45 .74],'LineWidth',1);
            f4=plot(x(1:ncal),-output(it).ubx_ig_out(1:ncal)+swl,'--','color',[0.87 .49 0],'LineWidth',1);
            f5=plot(x(1:ncal),output(it).ubx_sw_in(1:ncal)-output(it).ubx_sw_out(1:ncal)+output(it).ubx_ig_in(1:ncal)-output(it).ubx_ig_out(1:ncal)...
                +swl,'-','color',[0 0 0],'LineWidth',2);

%             plot(input.be(:,1),input.be(:,2),'-','color',[.5 .5 .5],'LineWidth',0.5);
%             plot([x(1) x(ncal)], [swl swl],'--','color',[.5 .5 .5],'LineWidth',0.5);

            legend([f1 f2 f3 f4 f5],'ub_{sw}+','ub_{sw}-','ub_{ig}+','ub_{ig}-','ub_{all}','location','southwest');
            legend('boxoff');


         subplot(2,3,6); hold on;
%             xlim([-80,0]); ylim([-0.1,0.2]);
            xlabel('{Number of iterations}','FontSize',8,'Fontname','Helvetica');
            ylabel('RMSPE','Fontsize',8,'Fontname','Helvetica');
            set(gca,'FontSize',8,'Fontname','Helvetica');box('on');  
            f1=plot(output.Errors,'k-o','LineWidth',1);



end
