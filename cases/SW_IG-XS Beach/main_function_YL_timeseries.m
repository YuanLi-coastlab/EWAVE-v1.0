clearvars;
close all;

%% Locate the case directory and EWAVE source files
case_dir = fileparts(mfilename('fullpath'));
if isempty(case_dir)
    error('EWAVE:ScriptLocationUnknown','Run the saved main_function_YL_timeseries.m script.');
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
% Keep plot_EWAVE_timeseries.m alongside this script or on the MATLAB path.
addpath(case_dir,'-begin');
addpath(genpath(src_dir),'-begin');

%% Model inputs
casename = 'xsw';
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
setup.d50 = '0.0003';                % Bed roughness scale for wave friction [m].
setup.Cf = '0.01';                 % Mean-current bottom friction coefficient for bed shear stress.
setup.Tig = '4';                  % IG-to-short-wave peak-period ratio used if the IG period is missing.
setup.dx = '4';                 % Computational grid spacing [m].
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
% Read the boundary-condition time IDs.
time=input.bc(:,1);
% The fixed gauge positions and time IDs are matched inside the plotting function.
% No assumption about 'ngauge consecutive records per time' is needed.
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


%% TIME-SERIES VALIDATION: TOTAL IG WAVE HEIGHT ONLY
% RMSE [m] and signed relative bias [%] for paired gauge-time records.
% The existing model configuration above is unchanged (including reflection).
% Wave heights in input.obs were converted to RMS conventions above.
% All IG gauge time series share one figure; the bottom panel shows
% the beach profile and every gauge position.

validation_options = struct();
validation_options.minutesPerStep = 12.8; % Used for time-axis fallback ONLY
validation_options.gaugeNumbers = 1:4; % Plot and assess Gauge 1-4 only; exclude Gauge 5
validation_options.minDepth = 0.5;       % Exclude depth < 0.5 m in plot and statistics
validation_options.showTideShade = true;  % Grey transparent tidal-level shading
validation_options.tideShadeAlpha = 0.17; % Transparency strength (0 to 1)
validation_options.timeLabel = 'Elapsed time (h)';
validation_options.saveFigures = true;
validation_options.showFigures = true;
validation_options.outputDir = fullfile(case_dir,['validation_' casename '_timeseries']);
% Matching uses values in obs(:,t) and bc(:,t), NOT presumed gauge row order.
% If their time units differ, convert them to the same units before this call.
validation_options.timeTolerance = 1e-6;

[validation_summary, validation, series] = ...
    plot_EWAVE_timeseries(output, input, casename, validation_options);


%% HIGH-TIDE CROSS-SHORE TOTAL IG-WAVE HEIGHT (G0 INCLUDED)
% This plot is independent of the G1-G4 time-series validation above.
% The highest observed/modelled still-water level in bc_xsw.inp is selected.
% Set highTideOptions.recordIndex to a positive integer for a specific record.
% G0 is identified by the offshore boundary coordinate (be(1,1)), not by
% assuming that observation rows have a particular gauge ordering.
% If obs_xsw.inp does not contain G0, its measured total IG height is taken
% from the prescribed total IG boundary condition (bc_xsw.inp).
highTideOptions = struct();
highTideOptions.recordIndex = [29];      % []: automatic highest-SWL record; or a specific BC row number
highTideOptions.minDepth = 0.5;      % Same minimum water depth as the time-series validation (m)
highTideOptions.timeTolerance = validation_options.timeTolerance;
highTideOptions.maxGaugeNumber = 5;  % Show G0-G5 when present; G5 may be too shallow
highTideOptions.showFigures = true;
highTideOptions.saveFigures = true;
highTideOptions.outputDir = validation_options.outputDir;

[highTideInfo, highTideProfile] = ...
    make_Xisha_high_tide_IG_profile(output, input, highTideOptions);


%% G0 TIME-SERIES COMPARISON (BOUNDARY STATION)
% Compare the measured total IG-wave height used at the offshore boundary
% with the modelled total IG-wave height at G0 over the full record.
% This is a boundary-consistency check rather than an independent validation.
g0Options = struct();
g0Options.timeTolerance = validation_options.timeTolerance;
g0Options.showFigures = true;
g0Options.saveFigures = true;
g0Options.outputDir = validation_options.outputDir;
g0Options.showTideShade = true;
g0Options.tideShadeAlpha = 0.17;
g0Options.timeLabel = validation_options.timeLabel;
g0Options.minutesPerStep = validation_options.minutesPerStep;

[g0_summary, g0_series] = make_Xisha_G0_timeseries(output, input, casename, g0Options);


function [info, profile] = make_Xisha_high_tide_IG_profile(output, input, opts)
% Compare observed and modelled total IG wave heights at one high-tide time.
% Plots the continuous model transect plus G0-G5 observations, and topography.
% Requires output(it).Hig_tt, output(it).ncal, input.be/bc/obs/index.
idx = input.index;
if ~exist(opts.outputDir,'dir'), mkdir(opts.outputDir); end
swl = input.bc(:,idx.swl);
times = input.bc(:,idx.t);
if numel(output)~=numel(times)
    error('The number of model outputs and boundary records does not match.');
end
if isempty(opts.recordIndex)
    eligible = find(isfinite(swl) & isfinite(times));
    if isempty(eligible), error('No finite SWL and time values in boundary input.'); end
    [~,p] = max(swl(eligible));
    it = eligible(p);
else
    it = opts.recordIndex;
    if ~isscalar(it) || ~isfinite(it) || it~=fix(it) || it<1 || it>numel(output)
        error('highTideOptions.recordIndex must be a valid BC row number.');
    end
end
if ~isfield(output(it),'Hig_tt') || ~isfield(output(it),'ncal')
    error('output(%d) has no Hig_tt or ncal field.',it);
end

% Model x-coordinates follow the same grid as the beach profile.
ncal = min([size(input.be,1),numel(output(it).Hig_tt),double(output(it).ncal)]);
if ncal<2, error('The selected record has fewer than two wet model nodes.'); end
modelX = input.be(1:ncal,1);
modelH = output(it).Hig_tt(1:ncal);
modelH = modelH(:);
validModel = isfinite(modelX) & isfinite(modelH) & modelH>=0;
modelX = modelX(validModel); modelH = modelH(validModel);
[modelX,ord] = sort(modelX); modelH = modelH(ord);
[modelX,uniq] = unique(modelX,'stable'); modelH = modelH(uniq);
if numel(modelX)<2, error('Not enough valid model IG wave-height values.'); end
beX = input.be(:,1);
beZ = input.be(:,2);
[gx,ix] = unique(beX,'stable');
gz = beZ(ix);
[gx,order] = sort(gx); gz = gz(order);

% Define G0 from the prescribed offshore grid coordinate. Other gauges are
% sorted landward. This avoids the original time-series plot's gauge-index
% assumption when the observations include G0.
g0x = beX(1);
gaugeCoordsObs = unique(input.obs(isfinite(input.obs(:,idx.x)),idx.x));
gaugeCoordsObs = sort(gaugeCoordsObs(:));
dx = min(abs(diff(gx)));
if isempty(dx) || ~isfinite(dx), dx = 1; end
coordTolerance = max(1e-6,min(0.5*dx,0.5));
% Include observation positions inside the model bathymetry domain only.
gaugeCoordsObs = gaugeCoordsObs(gaugeCoordsObs>=min(gx)-coordTolerance & ...
                                  gaugeCoordsObs<=max(gx)+coordTolerance);
otherGauges = gaugeCoordsObs(abs(gaugeCoordsObs-g0x)>coordTolerance & ...
                             gaugeCoordsObs>g0x);
otherGauges = otherGauges(1:min(numel(otherGauges),opts.maxGaugeNumber));
gaugeX = [g0x; otherGauges(:)];
if numel(gaugeX)<2
    warning('Only G0 was located. Check obs(:,x), be(:,x) and the x-coordinate datum.');
end
nGauge = numel(gaugeX);
gaugeNumber = (0:nGauge-1)';
gaugeLabels = cell(nGauge,1);
for ig = 1:nGauge, gaugeLabels{ig} = sprintf('G%d',gaugeNumber(ig)); end

% Match OBSERVATIONS to the chosen time using their actual time values;
% do not assume that gauge data are stored in consecutive rows.
matchTime = isfinite(input.obs(:,idx.t)) & ...
    abs(input.obs(:,idx.t)-times(it))<=opts.timeTolerance;
obsX = input.obs(:,idx.x);
obsH = input.obs(:,idx.Hig_tt);
obsAt = nan(nGauge,1);
source = repmat({'Not available'},nGauge,1);
for ig=1:nGauge
    hits = find(matchTime & isfinite(obsX) & ...
        abs(obsX-gaugeX(ig))<=coordTolerance & isfinite(obsH) & obsH>=0);
    if ~isempty(hits)
        obsAt(ig) = mean(obsH(hits));
        source{ig} = 'Observed gauge';
        if numel(hits)>1
            warning('Time %.6g, %s: %d matching rows; using their mean.', ...
                times(it),gaugeLabels{ig},numel(hits));
        end
    end
end
% G0 is the prescribed boundary station. Fall back to the measured BC
% total IG height if no G0 record is contained in the observation file.
if ~isfinite(obsAt(1)) && isfield(idx,'Hig_tt')
    boundaryHeight = input.bc(it,idx.Hig_tt);
    if isfinite(boundaryHeight) && boundaryHeight>0
        obsAt(1) = boundaryHeight;
        source{1} = 'Measured total IG boundary (BC file)';
    end
end

modelAt = interp1(modelX,modelH,gaugeX,'linear',NaN);
gaugeZ = interp1(gx,gz,gaugeX,'linear',NaN);
gaugeDepth = swl(it)-gaugeZ;
keep = isfinite(obsAt) & isfinite(modelAt) & ...
    isfinite(gaugeDepth) & gaugeDepth>=opts.minDepth;
% G0 is always identified and displayed when data exist; it is NOT an
% independent validation point because its observed height sets the BC.
independent = keep & gaugeNumber>=1 & gaugeNumber<=4;
nComp = sum(independent);
if nComp>0
    rmse = sqrt(mean((modelAt(independent)-obsAt(independent)).^2));
    if sum(obsAt(independent))>0
        rb = 100*sum(modelAt(independent)-obsAt(independent))/sum(obsAt(independent));
    else
        rb = NaN;
    end
else
    rmse=NaN; rb=NaN;
end

profile = table(gaugeNumber,gaugeLabels,gaugeX,gaugeZ,gaugeDepth, ...
    obsAt,modelAt,modelAt-obsAt,keep,independent,source, ...
    'VariableNames',{'GaugeNumber','Gauge','x_m','bedElevation_m', ...
    'waterDepth_m','IG_observed_m','IG_modelled_m','ModelMinusObs_m', ...
    'IncludedInPlot','IncludedInIndependentMetrics','ObservationSource'});
info = struct('recordIndex',it,'recordTime',times(it), ...
    'stillWaterLevel_m',swl(it),'RMSE_G1toG4_m',rmse, ...
    'RB_G1toG4_percent',rb,'N_G1toG4',nComp, ...
    'G0Source',source{1});

fprintf('\nXISHA BAY HIGH-TIDE IG CROSS-SHORE PROFILE\n');
fprintf('BC record %d | time ID = %.8g | SWL = %.4f m\n',it,times(it),swl(it));
fprintf('Independent paired points G1-G4: N=%d, RMSE=%.4f m, RB=%.2f%%\n', ...
    nComp,rmse,rb);
disp(profile);
writetable(profile,fullfile(opts.outputDir,'IG_high_tide_profile_G0_G5.csv'));
save(fullfile(opts.outputDir,'IG_high_tide_profile_G0_G5.mat'), ...
    'info','profile','opts');

if opts.showFigures, visible='on'; else, visible='off'; end
fig = figure('Color','w','Visible',visible,'Position',[130 100 1000 710]);

ax1 = subplot(2,1,1);
hModel = plot(ax1,modelX,modelH,'-','LineWidth',2, ...
    'Color',[0.13 0.42 0.70]);
hold(ax1,'on');
plotMask = keep & gaugeNumber>0;
hMeasured = plot(ax1,gaugeX(plotMask),obsAt(plotMask),'ko', ...
    'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.2);
% Special symbol makes the prescribed boundary G0 obvious.
hG0 = plot(ax1,gaugeX(1),obsAt(1),'s','MarkerSize',8, ...
    'MarkerEdgeColor',[0.68 0.16 0.16], ...
    'MarkerFaceColor',[1 0.87 0.84],'LineWidth',1.2);
if ~isfinite(obsAt(1)), set(hG0,'Visible','off'); end
for ig=1:nGauge
    if keep(ig)
        text(ax1,gaugeX(ig),obsAt(ig),['  ' gaugeLabels{ig}], ...
            'FontSize',9,'VerticalAlignment','bottom');
    end
end
xlabel(ax1,'Cross-shore distance (m)');
ylabel(ax1,'Total IG wave height, H_{rms} (m)','Interpreter','tex');
title(ax1,sprintf('Xisha Bay: high tide, time ID %.8g, SWL = %.3f m', ...
      times(it),swl(it)),'FontWeight','normal');
legend(ax1,[hModel,hMeasured,hG0], ...
    {'Modelled total IG','Measured G1-G5','Measured G0 (boundary)'}, ...
    'Location','best','Box','off');
if nComp>0
    text(ax1,0.98,0.95,sprintf('G1-G4: RMSE = %.3f m, RB = %.1f%% (N = %d)', ...
        rmse,rb,nComp),'Units','normalized','HorizontalAlignment','right', ...
        'VerticalAlignment','top','FontSize',9);
end
grid(ax1,'on'); box(ax1,'on');

ax2 = subplot(2,1,2);
plot(ax2,gx,gz,'k-','LineWidth',1.5); hold(ax2,'on');
plot(ax2,[min(gx) max(gx)],swl(it)*[1 1],'--', ...
    'Color',[0.3 0.55 0.75],'LineWidth',1.2);
plot(ax2,gaugeX,gaugeZ,'v','MarkerFaceColor',[0.90 0.46 0.37], ...
    'MarkerEdgeColor','k','MarkerSize',7,'LineStyle','none');
for ig=1:nGauge
    if isfinite(gaugeZ(ig))
        text(ax2,gaugeX(ig),gaugeZ(ig),['  ' gaugeLabels{ig}], ...
            'FontSize',9,'VerticalAlignment','bottom');
    end
end
xlabel(ax2,'Cross-shore distance (m)');
ylabel(ax2,'Elevation (m)');
title(ax2,'Beach profile, gauge positions and selected high-tide level', ...
    'FontWeight','normal');
legend(ax2,{'Beach profile','Still-water level','Gauges G0-G5'}, ...
    'Location','best','Box','off');
grid(ax2,'on'); box(ax2,'on');
linkaxes([ax1 ax2],'x');
xlim(ax1,[min(gx) max(gx)]);

if opts.saveFigures
    saveas(fig,fullfile(opts.outputDir,'IG_high_tide_crossshore_G0_G5.png'));
    savefig(fig,fullfile(opts.outputDir,'IG_high_tide_crossshore_G0_G5.fig'));
    try
        print(fig,fullfile(opts.outputDir,'IG_high_tide_crossshore_G0_G5.pdf'), ...
            '-dpdf','-painters');
    catch ME
        warning('Could not export high-tide profile as PDF: %s',ME.message);
    end
end
if ~opts.showFigures, close(fig); end
end


function [summary, series] = make_Xisha_G0_timeseries(output, input, casename, opts)
% Plot the total IG-wave-height time series at G0 (offshore boundary).
% Uses observed G0 values from obs file when available; otherwise falls back
% to the measured total IG boundary prescribed in bc_xsw.inp.
idx = input.index;
if ~exist(opts.outputDir,'dir'), mkdir(opts.outputDir); end
ntime = numel(output);
bcTime = input.bc(:,idx.t);
if ntime~=numel(bcTime)
    error('The number of model outputs and BC records does not match.');
end

% Time vector in hours for plotting; use actual t values when they are
% monotonic and finite, otherwise use a uniform step fallback.
if all(isfinite(bcTime)) && numel(unique(bcTime))>1
    tPlot = bcTime - bcTime(1);
    xLabel = 'Time relative to first record';
    timeUnit = 'BC time units';
else
    tPlot = (0:ntime-1)' * opts.minutesPerStep / 60;
    xLabel = opts.timeLabel;
    timeUnit = 'hours';
end

% Identify G0 from the offshore boundary coordinate.
g0x = input.be(1,1);
beX = input.be(:,1);
dx = min(abs(diff(unique(beX))));
if isempty(dx) || ~isfinite(dx), dx = 1; end
coordTolerance = max(1e-6,min(0.5*dx,0.5));

obsTime = input.obs(:,idx.t);
obsX = input.obs(:,idx.x);
obsH = input.obs(:,idx.Hig_tt);
obsG0 = nan(ntime,1);
obsSource = repmat({''},ntime,1);
for it=1:ntime
    hits = find(isfinite(obsTime) & abs(obsTime-bcTime(it))<=opts.timeTolerance & ...
                isfinite(obsX) & abs(obsX-g0x)<=coordTolerance & ...
                isfinite(obsH) & obsH>=0);
    if ~isempty(hits)
        obsG0(it) = mean(obsH(hits));
        obsSource{it} = 'Observed G0 in obs file';
    else
        boundaryHeight = input.bc(it,idx.Hig_tt);
        if isfinite(boundaryHeight) && boundaryHeight>0
            obsG0(it) = boundaryHeight;
            obsSource{it} = 'Measured total IG boundary (BC file)';
        else
            obsSource{it} = 'Not available';
        end
    end
end

modelG0 = nan(ntime,1);
swl = nan(ntime,1);
for it=1:ntime
    swl(it) = input.bc(it,idx.swl);
    if isfield(output(it),'Hig_tt') && ~isempty(output(it).Hig_tt)
        modelG0(it) = output(it).Hig_tt(1);
    end
end

valid = isfinite(obsG0) & isfinite(modelG0) & obsG0>=0 & modelG0>=0;
N = sum(valid);
if N>0
    rmse = sqrt(mean((modelG0(valid)-obsG0(valid)).^2));
    if sum(obsG0(valid))>0
        rb = 100*sum(modelG0(valid)-obsG0(valid))/sum(obsG0(valid));
    else
        rb = NaN;
    end
else
    rmse = NaN; rb = NaN;
end
summary = table(N,rmse,rb,'VariableNames',{'N','RMSE_m','RelativeBias_percent'});
series = table((1:ntime)',bcTime,tPlot,swl,obsG0,modelG0,obsSource,valid, ...
    'VariableNames',{'Record','TimeID','PlotTime','SWL_m','IG_observed_m', ...
    'IG_modelled_m','ObservationSource','ValidPair'});

fprintf('G0 TOTAL IG TIME SERIES');
fprintf('Valid paired records: N=%d, RMSE=%.4f m, RB=%.2f%%',N,rmse,rb);
disp(summary);
writetable(series,fullfile(opts.outputDir,'IG_G0_timeseries.csv'));
save(fullfile(opts.outputDir,'IG_G0_timeseries.mat'),'summary','series','opts');

if opts.showFigures, visible='on'; else, visible='off'; end
fig = figure('Color','w','Visible',visible,'Position',[150 120 980 420]);
ax = axes('Parent',fig); hold(ax,'on');
% Optional grey shading for tidal level, scaled to current y-range later.
if opts.showTideShade && any(isfinite(swl))
    yyaxis(ax,'right');
    area(ax,tPlot,swl,'LineStyle','none','FaceColor',[0.6 0.6 0.6], ...
        'FaceAlpha',opts.tideShadeAlpha);
    ylabel(ax,'Still-water level (m)');
    set(ax,'YColor',[0.45 0.45 0.45]);
    yyaxis(ax,'left');
end
hObs = plot(ax,tPlot,obsG0,'ko','MarkerSize',5,'LineWidth',1.0, ...
    'MarkerFaceColor','w');
hMod = plot(ax,tPlot,modelG0,'-','LineWidth',1.8,'Color',[0.13 0.42 0.70]);
if strcmp(timeUnit,'hours')
    xlabel(ax,opts.timeLabel);
else
    xlabel(ax,xLabel);
end
ylabel(ax,'Total IG wave height at G0, H_{rms} (m)','Interpreter','tex');
title(ax,'Xisha Bay: G0 total IG-wave-height time series','FontWeight','normal');
legend(ax,[hObs,hMod],{'Measured G0 / prescribed boundary','Modelled G0'}, ...
    'Location','best','Box','off');
text(ax,0.98,0.95,sprintf('RMSE = %.3f m, RB = %.1f%% (N = %d)',rmse,rb,N), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',9);
text(ax,0.02,0.95,'G0 is the boundary station; this figure checks boundary consistency, not independent validation.', ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','top','FontSize',8);
grid(ax,'on'); box(ax,'on');
if opts.saveFigures
    saveas(fig,fullfile(opts.outputDir,'IG_G0_timeseries.png'));
    savefig(fig,fullfile(opts.outputDir,'IG_G0_timeseries.fig'));
    try
        print(fig,fullfile(opts.outputDir,'IG_G0_timeseries.pdf'),'-dpdf','-painters');
    catch ME
        warning('Could not export G0 time series as PDF: %s',ME.message);
    end
end
if ~opts.showFigures, close(fig); end
end
