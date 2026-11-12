%% WME coefficient sensitivity: pooled RMSE and relative bias (RB)
% Based on Sensitive.m; only the statistics, plotting and input checking
% are revised. The original EWAVE solver and physical setup are unchanged.
%
% For each alpha_sw / alpha_r pair, all valid observation-model pairs
% across the selected time records are pooled BEFORE calculating metrics:
%   RMSE = sqrt( sum((M-O).^2) / N )                        [m]
%   RB   = 100 * sum(M-O) / sum(O)                         [%]
% No median of individual-record RMSEs; no R^2 or RMSPE.
% The two WME partition coefficients are independent (not required to sum to 1).
% MATLAB R2016b+ (script with local functions).

clearvars; close all; clc;

%% 0. Paths and case (edit for your computer)
caseDir = 'C:\Users\DELL\Desktop\EWAVE_v1.0_1003\cases\beach\SW_IG_RF_Globex JPO2015';
srcDir  = 'C:\Users\DELL\Desktop\EWAVE_v1.0_1003\src';
addpath(caseDir);
addpath(genpath(srcDir));

casename       = 'A1';
input.be_file  = fullfile(caseDir, ['be_'  casename '.inp']);
input.bc_file  = fullfile(caseDir, ['bc_'  casename '.inp']);
input.obs_file = fullfile(caseDir, ['obs_' casename '.inp']);

input.index.t = 1; input.index.x = 2; input.index.swl = 3;
input.index.mwl = 4; input.index.spectype = 5;
input.index.Hsw_tt = 6; input.index.Hsw_in = 7;
input.index.Hsw_out = 8; input.index.Tsw = 9;
input.index.angle_sw = 10;
input.index.Hig_tt = 11; input.index.Hig_in = 12;
input.index.Hig_out = 13; input.index.Tig = 14;
input.index.angle_ig = 15; input.index.Lr_sw = 16;
input.index.Lr_ig = 17; input.index.V_wind = 18;

%% 1. Model configuration: retained from the uploaded Sensitive.m
setup.wave_type = 'random';
setup.B_formula = '0.125';
setup.c_formula = 'nonlinear';
setup.Db_formula = 'B1998';
setup.gamma_formula = 'zl2020';
setup.reflected_wave = 'yes';
setup.Kr_formula = {'M1983','ZV2008','corrected','1'};
setup.ig_wave = 'yes';
setup.mwl_calculation = 'yes';
setup.surface_roller = 'yes';
setup.beta_formula = 'Z2017';
setup.roller_IG_transfer = '0.7';  % alpha_r (reference)
setup.SW_IG_transfer = '0.3';      % alpha_sw (reference)
setup.rhor = '1000';
setup.s0_method = 'offshore';
setup.d50 = '0.0';
setup.Cf = '0.0';
setup.Tig = '10';
setup.dx = '0.1';
setup.kh_cut = '0.1';
setup.relax = '0';
setup.converg = '0.001';
setup.N_iteration = '20';
setup.mwl_correction = 'no';
setup.current_advection = 'no';
setup.igbc_LHS62 = {'yes','1'};
% These switches remain OFF, as in your current WME-only calculations.
% If the current solver has removed these fields, unused fields do no harm.
setup.RS_work_swin = 'no';
setup.RS_work_swout = 'no';
setup.RS_work_roller = 'no';

%% 2. Sensitivity parameters and plotting options
alpha_sw_vec = 0.1:0.1:0.7;
alpha_r_vec  = 0.3:0.1:0.9;

% Use 'mm' for the laboratory-scale colour bars or 'm' to follow the
% unit used in the paper's statistical definition. CSV always stores m.
plotRMSEUnit = 'mm';

% If a future field-case sensitivity test needs a minimum local depth,
% set the threshold here. For GLOBEX/Boers leave it at 0 (no depth cutoff).
minWaterDepth = 0;  % m

% If one solver run fails, record its parameter pair as NaN and continue.
% Set false to stop immediately and inspect the solver error instead.
continueAfterSolverError = false;

saveDir = fullfile(pwd, ['sensitivity_' casename '_RMSE_RB']);
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% 3. Read inputs and retain the original wave-height conversion
input.be  = read_numeric_input(input.be_file);
input.bc  = read_numeric_input(input.bc_file);
input.obs = read_numeric_input(input.obs_file);

if size(input.be,2) < 2 || size(input.bc,2) < 18 || size(input.obs,2) < 18
    error('Unexpected input column count in be/bc/obs file.');
end
if any(diff(input.be(:,1)) <= 0)
    error('Beach-profile x coordinates must be strictly increasing.');
end

x = (input.be(1,1):str2double(setup.dx):input.be(end,1)).';
zb = interp1(input.be(:,1),input.be(:,2),x,'linear');
input.be = [x,zb];

% spectype = 1: RMS wave height; spectype = 2: significant height.
% Preserve the conversion used in the original script.
waveCols = [input.index.Hsw_tt input.index.Hsw_in input.index.Hsw_out, ...
            input.index.Hig_tt input.index.Hig_in input.index.Hig_out];
for source = {'bc','obs'}
    field = source{1};
    raw = input.(field);
    typ = raw(:,input.index.spectype);
    if any(~isfinite(typ) | typ <= 0)
        error('%s: spectype must be finite and > 0 for all rows.',field);
    end
    raw(:,waveCols) = bsxfun(@rdivide,raw(:,waveCols),sqrt(typ));
    raw(:,input.index.spectype) = 1;
    input.(field) = raw;
end

time = input.bc(:,input.index.t);
nTime = numel(time);
if nTime < 1 || mod(size(input.obs,1),nTime) ~= 0
    error('Observation rows must be grouped into equal-size time blocks.');
end
ngauge = size(input.obs,1)/nTime;
if ngauge < 1, error('No observation rows found.'); end

% Original code assumes observations ordered in nTime blocks of ngauge rows.
% Verify that assumption, rather than silently comparing different times.
for it = 1:nTime
    idx = (it-1)*ngauge + (1:ngauge);
    tObs = input.obs(idx,input.index.t);
    if any(~isfinite(tObs)) || any(abs(tObs-time(it)) > ...
            1e-8*max(1,abs(time(it))))
        warning(['Observation time labels differ from the boundary record ', ...
            'at block %d. Original fixed-block ordering is retained; ', ...
            'please confirm the rows correspond to the same time.'],it);
    end
    if isnan(input.bc(it,input.index.Tig))
        input.bc(it,input.index.Tig) = str2double(setup.Tig)* ...
            input.bc(it,input.index.Tsw);
    end
    if isnan(input.bc(it,input.index.Hsw_in))
        input.bc(it,input.index.Hsw_in) = input.bc(it,input.index.Hsw_tt);
    end
    if isnan(input.bc(it,input.index.Hig_in))
        input.bc(it,input.index.Hig_in) = input.bc(it,input.index.Hig_tt);
    end
    input.bc(it,isnan(input.bc(it,:))) = 0;
end

%% 4. Sweep: pool ALL valid pairs for each component and coefficient pair
nSW = numel(alpha_sw_vec);
nR  = numel(alpha_r_vec);
compVars = {'Hig_tt','Hig_in','Hig_out'};
nComp = numel(compVars);

% Matrix size: [nR x nSW x 3]; dimensions correspond to alpha_r,
% alpha_sw, and (total, incoming, outgoing IG wave height).
RMSE = nan(nR,nSW,nComp);  % m
RB   = nan(nR,nSW,nComp);  % percent, signed
N    = zeros(nR,nSW,nComp); % number of VALID paired values

setup0 = setup;
fprintf('\nWME parameter sensitivity: %s\n',casename);
fprintf('Coefficient pairs: %d; time records: %d; gauges/block: %d\n', ...
    nR*nSW,nTime,ngauge);

for ir = 1:nR
    for isw = 1:nSW
        setup = setup0;
        setup.SW_IG_transfer     = num2str(alpha_sw_vec(isw),'%.12g');
        setup.roller_IG_transfer = num2str(alpha_r_vec(ir),'%.12g');

        sumErr   = zeros(1,nComp);
        sumErr2  = zeros(1,nComp);
        sumObs   = zeros(1,nComp);
        nValid   = zeros(1,nComp);
        failed = false;

        for it = 1:nTime
            try
                out = EWAVE_functions_2020d( ...
                    setup,input.be,input.bc(it,:),input.index);
            catch ME
                if ~continueAfterSolverError, rethrow(ME); end
                warning('Solver failed at alpha_sw=%.3f, alpha_r=%.3f: %s', ...
                    alpha_sw_vec(isw),alpha_r_vec(ir),ME.message);
                failed = true;
                break
            end

            ncal = out.ncal;
            if ~isscalar(ncal) || ~isfinite(ncal) || ncal < 2 || ...
                    ncal > size(input.be,1)
                error('Invalid out.ncal for alpha_sw=%.3f, alpha_r=%.3f.', ...
                      alpha_sw_vec(isw),alpha_r_vec(ir));
            end
            xModel = input.be(1:ncal,1);
            idx = (it-1)*ngauge + (1:ngauge);
            xObs = input.obs(idx,input.index.x);

            % Optional water-depth screening based on still-water level
            % and bathymetry; disabled by default for laboratory cases.
            if minWaterDepth > 0
                zObs = interp1(input.be(:,1),input.be(:,2),xObs,'linear',NaN);
                depthObs = input.bc(it,input.index.swl) - zObs;
                wet = isfinite(depthObs) & depthObs >= minWaterDepth;
            else
                wet = true(size(xObs));
            end

            for k = 1:nComp
                name = compVars{k};
                O = input.obs(idx,input.index.(name));
                fieldValues = out.(name);
                M = interp1(xModel,fieldValues(1:ncal),xObs,'linear',NaN);
                O = O(:); M = M(:);

                % Exclude NaN, Inf, negative wave heights, outside-domain
                % interpolation and (optional) too-shallow observations.
                valid = wet(:) & isfinite(O) & isfinite(M) & O>=0 & M>=0;
                if ~any(valid), continue; end
                err = M(valid)-O(valid);
                sumErr(k)  = sumErr(k)+sum(err);
                sumErr2(k) = sumErr2(k)+sum(err.^2);
                sumObs(k)  = sumObs(k)+sum(O(valid));
                nValid(k)  = nValid(k)+numel(err);
            end
        end

        if failed, continue; end
        for k = 1:nComp
            N(ir,isw,k) = nValid(k);
            if nValid(k) > 0
                RMSE(ir,isw,k) = sqrt(sumErr2(k)/nValid(k));
            end
            if sumObs(k) > 0
                RB(ir,isw,k) = 100*sumErr(k)/sumObs(k);
            end
        end
        fprintf('alpha_sw=%.2f, alpha_r=%.2f | RMSE(mm): %.3f %.3f %.3f | RB(%%): %+.2f %+.2f %+.2f\n', ...
                alpha_sw_vec(isw),alpha_r_vec(ir), ...
                1000*RMSE(ir,isw,1),1000*RMSE(ir,isw,2),1000*RMSE(ir,isw,3), ...
                RB(ir,isw,1),RB(ir,isw,2),RB(ir,isw,3));
    end
end
setup = setup0;

%% 5. Export the full 2D parameter scan (each row is one coefficient pair)
[X,Y] = meshgrid(alpha_sw_vec,alpha_r_vec);
results = table(X(:),Y(:), ...
    reshape(RMSE(:,:,1),[],1),reshape(RB(:,:,1),[],1),reshape(N(:,:,1),[],1), ...
    reshape(RMSE(:,:,2),[],1),reshape(RB(:,:,2),[],1),reshape(N(:,:,2),[],1), ...
    reshape(RMSE(:,:,3),[],1),reshape(RB(:,:,3),[],1),reshape(N(:,:,3),[],1), ...
    'VariableNames',{'alpha_sw','alpha_r', ...
    'RMSE_total_m','RB_total_percent','N_total', ...
    'RMSE_incoming_m','RB_incoming_percent','N_incoming', ...
    'RMSE_outgoing_m','RB_outgoing_percent','N_outgoing'});

csvFile = fullfile(saveDir,['sensitivity_' casename '_RMSE_RB.csv']);
writetable(results,csvFile);
matFile = fullfile(saveDir,['sensitivity_' casename '_RMSE_RB.mat']);
save(matFile,'RMSE','RB','N','X','Y','alpha_sw_vec','alpha_r_vec', ...
    'compVars','results','casename','setup0','minWaterDepth');

%% 6. Plot: 2 rows x 3 columns (RMSE, then signed relative bias)
if strcmpi(plotRMSEUnit,'mm')
    rmseScale = 1000;
    rmseAxisLabel = 'RMSE (mm)';
elseif strcmpi(plotRMSEUnit,'m')
    rmseScale = 1;
    rmseAxisLabel = 'RMSE (m)';
else
    error('plotRMSEUnit must be ''mm'' or ''m''.');
end

fig = figure('Color','w','Units','pixels','Position',[80 80 1380 750]);
panelNames = {'Total IG wave height', ...
              'Incoming IG wave height', ...
              'Outgoing IG wave height'};
alphaRefSw = str2double(setup0.SW_IG_transfer);
alphaRefR  = str2double(setup0.roller_IG_transfer);

for k = 1:nComp
    % Upper row: pooled RMSE for each IG component
    ax = subplot(2,3,k);
    Z = RMSE(:,:,k)*rmseScale;
    contourf(ax,X,Y,Z,12,'LineStyle','none');
    colormap(ax,parula(256));
    cb = colorbar(ax); cb.Label.String = rmseAxisLabel;
    style_sensitivity_axes(ax,alpha_sw_vec,alpha_r_vec, ...
        alphaRefSw,alphaRefR);
    title(ax,sprintf('(%s) %s',char('a'+k-1),panelNames{k}));

    % Lower row: signed relative bias, 0 centred in diverging colour scale
    ax = subplot(2,3,k+3);
    Z = RB(:,:,k);
    contourf(ax,X,Y,Z,12,'LineStyle','none');
    colormap(ax,blue_white_red(256));
    maxAbs = max(abs(Z(isfinite(Z))));
    if isempty(maxAbs) || maxAbs < 1, maxAbs = 1; end
    caxis(ax,[-maxAbs maxAbs]);
    cb = colorbar(ax); cb.Label.String = 'RB (%)';
    style_sensitivity_axes(ax,alpha_sw_vec,alpha_r_vec, ...
        alphaRefSw,alphaRefR);
    title(ax,sprintf('(%s) %s',char('d'+k-1),panelNames{k}));
end

set(findall(fig,'Type','axes'),'FontName','Times New Roman','FontSize',11);
% Plot uses a star for the original reference settings (0.3,0.7).
annotation(fig,'textbox',[0.23 0.006 0.60 0.025], ...
    'String',sprintf('Black star: reference parameters (\\alpha_{sw} = %.2f, \\alpha_r = %.2f)', ...
    alphaRefSw,alphaRefR), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'Interpreter','tex','FontName','Times New Roman','FontSize',10);

figFile = fullfile(saveDir,['sensitivity_' casename '_RMSE_RB.fig']);
pngFile = fullfile(saveDir,['sensitivity_' casename '_RMSE_RB.png']);
savefig(fig,figFile);
print(fig,pngFile,'-dpng','-r300');

%% 7. Report results for the reference pair (if present in the scan)
[~,isRef] = min(abs(alpha_sw_vec-alphaRefSw));
[~,irRef] = min(abs(alpha_r_vec-alphaRefR));
if abs(alpha_sw_vec(isRef)-alphaRefSw)<1e-9 && ...
   abs(alpha_r_vec(irRef)-alphaRefR)<1e-9
    fprintf('\nREFERENCE: alpha_sw=%.2f, alpha_r=%.2f\n',alphaRefSw,alphaRefR);
    for k = 1:nComp
        fprintf('%s: N=%d, RMSE=%.6f m, RB=%+.2f%%\n',compVars{k}, ...
            N(irRef,isRef,k),RMSE(irRef,isRef,k),RB(irRef,isRef,k));
    end
end
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', ...
    csvFile,matFile,figFile,pngFile);

%% Local helper functions
function A = read_numeric_input(filename)
    if ~exist(filename,'file')
        error('Cannot find file: %s',filename);
    end
    data = importdata(filename);
    if isstruct(data)
        if ~isfield(data,'data') || isempty(data.data)
            error('No numeric data found in: %s',filename);
        end
        A = data.data;
    else
        A = data;
    end
    if ~isnumeric(A) || isempty(A)
        error('Input is not a numeric matrix: %s',filename);
    end
end

function style_sensitivity_axes(ax,aSW,aR,alphaSwRef,alphaRRef)
    hold(ax,'on');
    plot(ax,alphaSwRef,alphaRRef,'kp', ...
        'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1.4);
    xlabel(ax,'\alpha_{sw}','Interpreter','tex');
    ylabel(ax,'\alpha_r','Interpreter','tex');
    set(ax,'XLim',[min(aSW) max(aSW)],'YLim',[min(aR) max(aR)], ...
        'XTick',aSW,'YTick',aR,'TickDir','out','Box','on', ...
        'Layer','top');
    axis(ax,'square');
end

function cmap = blue_white_red(n)
    n1 = ceil(n/2); n2 = n-n1;
    cmap = [linspace(0.16,1,n1)',linspace(0.37,1,n1)',ones(n1,1); ...
            ones(n2,1),linspace(1,0.23,n2)',linspace(1,0.24,n2)'];
end
