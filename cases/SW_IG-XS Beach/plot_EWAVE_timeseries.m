function [validation_summary, validation, series] = plot_EWAVE_timeseries(output, input, casename, options)
% Validation of TOTAL IG wave height only: all gauges in one figure,
% with grey transparent tidal-level shading (right axis) in each selected gauge
% panel and a shared beach profile with gauges at the bottom.
% Metrics pool all valid gauge-time observation/model pairs:
%   RMSE = sqrt(mean((model - observed).^2)) [m]
%   RB = 100*sum(model-observed)/sum(observed) [%]
% Each gauge also receives its own RMSE and RB.
% Model/observation wave heights must use the same conventions (RMS here).
%
% This function assumes observation-time column and model BC time column
% use identical units and origins. No interpolation in time is performed.

if nargin < 4, options = struct(); end
options = optDefault(options, 'minDepth', 0.5);
options = optDefault(options, 'gaugeNumbers', 1:4); % Only plot/assess Gauge 1-4
options = optDefault(options, 'showTideShade', true);
options = optDefault(options, 'tideShadeColour', [0.50 0.50 0.50]);
options = optDefault(options, 'tideShadeAlpha', 0.17);
options = optDefault(options, 'timeTolerance', 1e-6);
options = optDefault(options, 'minutesPerStep', 12.8);
options = optDefault(options, 'timeLabel', 'Elapsed time (h)');
options = optDefault(options, 'saveFigures', true);
options = optDefault(options, 'showFigures', true);
options = optDefault(options, 'outputDir', fullfile(pwd,['validation_' casename '_timeseries']));
if ~exist(options.outputDir,'dir'), mkdir(options.outputDir); end

idx = input.index;
req = {'t','x','Hig_tt','swl'};
for k=1:numel(req)
    if ~isfield(idx,req{k}), error('Missing input.index.%s',req{k}); end
end
if ~isfield(input,'obs') || isempty(input.obs), error('Missing observation data'); end
if ~isfield(input,'bc') || isempty(input.bc), error('Missing boundary data'); end
if ~isfield(input,'be') || size(input.be,2)<2, error('Missing beach profile'); end

obsTime = input.obs(:,idx.t);
obsX = input.obs(:,idx.x);
obsH = input.obs(:,idx.Hig_tt);
modelTime = input.bc(:,idx.t);
if numel(modelTime) ~= numel(output)
    error('Number of model times (%d) differs from output records (%d).', numel(modelTime),numel(output));
end
validCoord = isfinite(obsX);
allGauges = unique(obsX(validCoord));
selected = options.gaugeNumbers(:);
if any(~isfinite(selected) | selected < 1 | selected ~= fix(selected)) || ...
        any(selected > numel(allGauges))
    error('Requested gauge numbers must be integers between 1 and %d.',numel(allGauges));
end
gauges = allGauges(selected);
nG = numel(gauges);
if nG == 0, error('No finite observation gauge coordinates'); end

% Map observation times to model time records, within configured tolerance.
% Exact/nearest record matching avoids merging rows by assumed gauge order.
modelIdx = nan(size(obsTime));
for j=1:numel(obsTime)
    if ~isfinite(obsTime(j)), continue; end
    [dist, ii] = min(abs(modelTime - obsTime(j)));
    if isfinite(dist) && dist <= options.timeTolerance, modelIdx(j) = ii; end
end
if ~any(isfinite(modelIdx))
    error(['No matching observation/model times. Check the units/origin of ', ...
           'input.obs(:,t), input.bc(:,t), and options.timeTolerance.']);
end
if sum(isfinite(modelIdx)) < sum(isfinite(obsTime))
    warning('%d observation rows have no matching model time and will be excluded.', ...
        sum(isfinite(obsTime) & ~isfinite(modelIdx)));
end

M = nan(size(obsH));
localDepth = nan(size(obsH));
% Cache interpolations at each time record.
for it=1:numel(output)
    rows = find(modelIdx == it & isfinite(obsX));
    if isempty(rows), continue; end
    if ~isfield(output(it),'Hig_tt') || ~isfield(output(it),'ncal')
        error('output(%d) requires fields Hig_tt and ncal.',it);
    end
    ncal = min(output(it).ncal, size(input.be,1));
    h = output(it).Hig_tt(:);
    ncal = min(ncal,numel(h));
    if ncal < 2, continue; end
    xx = input.be(1:ncal,1);
    zz = input.be(1:ncal,2);
    hh = h(1:ncal);
    % Bathymetry coordinate must match the model output coordinate.
    [xx, order] = sort(xx);
    zz = zz(order); hh = hh(order);
    [xx, uniqueIdx] = unique(xx, 'stable');
    zz = zz(uniqueIdx); hh = hh(uniqueIdx);
    if numel(xx)<2, continue; end
    M(rows) = interp1(xx,hh,obsX(rows),'linear',NaN);
    % Use model boundary SWL; include any prescribed tidal elevation.
    waterLevel = input.bc(it,idx.swl);
    % Depth is based on the full prescribed beach profile, NOT ncal.
    % ncal truncates the wave solution at the terminal wet model node.
    localDepth(rows) = waterLevel - interp1(input.be(:,1), ...
        input.be(:,2),obsX(rows),'linear',NaN);
end

% Preserve prior 0.5 m minimum-depth criterion; no arbitrary zero-depth filling.
valid = isfinite(obsH) & isfinite(M) & isfinite(localDepth) & ...
        obsH >= 0 & M >= 0 & localDepth >= options.minDepth & ...
        isfinite(modelIdx) & isfinite(obsX) & ismember(obsX,gauges);
if any(isfinite(M) & M < 0)
    warning('%d negative IG model wave heights excluded.',sum(isfinite(M) & M < 0));
end

N = sum(valid);
[overallRmse,overallRb] = metrics(obsH(valid),M(valid));
validation_summary = table({'Hig_tt'}, N, overallRmse, overallRb, ...
    'VariableNames',{'Variable','N','RMSE_m','RelativeBias_percent'});

nPer = zeros(nG,1); rmsePer = nan(nG,1); rbPer = nan(nG,1);
for ig=1:nG
    mask = valid & obsX == gauges(ig);
    nPer(ig) = sum(mask);
    [rmsePer(ig),rbPer(ig)] = metrics(obsH(mask),M(mask));
end
perGauge = table(gauges,nPer,rmsePer,rbPer, ...
    'VariableNames',{'x_m','N','RMSE_m','RelativeBias_percent'});
validation = struct('perGauge',perGauge,'valid',valid,'minDepth',options.minDepth, ...
    'timeTolerance',options.timeTolerance);
series = table(obsTime,obsX,obsH,M,localDepth,valid,modelIdx, ...
    'VariableNames',{'time','x_m','Hig_observed_m','Hig_modelled_m', ...
                     'localDepth_m','valid','modelTimeIndex'});

fprintf('\nIG HEIGHT VALIDATION: %s (total IG only)\n',casename);
disp(validation_summary);
disp(perGauge);
writetable(validation_summary,fullfile(options.outputDir,'IG_validation_summary.csv'));
writetable(perGauge,fullfile(options.outputDir,'IG_validation_per_gauge.csv'));
writetable(series,fullfile(options.outputDir,'IG_validation_pairs.csv'));
save(fullfile(options.outputDir,'IG_validation.mat'), ...
    'validation_summary','validation','series','options');

% Display time in hours relative to first BC record if BC times are
% minute indices with fixed steps (existing original-script convention).
% If BC/obs time columns already represent elapsed hours, set
% options.minutesPerStep = 0 to plot the original numeric time values.
t0 = modelTime(1);
if options.minutesPerStep > 0
    plotHours = (obsTime-t0)*options.minutesPerStep/60;
    tideHours = (modelTime-t0)*options.minutesPerStep/60;
else
    plotHours = obsTime-t0;
    tideHours = modelTime-t0;
end
% The prescribed still-water level at each BC time is the tidal signal.
tideLevels = input.bc(:,idx.swl);
tideHours = tideHours(:);
tideLevels = tideLevels(:);
tideOK = isfinite(tideHours) & isfinite(tideLevels);
if any(tideOK)
    tideLo = min(tideLevels(tideOK));
    tideHi = max(tideLevels(tideOK));
    tidePad = max(0.04*(tideHi-tideLo),0.05);
    tideLims = [tideLo-tidePad tideHi+tidePad];
else
    tideLims = [-1 1];
end

% One combined figure: ONLY selected IG gauge panels, followed
% by a common beach-profile panel with all gauge locations marked.
% Invalid/shallow records are kept as NaN on the plotting time axis so that
% a line is NOT drawn through periods with depth < options.minDepth.
beachX = input.be(:,1);
beachZ = input.be(:,2);
if options.showFigures, visibility='on'; else, visibility='off'; end
figHeight = max(680, min(2600, 180*(nG+1)));
fig = figure('Color','w','Visible',visibility, ...
    'Position',[100 50 1050 figHeight]);
axesList = gobjects(nG,1);
colorModel = [0.16 0.52 0.73];
colorObs = [0.72 0.88 0.93];

for ig=1:nG
    ax = subplot(nG+1,1,ig);
    axesList(ig) = ax;

    % Semi-transparent grey tidal-level area on the RIGHT vertical axis.
    % It is not a wave-height uncertainty band and is not depth-filtered.
    % With the secondary axis, the actual tidal range is shown in metres.
    if options.showTideShade && any(tideOK)
        yyaxis(ax,'right');
        hold(ax,'on');
        [tSorted,tid] = sort(tideHours(tideOK));
        swlSorted = tideLevels(tideOK);
        swlSorted = swlSorted(tid);
        % A light filled area from the minimum tidal level to SWL(t).
        % NaN boundary values are excluded, no artificial tide is filled.
        baseLevel = tideLo;
        hShade = patch(ax, [tSorted; flipud(tSorted)], ...
            [swlSorted; repmat(baseLevel,numel(swlSorted),1)], ...
            options.tideShadeColour, 'FaceAlpha',options.tideShadeAlpha, ...
            'EdgeColor','none', 'HandleVisibility','off'); %#ok<NASGU>
        hTide = plot(ax,tSorted,swlSorted,'-', ...
            'Color',[0.48 0.48 0.48],'LineWidth',0.8);
        ylim(ax,tideLims);
        ylabel(ax,'Tidal level (m)');
        ax.YAxis(2).Color = [0.42 0.42 0.42];
    end

    % Waves on the LEFT axis. Invalid/shallower than 0.5 m records
    % remain NaN to break model lines at dry or insufficiently wet times.
    yyaxis(ax,'left');
    hold(ax,'on');
    rows = find(obsX == gauges(ig) & isfinite(plotHours));
    [tPlot, order] = sort(plotHours(rows));
    rows = rows(order);
    oPlot = obsH(rows);
    mPlot = M(rows);
    oPlot(~valid(rows)) = NaN;
    mPlot(~valid(rows)) = NaN;
    hModel = plot(ax,tPlot,mPlot,'-','LineWidth',1.5,'Color',colorModel);
    hObs = plot(ax,tPlot,oPlot,'o','MarkerSize',3.5, ...
        'MarkerEdgeColor',[0.12 0.12 0.12], ...
        'MarkerFaceColor',colorObs);
    ylabel(ax,'H_{IG} (m)','Interpreter','tex');
    title(ax,sprintf('Gauge %d (x = %.2f m): RMSE = %.3f m, RB = %.2f%%, N = %d', ...
        selected(ig),gauges(ig),rmsePer(ig),rbPer(ig),nPer(ig)), ...
        'FontWeight','normal');
    if ig==1
        if options.showTideShade && any(tideOK)
            legend(ax,[hModel,hObs,hTide], ...
                {'Modelled IG wave height','Measured IG wave height','Tidal level (grey shading)'}, ...
                'Location','best','Box','off');
        else
            legend(ax,[hModel,hObs], ...
                {'Modelled IG wave height','Measured IG wave height'}, ...
                'Location','best','Box','off');
        end
    end
    if ig ~= nG
        set(ax,'XTickLabel',[]);
    else
        xlabel(ax,options.timeLabel);
    end
    grid(ax,'on'); box(ax,'on');
end
if nG > 1, linkaxes(axesList,'x'); end

axTopo = subplot(nG+1,1,nG+1);
plot(beachX,beachZ,'k-','LineWidth',1.6); hold on;
zGauges = interp1(beachX,beachZ,gauges,'linear',NaN);
plot(gauges,zGauges,'v','MarkerSize',7, ...
    'MarkerFaceColor',[0.89 0.35 0.39], 'MarkerEdgeColor','k', ...
    'LineStyle','none');
for ig=1:nG
    if isfinite(zGauges(ig))
        text(gauges(ig),zGauges(ig),sprintf('  G%d',selected(ig)), ...
            'FontSize',8,'VerticalAlignment','bottom');
    end
end
% SWL is drawn as a reference only, since it may vary over the time series.
swl0=input.bc(1,idx.swl);
plot([min(beachX) max(beachX)],[swl0 swl0],'--', ...
    'Color',[0.55 0.55 0.55],'LineWidth',1);
xlabel('Cross-shore distance (m)'); ylabel('z (m)');
title('Beach profile and measurement locations','FontWeight','normal');
legend('Beach profile','Gauges','SWL at first record', ...
    'Location','best','Box','off');
grid on; box on;

if options.saveFigures
    saveas(fig,fullfile(options.outputDir,'IG_all_gauges_with_profile.png'));
    savefig(fig,fullfile(options.outputDir,'IG_all_gauges_with_profile.fig'));
end
if ~options.showFigures, close(fig); end
end

function opts=optDefault(opts,key,val)
if ~isfield(opts,key) || isempty(opts.(key)), opts.(key)=val; end
end
function [r,b]=metrics(o,m)
r=NaN;b=NaN;
if isempty(o), return; end
r=sqrt(mean((m-o).^2));
if sum(o)>0, b=100*sum(m-o)/sum(o); end
end
