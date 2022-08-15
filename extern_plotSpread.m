function handles = extern_plotSpread(varargin)
%PLOTSPREAD plots distributions of points by spreading them around the y-axis
%
% SYNOPSIS: handles = plotSpread(data, propertyName, propertyValue, ...)
%           handles = plotSpread(ah, ...
%           deprecated:
%           handles = plotSpread(data,binWidth,spreadFcn,xNames,showMM,xValues)
%
% INPUT data: cell array of distributions or nDatapoints-by-mDistributions
%           array, or array with data that is indexed by either
%           distributionIdx or categoryIdx, or both.
%       distributionIdx: grouping variable that determines to which
%           distribution a data point belongs. Grouping is
%           resolved by calling grp2idx, and unless xNames have
%           been supplied, group names determine the x-labels.
%           If the grouping variable is numeric, group labels also
%           determine x-values, unless the parameter xValues has
%           been specified.
%       distributionColors : color identifier (string, cell array of
%           strings), or colormap, with a single color, or one color per
%           distribution (or per entry in distributionIdx). Colors the
%           distributions. Default: 'b'
%       distributionMarkers : string, or cell array of strings, with either
%           a single marker or one marker per distribution (or per entry in
%           distributionIdx). See linespec for admissible markers.
%           Default: '.'
%		categoryIdx: grouping variable that determines group membership for data
%			points across distributions. Grouping is resolved by calling
%           grp2idx.
%       categoryColors : color identifier (cell array of
%           strings), or colormap, with one color per category.
%           Colors the categories, and will override distributionColors.
%           Default is generated using distinguishable_colors by Timothy E.
%           Holy.
%       categoryMarkers : cell array of strings, with one marker per
%           category. See linespec for admissible markers. Will override
%           distributionMarkers. Default: ''
%       categoryLabels : cell array of strings with one label per category
%           (categories sorted in ascending order). Default: unique
%           category indices
%       binWidth : width of bins (along y) that control which data
%           points are considered close enough to be spread. Default: 0.1
%       markerSize: how big to make your marker ; 5 default
%       spreadFcn : cell array of length 2 with {name,param}
%           if name is 'lin', the spread goes linear with the number of
%             points inside the bin, until it reaches the maximum of 0.9 at
%             n==param.
%           if name is 'xp', the spread increases as 1-exp(log(0.9)*x).
%             param is empty
%           Default {'xp',[]}
%       spreadWidth : width, along the x-axis (y-axis if flipped) that can
%           at most be covered by the points. Default:
%           median(diff(sort(xValues))); 1 if no xValues have been supplied
%       showMM : if 1, mean and median are shown as red crosses and
%                green squares, respectively. Default: 0
%                2: only mean
%                3: only median
%                4: mean +/- standard error of the mean (no median)
%                5: mean +/- standard deviation (no median)
%       xNames : cell array of length nDistributions containing x-tick names
%               (instead of the default '1,2,3')
%       xValues : list of x-values at which the data should
%                 be plotted. Default: 1,2,3...
%       xMode  : if 'auto', x-ticks are spaced automatically. If 'manual',
%                there is a tick for each distribution. If xNames is
%                provided as input, xMode is forced to 'manual'. Default:
%                'manual'.
%       xyOri  : orientation of axes. Either 'normal' (=default), or
%                'flipped'. If 'flipped', the x-and y-axes are switched, so
%                that violin plots are horizontal. Consequently,
%                axes-specific properties, such as 'yLabel' are applied to
%                the other axis.
%       yLabel : string with label for y-axis. Default : ''
%       ah  : handles of axes into which to plot
%
% OUTPUT handles: 3-by-1 cell array with handles to distributions,
%          mean/median etc, and the axes, respectively
%
% REMARKS: plotSpread is useful for distributions with a small number of
%          data points. For larger amounts of data, distributionPlot is
%          more suited.
%
% EXAMPLES: data = {randn(25,1),randn(100,1),randn(300,1)};
%           figure,plotSpread(data,[],[],{'25 pts','100 pts','300 pts'})
%
%            data = [randn(50,1);randn(50,1)+3.5]*[1 1];
%            catIdx = [ones(50,1);zeros(50,1);randi([0,1],[100,1])];
%            figure
%            plotSpread(data,'categoryIdx',catIdx,...
%                 'categoryMarkers',{'o','+'},'categoryColors',{'r','b'})
%
% END
%
% created with MATLAB ver.: 7.9.0.3470 (R2009b) on Mac OS X  Version: 10.5.7 Build: 9J61
%
% created by: jonas
% DATE: 11-Jul-2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def.binWidth = 0.1;
def.spreadFcn = {'xp',[]};
def.xNames = [];
def.showMM = false;
def.xValues = [];
def.distributionIdx = [];
def.distributionColors = 'b';
def.distributionMarkers = '.';
def.xMode = 'manual';
def.xyOri = 'normal';
def.categoryIdx = [];
def.categoryColors = [];
def.categoryMarkers = '';
def.categoryLabels = '';
def.yLabel = '';
def.spreadWidth = [];
def.markerSize = 5;

% in development
def.individualLabels = false; % one category label across all distributions
%                               this should be smartly determined rather
%                               than hard-coded

%% CHECK INPUT

% check for axes handle
if ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')
    ah = varargin{1};
    data = varargin{2};
    varargin(1:2) = [];
    newAx = false;
else
    ah = gca;
    data = varargin{1};
    varargin(1) = [];
    % if the axes have children, it's not new (important for adjusting
    % limits below)
    newAx = isempty(get(ah,'Children'));
end

% optional arguments
parserObj = inputParser;
parserObj.FunctionName = 'plotSpread';
distributionIdx = [];distributionLabels = '';
if ~isempty(varargin) && ~ischar(varargin{1}) && ~isstruct(varargin{1})
    % old syntax
    parserObj.addOptional('binWidth',def.binWidth);
    parserObj.addOptional('markerSize',def.markerSize);
    parserObj.addOptional('spreadFcn',def.spreadFcn);
    parserObj.addOptional('xNames',def.xNames);
    parserObj.addOptional('showMM',def.showMM);
    parserObj.addOptional('xValues',def.xValues);
    
    parserObj.parse(varargin{:});
    opt = parserObj.Results;
    
    opt.distributionIdx = [];
    opt.distributionColors = def.distributionColors;
    opt.distributionMarkers = def.distributionMarkers;
    opt.xMode = def.xMode;
    opt.xyOri = def.xyOri;
    opt.categoryIdx = [];
    opt.categoryColors = def.distributionColors;
    opt.categoryMarkers = def.distributionMarkers;
    opt.yLabel = '';
    opt.spreadWidth = def.spreadWidth;
    opt.individualLabels = false;
    
    for fn = fieldnames(def)'
        if ~isfield(opt,fn{1})
            % Manually adding the new defaults means a lot fewer bugs
            error('please add option %s to old syntax',fn{1});
        end
        if isempty(opt.(fn{1}))
            opt.(fn{1}) = def.(fn{1});
        end
    end
    
else
    % new syntax
    defNames = fieldnames(def);
    for dn = defNames(:)'
        parserObj.addParamValue(dn{1},def.(dn{1}));
    end
    
    
    parserObj.parse(varargin{:});
    opt = parserObj.Results;
end

% We want data to be a vector, so that indexing with both groupIdx and
% distributionIdx becomes straightforward, and so that we can conveniently
% eliminate NaNs that otherwise could mess up grouping.
% Consequently, if data is a cell array, we convert it, and build a
% corresponding distributionIdx (allowing a user-supplied distributionIdx
% to override, though), and then we go and take care of groupIdx. Once all
% three indices have been built, NaN can be removed.

if iscell(data)
    % make sure data is all n-by-1
    data = cellfun(@(x)x(:),data,'UniformOutput',false);
    nData = length(data);
    nn = cellfun(@numel,data);
    % make vector
    data = cat(1,data{:});
    distributionIdx = repeatEntries((1:nData)',nn);
else
    % distributions in columns
    nData = size(data,2);
    distributionIdx = repeatEntries((1:nData)',size(data,1));
    data = data(:);
end



% distribution groups
if ~isempty(opt.distributionIdx)
    [distributionIdx,distributionLabels,vals] = grp2idx(opt.distributionIdx);
    % convert data to cell array
    nData = length(distributionLabels);
    % if not otherwise provided, use group labels for xnames
    if isempty(opt.xNames)
        opt.xNames = distributionLabels;
        if ~iscell(opt.xNames)
            opt.xNames = num2cell(opt.xNames);
        end
    end
    if isnumeric(vals) && isempty(opt.xValues)
        opt.xValues = vals;
    end
end

if ~isempty(opt.xNames)
    opt.xMode = 'manual';
end


% distribution colors&markers
if ischar(opt.distributionColors)
    opt.distributionColors = {opt.distributionColors};
end
if iscell(opt.distributionColors)
    if length(opt.distributionColors) == 1
        % expand
        opt.distributionColors = repmat(opt.distributionColors,nData,1);
    elseif length(opt.distributionColors) ~= nData
        error('please submit one color per distribution (%i dist, %i colors)',nData,length(opt.distributionColors));
    end
    
else
    if size(opt.distributionColors,2) ~= 3
        error('please specify colormap with three columns')
    end
    if size(opt.distributionColors,1) == 1
        opt.distributionColors = repmat(opt.distributionColors,nData,1);
    elseif size(opt.distributionColors,1) ~= nData
        error('please submit one color per distribution (%i dist, %i colors)',nData,size(opt.distributionColors,1));
    end
    
    % create a cell array
    opt.distributionColors = mat2cell(opt.distributionColors,ones(nData,1),3);
end

if ischar(opt.distributionMarkers)
    opt.distributionMarkers = {opt.distributionMarkers};
end
if length(opt.distributionMarkers) == 1
    % expand
    opt.distributionMarkers = repmat(opt.distributionMarkers,nData,1);
elseif length(opt.distributionMarkers) ~= nData
    error('please submit one color per distribution (%i dist, %i colors)',nData,length(opt.distributionMarkers));
end


stdWidth = 1;
if isempty(opt.xValues)
    opt.xValues = 1:nData;
end


if isempty(opt.spreadWidth) 
    % scale width
    tmp = median(diff(sort(opt.xValues)));
    if ~isnan(tmp)
        stdWidth = tmp;
    end
else
    stdWidth = opt.spreadWidth;
end

if ~ischar(opt.xyOri) || ~any(ismember(opt.xyOri,{'normal','flipped'}))
    error('option xyOri must be either ''normal'' or ''flipped'' (is ''%s'')',opt.xyOri);
end


% check for categoryIdx/colors/markers
% If there are categories, check colors/markers individually first,
% then check whether any of them at all have been supplied, and
% if not, override distributionColors with default categoryColors

if isempty(opt.categoryIdx)
    categoryIdx = ones(size(distributionIdx));
    nCategories = 1;
    categoryLabels = '';
else
    [categoryIdx,categoryLabels] = grp2idx(opt.categoryIdx(:));
    nCategories = max(categoryIdx);
end
if ~isempty(opt.categoryLabels)
    categoryLabels = opt.categoryLabels;
elseif ~iscell(categoryLabels)
    categoryLabels = num2cell(categoryLabels);
end

% plotColors, plotMarkers, plotLabels: nDist-by-nCat arrays
plotColors = repmat(opt.distributionColors(:),1,nCategories);
plotMarkers= repmat(opt.distributionMarkers(:),1,nCategories);

if isempty(distributionLabels)
    distributionLabels = opt.xNames;
    if isempty(distributionLabels)
        distributionLabels = cellstr(num2str(opt.xValues(:)));
    end
end

if nCategories == 1
    plotLabels = distributionLabels(:);
else
    plotLabels = cell(nData,nCategories);
    for iData = 1:nData
        for iCategory = 1:nCategories
            if opt.individualLabels
            plotLabels{iData,iCategory} = ...
                sprintf('%s-%s',num2str(distributionLabels{iData}),...
                num2str(categoryLabels{iCategory}));
            else
                plotLabels{iData,iCategory} = ...
                sprintf('%s',...
                num2str(categoryLabels{iCategory}));
            end
        end
    end
    
end




categoryIsLabeled = false;
if nCategories > 1
    % if not using defaults for categoryColors: apply them
    if ~any(strcmp('categoryColors',parserObj.UsingDefaults))
        if iscell(opt.categoryColors)
            if length(opt.categoryColors) ~= nCategories
                error('please supply one category color per category')
            end
            plotColors = repmat(opt.categoryColors(:)',nData,1);
            categoryIsLabeled = true;
        else
            if all(size(opt.categoryColors) ~= [nCategories,3])
                error('please supply a #-of-categories-by-3 color array')
            end
            plotColors = repmat( mat2cell(opt.categoryColors,ones(nCategories,1),3)', nData,1);
            categoryIsLabeled = true;
        end
    end
    
    if ~any(strcmp('categoryMarkers',parserObj.UsingDefaults))
        if length(opt.categoryMarkers) ~= nCategories
            error('please supply one category marker per category')
        end
        if ~iscell(opt.categoryMarkers)
            error('please supply a list of markers as cell array')
        end
        plotMarkers = repmat(opt.categoryMarkers(:)',nData,1);
        categoryIsLabeled = true;
    end
    
    if ~categoryIsLabeled
        % use distinguishable_colors to mark categories
        
        plotColors = repmat( mat2cell(...
            distinguishable_colors(nCategories),...
            ones(nCategories,1),3)', nData,1);
        
    end
    
end


% remove NaNs from data
badData = ~isfinite(data) | ~isfinite(distributionIdx) | ~isfinite(categoryIdx);
data(badData) = [];
distributionIdx(badData) = [];
categoryIdx(badData) = [];




%% TRANSFORM DATA
% Here, I try to estimate what the aspect ratio of the data is going to be
fh = figure('Visible','off');
if ~isempty(data)
    minMax = [min(data);max(data)];
else
    minMax = [0 1];
end
switch opt.xyOri
    case 'normal'
        plot([0.5;nData+0.5],minMax,'o');
    case 'flipped'
        plot(minMax,[0.5;nData+0.5],'o');
        
end
aspectRatio = get(gca,'DataAspectRatio');
close(fh);

tFact = aspectRatio(2)/aspectRatio(1);
if strcmp(opt.xyOri,'flipped')
    tFact = 1/tFact;
end

%% SPREAD POINTS
% assign either nData, or xValues number of values, in case we're working
% with group-indices
[m,md,sem,sd] = deal(nan(max(nData,length(opt.xValues)),1));
% make sure xValues are not something weird
opt.xValues = double(opt.xValues);

    
% augment data to make n-by-2
data(:,2) = 0;
for iData = 1:nData
    currentDataIdx = distributionIdx==iData;
    currentData = data(currentDataIdx,1);
    
    if ~isempty(currentData)
        
        % transform and sort
        currentData = currentData / tFact;
        %currentData = sort(currentData);
        
        % add x
        currentData = [ones(size(currentData))*opt.xValues(iData),currentData]; %#ok<AGROW>
        
        % step through the data in 0.1 increments. If there are multiple
        % entries, spread along x
        for y = min(currentData(:,2)):opt.binWidth:max(currentData(:,2))
            % find values
            valIdx = find(currentData(:,2) >= y & currentData(:,2) < y+opt.binWidth);
            nVal = length(valIdx);
            if nVal > 1
                % spread
                switch opt.spreadFcn{1}
                    case 'xp'
                        spreadWidth = stdWidth*0.9*(1-exp(log(0.9)*(nVal-1)));
                    case 'lin'
                        spreadWidth = stdWidth*0.9*min(nVal-1,opt.spreadFcn{2})/opt.spreadFcn{2};
                end
                spreadDist = spreadWidth / (nVal - 1);
                if isEven(nVal)
                    offset = spreadDist / 2;
                else
                    offset = eps;
                end
                for v = 1:nVal
                    currentData(valIdx(v),1) = opt.xValues(iData) + offset;
                    % update offset
                    offset = offset - sign(offset) * spreadDist * v;
                end
            end
        end
        
        % update data
        currentData(:,2) = data(currentDataIdx,1);
        data(currentDataIdx,:) = currentData;
        
        
        if opt.showMM > 0
            m(iData) = nanmean(currentData(:,2));
            md(iData) = nanmedian(currentData(:,2));
            sd(iData) = nanstd(currentData(:,2));
            sem(iData) = sd(iData)/sqrt(sum(isfinite(currentData(:,2))));
        end
    end % test isempty
end


%% plot
edgeColor = [1 1 1]*0.9;
set(ah,'NextPlot','add')
ph = NaN(nData,nCategories);
for iData = 1:nData
    for iCategory = 1:nCategories
        currentIdx = distributionIdx == iData & categoryIdx == iCategory;
        if any(currentIdx)
            switch opt.xyOri
                case 'normal'
                    ph(iData,iCategory) = plot(ah,data(currentIdx,1),...
                        data(currentIdx,2),...
                        'marker',plotMarkers{iData,iCategory},...
                        'MarkerFaceColor',plotColors{iData,iCategory},...
                        'Color',edgeColor, ...
                        'lineStyle','none',...
                        'MarkerSize',opt.markerSize, ...
                        'DisplayName',plotLabels{iData,iCategory});
                case 'flipped'
                    ph(iData,iCategory) = plot(ah,data(currentIdx,2),...
                        data(currentIdx,1),...
                        'marker',plotMarkers{iData,iCategory},...
                        'MarkerFaceColor',plotColors{iData,iCategory},...
                        'Color',edgeColor, ...
                        'MarkerSize', opt.markerSize, ...
                        'lineStyle','none',...
                        'DisplayName',plotLabels{iData,iCategory});
            end
        end
    end
end



% if ~empty, use xNames
switch opt.xyOri
    case 'normal'
        switch opt.xMode
            case 'manual'
                set(ah,'XTick',opt.xValues);
                if ~isempty(opt.xNames)
                    set(ah,'XTickLabel',opt.xNames)
                end
            case 'auto'
                % no need to do anything
        end
        
        % have plot start/end properly
        minX = min(opt.xValues)-stdWidth;
        maxX = max(opt.xValues)+stdWidth;
        if ~newAx
            oldLim = xlim;
            minX = min(minX,oldLim(1));
            maxX = max(maxX,oldLim(2));
        end
        xlim([minX,maxX])
        
        ylabel(ah,opt.yLabel)
        
    case 'flipped'
        switch opt.xMode
            case 'manual'
                set(ah,'YTick',opt.xValues);
                if ~isempty(opt.xNames)
                    set(ah,'YTickLabel',opt.xNames)
                end
            case 'auto'
                % no need to do anything
        end
        
        % have plot start/end properly (for ease of copying, only switch
        % xlim to ylim
        minX = min(opt.xValues)-stdWidth;
        maxX = max(opt.xValues)+stdWidth;
        if ~newAx
            oldLim = ylim;
            minX = min(minX,oldLim(1));
            maxX = max(maxX,oldLim(2));
        end
        ylim([minX,maxX])
        
        xlabel(ah,opt.yLabel);
        
end

% ## in development
if ~opt.individualLabels
       % hack: add legend entry only once per category
       goodH = ishandle(ph);
       for iCategory = 1:nCategories
           for iData = find(goodH(:,iCategory),1,'first')+1:nData
       if goodH(iData,iCategory)
           set(get(get(ph(iData,iCategory),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
       end
           end
       end
       
end


% add mean/median
mh = [];mdh=[];
if opt.showMM
    % plot mean, median. Mean is filled red circle, median is green square
    % I don't know of a very clever way to flip xy and keep everything
    % readable, thus it'll be copy-paste
    switch opt.xyOri
        case 'normal'
            if any(opt.showMM==[1,2])
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
            end
            if any(opt.showMM==[1,3])
%                mdh = plot(ah,opt.xValues,md,'sg','MarkerSize',12);
                for ii=1:length(md)
                    mdh(ii) = plot(ah,opt.xValues(ii)+[-0.25 0.25],md(ii)*[1 1],'k-','LineWidth', 5);
                end
            end
            if opt.showMM == 4
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,opt.xValues,m,sem);
            end
            if opt.showMM == 5
                mh = plot(ah,opt.xValues,m,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,opt.xValues,m,sd);
            end
        case 'flipped'
            if any(opt.showMM==[1,2])
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
            end
            if any(opt.showMM==[1,3])
                mdh = plot(ah,md,opt.xValues,'sg','MarkerSize',12);
            end
            if opt.showMM == 4
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,m,opt.xValues,[sem,NaN(size(sem))]);
            end
            if opt.showMM == 5
                mh = plot(ah,m,opt.xValues,'+r','Color','r','MarkerSize',12);
                mdh = myErrorbar(ah,m,opt.xValues,[sd,NaN(size(sd))]);
            end
    end
end

%==========================
%% CLEANUP & ASSIGN OUTPUT
%==========================

if nargout > 0
    handles{1} = ph;
    handles{2} = [mh;mdh];
    handles{3} = ah;
end



function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);

% Copyright 2010-2011 by Timothy E. Holy

  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end

  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [dummy,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end



function out = isEven(in)
%ISEVEN checks whether a number is even
%
% SYNOPSIS out = isEven(in)
%
% INPUT    in :  input (array) of numbers to be tested. 
% OUTPUT   out:  array of size(in) with 
%                   1 for even integers and zero
%                   0 for odd integers
%                 NaN for non-integers
%                out is a logical array as long as the input is all integers.
%
% c: jonas 5/05
% Last modified 11/24/2009 - Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = mod(in+1, 2);
% Set NaN for non-integer data, because they are neither odd or even
out((out ~= 0) & (out ~= 1)) = NaN;

% since doubles cannot be used for logical indexing, we should convert to
% logicals if possible. 
if all(isfinite(out(:)))
    out = logical(out);
end

function hh = myErrorbar(varargin)
%MYERRORBAR Adds errorbars to existing plot (unlike errorbar.m, which creates a new plot, and allows only bars for y values)
%   MYERRORBAR(X,Y,L,U) adds error bars to the graph of vector X vs. vector Y with
%   error bars specified by the vectors L and U.  L and U contain the
%   lower and upper error ranges for each point in Y.  Each error bar
%   is L(i) + U(i) long and is drawn a distance of U(i) above and L(i)
%   below the points in (X,Y). If X,Y,L and U are matrices then each column
%   produces a separate line.
%   If L,U are the same size as X, Y, only error bars for Y will be plotted.
%   If L,U are twice the size of X,Y (or have twice the number of columns for
%   matrices), the first half of L, U specifies error bar lengths for X and the
%   second half specifies error bars for Y
%
%   MYERRORBAR(X,Y,E) or MYERRORBAR(Y,E) plots error bars [Y-E Y+E].
%
%   MYERRORBAR(AX,...), where AX is an axis handle, plots errorbars into
%                       axes AX
%
%   H = MYERRORBAR(...) returns a vector of line handles.
%
%   The tag of the errorbar-lines is: errorBar
%
%   For example,
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      myErrorbar(x,y,e)
%   draws symmetric error bars of unit standard deviation for y values.
%      myErrorbar(x,y,[e,e])
%   draws symmetric error bars of unit standard deviation for x and y
%   values.
%
%   Based on the matlab-function errorbar as revised by Claude Berney
%   c: jonas, 06-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==================
% check input
%==================

if nargin < 2
    error('not enough input arguments!')
end

% check if the first input argument is a handle
if length(varargin{1}) == 1 && ishandle(varargin{1}) && strcmpi(get(varargin{1},'Type'),'axes')
    axesH = varargin{1};
    % remove axis handle
    varargin(1) = [];
else
    axesH = gca;
end

% there could be
% y,e
% x,y,e
% x,y,l,u

switch length(varargin)
    case 2
        % y, e
        y = varargin{1};
        y = y(:);
        lengthY = length(y);
        x = [1:lengthY]';
        
        e = varargin{2};
        % check for 2 dimension errorbars
        e = e(:);
        if length(e) == 2*lengthY
            e = reshape(e,lengthY,2);
        end
        [l,u] = deal(e);
        
    case 3
        % x,y,e
        x = varargin{1};
        x = x(:);
        y = varargin{2};
        y = y(:);
        lengthY = length(y);
        
        e = varargin{3};
        % check for 2 dimension errorbars
        e = e(:);
        if length(e) == 2*lengthY
            e = reshape(e,lengthY,2);
        end
        [l,u] = deal(e);
        
    case 4
        % x,y,l,u
        % x,y,e
        x = varargin{1};
        x = x(:);
        y = varargin{2};
        y = y(:);
        lengthY = length(y);
        
        l = varargin{3};
        % check for 2 dimension errorbars
        l = l(:);
        if length(l) == 2*lengthY
            l = reshape(l,lengthY,2);
        end
        u = varargin{4};
        % check for 2 dimension errorbars
        u = u(:);
        if length(u) == 2*lengthY
            u = reshape(u,lengthY,2);
        end
        
        if ~all(size(u)==size(l))
            error('l, u have to be the same size!')
        end
        
end % switch number of inputs


u = abs(u);
l = abs(l);

if ischar(x) || ischar(y) || ischar(u) || ischar(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y))
    error('The sizes of X and Y must be the same.');
end

if isequal([1 2].*size(x),size(l)) && isequal([1 2].*size(x),size(u))
    xyBars = 1;
elseif isequal(size(x),size(l)) && isequal(size(x),size(u))
    xyBars = 0;
else
    error('The sizes of L and U must be equal to or twice the size of X, Y')
end

%=======================


% Plot graph and bars
hold_state = ishold;
hold on;


%find color of current plot
dataH = get(axesH,'Children');
myLineH = dataH(1);
% support also bar plots
if strcmp(get(myLineH,'Type'),'hggroup')
    latestColor = get(myLineH,'EdgeColor'); %new children are added on top!
else
    latestColor = get(myLineH,'Color'); %new children are added on top!
end

tee=0;
if ~strcmp('log',get(axesH,'XScale'))
    tee = (max(x(:))-min(x(:)))/100;  % make tee .02 x-distance for error bars
    tee = min(tee,0.3*nanmedian(diff(unique(x(:))))); % or at most 0.3*deltaX
    xl = x - tee;
    xr = x + tee;
end
if strcmp('log',get(axesH,'XScale'))
    tee = (max(log(x(:)))-min(log(x(:))))/100;  % make tee .02 x-distance for error bars
    tee = min(tee,0.3*nanmedian(diff(unique(log(x(:)))))); % or at most 0.3*deltaX
    
    xl = x *exp(tee);
    xr = x *exp(-tee);
end

if xyBars
    if ~strcmp('log',get(axesH,'YScale'))
        tee = (max(y(:))-min(y(:)))/100;  % make tee .02 y-distance for error bars
        tee = min(tee,0.3*nanmedian(diff(unique(y(:))))); % or at most 0.3*deltaY
        
        yl = y - tee;
        yr = y + tee;
    end
    if strcmp('log',get(axesH,'YScale'))
        tee = (max(log(y(:)))-min(log(y(:))))/100;  % make tee .02 y-distance for error bars
        tee = min(tee,0.3*nanmedian(diff(unique(log(y(:)))))); % or at most 0.3*deltaX
        
        yl = y *exp(tee);
        yr = y *exp(-tee);
    end
end

%specify coordinates to plot error bars
if xyBars
    xtop = x + u(:,1:size(x,2));
    xbot = x - l(:,1:size(x,2));
    ytop = y + u(:,size(x,2)+1:end);
    ybot = y - l(:,size(x,2)+1:end);
else
    ytop = y + u;
    ybot = y - l;
end
n = size(y,2);

% build up nan-separated vector for bars
xb = zeros(lengthY*9,n);
xb(1:9:end,:) = x;
xb(2:9:end,:) = x;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xl;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(lengthY*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = ytop;
yb(5:9:end,:) = ytop;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ybot;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;

h = [line(xb,yb,'parent',axesH,'Color',latestColor)];

if xyBars
    
    xb(1:9:end,:) = xtop;
    xb(2:9:end,:) = xbot;
    xb(3:9:end,:) = NaN;
    xb(4:9:end,:) = xtop;
    xb(5:9:end,:) = xtop;
    xb(6:9:end,:) = NaN;
    xb(7:9:end,:) = xbot;
    xb(8:9:end,:) = xbot;
    xb(9:9:end,:) = NaN;
    
    yb(1:9:end,:) = y;
    yb(2:9:end,:) = y;
    yb(3:9:end,:) = NaN;
    yb(4:9:end,:) = yl;
    yb(5:9:end,:) = yr;
    yb(6:9:end,:) = NaN;
    yb(7:9:end,:) = yl;
    yb(8:9:end,:) = yr;
    yb(9:9:end,:) = NaN;
    
    h = [h;line(xb,yb,'parent',axesH,'Color',latestColor)];
    
end

%set the tag of all errorBar-objects to 'errorBar'
set(h,'Tag','errorBar');

% make sure errorbar doesn't produce a legend entry
for lineH = h'
    set(get(get(lineH,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
end


if ~hold_state, hold off; end

if nargout>0, hh = h; end


function out = repeatEntries(val,kTimes)
%REPEATENTRIES fills a matrix with k repeats the rows of the input matrix
%
% SYNOPSIS out = repeatEntries(val,kTimes)
%
% INPUT    val    : matrix (or vectors) containing the rows to repeat (works for strings, too)
%          kTimes : number of repeats of each row (scalar or vector of size(vlaues,1))
%
% OUTPUT   out    : matrix of size [sum(kTimes) size(values,2)] containing
%                   repeated entries specified with k
%
% EXAMPLES     repeatEntries([1;2;3;4],[2;3;1;1]) returns [1;1;2;2;2;3;4]
%
%              repeatEntries([1;2;3;4],2) returns [1;1;2;2;3;3;4;4]
%
% c: jonas, 2/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note: in case we need to speed this up: adapt the code below
% nn = cellfun(@numel,points);
% a = find(nn);
% index = zeros(sum(nn),1);
% index([1;cumsum(nn(a(1:end-1)))+1])=1;
% 
% % get the indices
% ii = a(cumsum(index));

%===========
% test input
%===========

% nargin
if nargin ~= 2 || isempty(val) || isempty(kTimes)
    error('two non-empty input arguments are needed!')
end

% size
valSize = size(val);
if length(valSize)>2
    error('only 2D arrays supported for val')
end



% decide whether we have scalar k
numK = length(kTimes);
if numK == 1
    scalarK = 1;
elseif numK ~= valSize(1)
    error('vector k must have the same length as the number of rows in val or be a scalar')
else
    % check again whether we could use scalar k
    if all(kTimes(1) == kTimes)
        scalarK = 1;
        kTimes = kTimes(1);
    else
        scalarK = 0;
    end
end

% do not care about size of k: we want to make a col vector out of it - and
% this vector should only contain nonzero positive integers
kTimes = round(kTimes(:));
% if there are any negative values or zeros, remove the entry
if scalarK && kTimes < 1
    out = [];
    return
end
if ~scalarK
    badK = kTimes < 1;
    kTimes(badK) = [];
    val(badK,:) = [];
    % update valSize
    valSize = size(val);
    if any(valSize==0)
        out = [];
        return
    end
end
%kTimes = max(kTimes,ones(size(kTimes)));


%============
% fill in out
%============

% first the elegant case: scalar k
if scalarK

    % build repeat index matrix idxMat
    idxMat = meshgrid( 1:valSize(1), 1:kTimes(1) );
    idxMat = idxMat(:); % returns [1;1...2;2;... etc]

    out = val(idxMat,:);

    % second: the loop
else

    % init out, init counter
    if iscell(val)
        out = cell(sum(kTimes) , valSize(2));
    else
    out = zeros( sum(kTimes), valSize(2) );
    end
    endct = 0;

    if valSize(2) == 1

        % vector: fill directly

        % loop and fill
        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            out(startct:endct,:) = val(i);
        end % for i=1:valSize(1)

    else

        % matrix: fill via index list

        idxMat = zeros(sum(kTimes),1);

        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            idxMat(startct:endct) = i;
        end % for i=1:valSize(1)
        out = val(idxMat,:);

    end

    % check for strings and transform if necessary
    if ischar(val)
        out = char(out);
    end

end % if doScalar


