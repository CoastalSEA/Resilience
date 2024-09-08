function resilience_plotting(isdefault)
%
%-------header-------------------------------------------------------------
% NAME
%   resilience_plotting.m
% PURPOSE
%   function to analyse and plot O'Neill et al, 2018 data and Fanning et
%   al, 2021 data.
% USAGE
%   resilience_plotting(isdefault);
% NOTES 
%   O'Neill D W, Fanning A L, Lamb W F and Steinberger J K, 2018, 
%   A good life for all within planetary boundaries. Nature Sustainability,
%   https://doi.org/10.1038/s41893-018-0021-4. Data is from SI.
%
%   Fanning A L, O’Neill D W, Hickel J and Roux N, 2021, 
%   The social shortfall and ecological overshoot of nations.
%   https://doi.org/10.1038/s41893-021-00799-z.
%
%   Gapminder, 2022, Gapminder - Population v7 and other sources – with 
%   major processing by Our World in Data, Our World in Data, 
%   https://ourworldindata.org/grapher/population.
% SEE ALSO
%   scripts to read Excel files and save as mat files
%   read_oneill_data.m, read_fanning_data.m, read_pop_data.m
%
% Author: Ian Townend
% CoastalSEA (c) Apr 2024
%--------------------------------------------------------------------------
%
    userprompt = 'Select mat file>';
    [~, path]=uigetfile('*.mat',userprompt,'MultiSelect','off');
    if path==0, return; end    
    %path = 'D:\OneDrive\Work\Research\Resilience\';

    listxt = {'Global radar plot','Country radar plot','Exceedance plot',...
                'Scaled data (scores)','Weighted scores','Country scores & weights',...
                'Index timeseries'};
    selection = listdlg('ListString',listxt,"PromptString",'Select plot type:',...
                        'Name','Plot type','SelectionMode','single','ListSize',[180,200]);
    if selection==length(listxt)
        %load Fanning data with struct of data tables from 1992-2015
        load([path,'FanningData'],'BioPhysData','SocialData','-mat')
        load([path,'PopData1900-2021'],'PopData','-mat')   
    else      
        %load O'Neill data tables for 2011
        load([path,'ONeillData'],'BioPhysData','SocialData','-mat')
        load([path,'PopData1900-2021'],'PopData','-mat')    
        YrPopData = PopData.yr2011;
    end
    clear path

    switch selection
        case 1
            global_radar_plot(BioPhysData,SocialData);
        case 2 
            country_radar_plot(BioPhysData,SocialData);
        case 3
            [BioPhysData,SocialData,YrPopData] = sortPop(BioPhysData,SocialData,YrPopData);
            exceedance_plot(BioPhysData,SocialData,YrPopData);
        case 4
            [BioPhysData,SocialData,YrPopData] = sortPop(BioPhysData,SocialData,YrPopData);
            scaled_data(BioPhysData,SocialData,YrPopData,true);
        case 5
            [BioPhysData,SocialData,YrPopData] = sortPop(BioPhysData,SocialData,YrPopData);
            weighted_scores(BioPhysData,SocialData,YrPopData);
        case 6
            [BioPhysData,SocialData,YrPopData] = sortPop(BioPhysData,SocialData,YrPopData);
            country_selection(BioPhysData,SocialData,YrPopData);
        case 7
            index_ts(BioPhysData,SocialData,PopData);
    end
end
%--------------------------------------------------------------------------
% Radar plots
%% ------------------------------------------------------------------------
function global_radar_plot(BioPhysData,SocialData)
    %create radar plots for the mean of the biophysical and social data
    meanfunc = @(x) mean(x,'omitnan');
    AllBioPhys = varfun(meanfunc,BioPhysData, 'InputVariables', @isnumeric);
    AllBioPhys.Properties.VariableDescriptions = BioPhysData.Properties.VariableDescriptions ;
    AllSocial = varfun(meanfunc,SocialData, 'InputVariables', @isnumeric);
    AllSocial.Properties.VariableDescriptions = SocialData.Properties.VariableDescriptions ;
    Nrec = height(BioPhysData);

    titletxt = sprintf('Metric mean values for %d countires',Nrec);
    radar_plot(AllBioPhys,AllSocial,titletxt);
end

%%
function country_radar_plot(BioPhysData,SocialData)
    %create radar plots for the biophysical and social data sets by country  
    Country = BioPhysData.Properties.RowNames;
    ok = 1;
    while ok>0
        selection = listdlg("PromptString",'Select a country:',...
                               'SelectionMode','single','ListString',Country);    
        if isempty(selection), ok=0; continue; end
        SelBioPhys = BioPhysData(selection,:);
        SelSocial = SocialData(selection,:);

        titletxt = sprintf('Metric values for %s',Country{selection});
        radar_plot(SelBioPhys,SelSocial,titletxt);      
    end
end

%%
function radar_plot(BioPhysData,SocialData,titletxt)
    %construct radar plot of biophysical and social data
    hf = figure('Tag','PlotFig');
    ax = axes(hf); 
    blue = mcolor('dark blue');
    green = mcolor('green'); 

    % USAGE: varargout = spider_plot(P, varargin)
    %P - The data points used to plot the spider chart. The rows are the 
    %    groups of data and the columns are the data points. The axes labels 
    %    and axes limits are automatically generated if not specified.
    %varargout - Figure handle of spider plot
    %varargin - Name-Value pairs

    %extract and plot BioPhysical data set
    dataset = BioPhysData{1,:};
    dataset(dataset<0.01 | isnan(dataset)) = 0.01;
    axlabels = BioPhysData.Properties.VariableDescriptions;
    axlim = ceil(max(dataset));
    axint = max(1,floor(axlim));
    nrec = length(dataset);    
    dataset = [ones(size(dataset));dataset];  
      
    s1 = subplot(1,2,1,ax);
    spider_plot(dataset,'AxesHandle',s1,...
                        'AxesLabels',axlabels,...
                        'AxesInterval',axint,...
                        'AxesLimits',repmat([0;axlim],1,nrec),...
                        'AxesOffset',0,...
                        'AxesPrecision', 1,...
                        'AxesDisplay', 'one',...    %limit axes labels to one axis
                        'AxesLabelsEdge','none',... %border for axis labels off
                        'Color',[green;blue],...
                        'Marker','none',...
                        'FillOption','on',...
                        'FillTransparency',0.5);
    hl1 = legend('BioPhysical','Safe space','Location','southwest');
    hl1.Position(1:2) = [0.05,0.1];
    clear dataset axlabels axlim nrec

    %extract and plot Social data set
    dataset = SocialData{1,:};
    dataset(dataset<0.01 | isnan(dataset)) = 0.01;
    axlabels = SocialData.Properties.VariableDescriptions;
    % axlim = ceil(max(dataset));
    % axint = max(1,floor(axlim));    
    axlim = 1.2; axint = 5;                               %bespoke settings
    nrec = length(dataset);    
    dataset = [ones(size(dataset));dataset];

    axes_shaded_limits = {repmat([1;axlim],1,nrec)};

    s2 = subplot(1,2,2);
    spider_plot(dataset,'AxesHandle',s2,...
                        'AxesLabels',axlabels,...
                        'AxesInterval',axint,...
                        'AxesLimits',repmat([0;axlim],1,nrec),... 'AxesTickLabels',axticklabels,...                       
                        'AxesOffset',0,...
                        'AxesPrecision', 2,...
                        'AxesDisplay', 'one',...    %limit axes labels to one axis
                        'AxesLabelsEdge','none',... %border for axis labels off
                        'Color',[green;blue],...
                        'Marker','none',...
                        'FillOption','on',...
                        'FillTransparency',0.5,...                        
                        'AxesShaded', 'on',...
                        'AxesShadedLimits', axes_shaded_limits,...
                        'AxesShadedColor', {green},...
                        'AxesShadedTransparency', 0.5);
    
    hf = findobj(s2.Children,'Type','Patch','FaceColor',green);
    idx = cellfun(@numel,{hf(:).Vertices},'UniformOutput',false);
    idx = [idx{:}]==22;
    hf(idx).FaceColor = [1,1,1];

    hl2 = legend('Social','Safe space','Location','southwest');
    hl2.Position(1:2) = [0.55,0.1];
    sgtitle(titletxt);
end

%--------------------------------------------------------------------------
% Exceedance plots
%% ------------------------------------------------------------------------
function exceedance_plot(BioPhysData,SocialData,YrPopData)
    %create plot to reproduce variant of Figure 2 in paper
    blue = mcolor('dark blue');
    %get the number of exceedances for each country
    BPexcess = sum(BioPhysData{:,:}>=1,2); %remove population data
    SLexcess = sum(SocialData{:,:}<=1,2);    
    
    %plot figure
    hf = figure('Tag','PlotFig');
    ax = axes(hf);
    markersz = 1000*YrPopData.Population/max(YrPopData.Population);
    scatter(ax,BPexcess,SLexcess,markersz,'filled',...
        'MarkerEdgeColor',blue,'MarkerFaceAlpha',0.1)
    xlabel('Number of BioPhysical Limits exceeded')
    ylabel('Number of Social Limits not met')
    title('Number of Countries missing limits')
    subtitle('Circles represent population size')
    ax.XLim(1) = -0.5;
    ax.YLim(1) = -0.5;
    ax.YTick = 0:1:12; 

    %plot histogram
    figure('Tag','PlotFig');
    tiledlayout(3,2);
    nexttile([2,2]);
    ylabels = flipud(cellstr(num2str((0:width(SocialData))')));
    tbl = table(BPexcess,SLexcess,'VariableNames',{'BioPhysical','Social'});
    hh1 = heatmap(tbl,1,2,'ColorMethod','count','YDisplayData',ylabels);
    hh1.XLabel = 'Number of BioPhysical Limits exceeded';
    hh1.YLabel = 'Number of Social Limits not met';
    hh1.Title = 'Number of Countries missing limits';

    nexttile
    histogram(SLexcess);
    xticks(0:width(SocialData));
    xticklabels(flipud(ylabels));
    xlabel('No. Social Limits not met');
    ylabel('No. Countries');

    nexttile    
    histogram(BPexcess);
    xticks(0:width(BioPhysData));
    xticklabels(cellstr(num2str((0:width(BioPhysData))'))');
    xlabel('No. BioPhysical Limits exceeded');
    ylabel('No. Countries');
end

%--------------------------------------------------------------------------
% Scores and weights
%% ------------------------------------------------------------------------
function [BPscores,SLscores,options] = scaled_data(BioPhysData,SocialData,YrPopData,isplot,options)
    %scale data to range 0-1 based on threshold and max-min values
    if nargin<5 || ~(any(strcmp({'Logistic','Linear segments','Linear'},options.type)))
        options = selectScoreScale();
    end
                 
    BPscores = scaledataset(BioPhysData,false,options);%biophysical variables scale negatively
    SLscores = scaledataset(SocialData,true,options); %social variables scale positively
    if isplot
        nBPvar = width(BPscores);
        ok = 1;
        listxt = [BPscores.Properties.VariableDescriptions,SLscores.Properties.VariableDescriptions];
        promptxt = 'Select individual variable to plot:';
        while ok>0
            selection = listdlg('ListString',listxt,"PromptString",promptxt,...
                            'Name','Plot type','SelectionMode','single');
            if isempty(selection), ok=0; continue; end
            %plot figure
            hf = figure('Tag','PlotFig');
            ax = axes(hf); %#ok<LAXES> 
            if selection<nBPvar+1 %biophys variable
                xvar = BioPhysData.(selection);
                yvar = BPscores.(selection);
            else
                xvar = SocialData.(selection-nBPvar);
                yvar = SLscores.(selection-nBPvar);
            end
        
            blue = mcolor('dark blue');
            markersz = 1000*YrPopData.Population/max(YrPopData.Population);
            scatter(ax,xvar,yvar,markersz,'filled',...
                            'MarkerEdgeColor',blue,'MarkerFaceAlpha',0.1)
            xlabel('Variable')
            ylabel('Score')
            title(sprintf('Score for %s',listxt{selection}))
            subtitle('Circles represent population size')
        end
        
        thumbnailplot(BioPhysData,SocialData,BPscores,SLscores,YrPopData,options,true);
    end
end

%%
function datascaled = scaledataset(dataset,ispos,options)
    %rescale the biophysical or social datasets
    bpmax = varfun(@max,dataset);
    bpmin = varfun(@min,dataset);
    nvar = width(dataset);   
    datascaled = dataset;
    for i=1:nvar
        datascaled{:,i} = scalevar(dataset.(i),bpmin.(i),bpmax.(i),ispos,options);
    end
end

%%
function yvar = scalevar(var,xmin,xmax,ispos,options)
    %variable rescaling function on to 0-1 score
    %variables have alreand been scaled such that x0=1
    x1 = 0.8;
    if strcmp(options.type,'Logistic')        
        if ispos
            params = struct('A',0,'B',5,'C',1,'K',1,'nu',0.5,'M',0.5,'Q',2);
        else
            params = struct('A',0,'B',3,'C',1,'K',1,'nu',0.5,'M',0.5,'Q',2);
        end
        params = struct('A',0,'B',5,'C',1,'K',1,'nu',2,'M',x1,'Q',2);
        yvar = general_logistic(var,params,~ispos);        
    elseif strcmp(options.type,'Linear')
        if ispos
            yvar = (var-xmin)./(xmax-xmin);
        else
            yvar = (xmax-var)./(xmax-xmin);
        end
    else
        x2 = 1;
        yvar = zeros(size(var));   
        idx1 = var<=x1;
        idx2 = var>x1 & var<=x2;
        idx3 = ~(idx1 | idx2);
        if ispos                                         %function increasing     
            y0 = options.SL;
            yvar(idx1) = y0(1)*(var(idx1)-xmin)./(x1-xmin);
            yvar(idx2) = y0(1)+(var(idx2)-x1)./(x2-x1)*(y0(2)-y0(1));
            yvar(idx3) = y0(2)+(var(idx3)-x2)./(xmax-x2).*(1-y0(2));
        else                                             %function decreasing
            y0 = options.BP;
            yvar(idx1) = y0(1)+(x1-var(idx1))./(x1-xmin).*(1-y0(1));
            yvar(idx2) = y0(2)+(x2-var(idx2))./(x2-x1)*(y0(1)-y0(2));
            yvar(idx3) = (xmax-var(idx3))./(xmax-x2).*y0(2);
        end
    end
end

%%
function [BPweightscores,SLweightscores] = weighted_scores(BioPhysData,SocialData,YrPopData)
    %set weights, apply to both datasets and plot scores and weighted
    %scores as radial plot and thumbnail plots

    %get the user to define weights for the variables
%     BPweights = setVariableWeights(BioPhysData);
%     SLweights = setVariableWeights(SocialData);
%     if isempty(BPweights) || isempty(SLweights), return; end %user cancelled
    BPweights = table(.2,.1,.1,.1,.2,.15,.15,'VariableNames',BioPhysData.Properties.VariableNames);
    SLweights = table(0.05,.15,.15,.1,.1,.1,.1,0.05,0.05,.1,0.05,'VariableNames',SocialData.Properties.VariableNames);
    %apply weights to biophysical scores-----------------------------------
    [BPscores,SLscores,options] = scaled_data(BioPhysData,SocialData,YrPopData,false);
    varnames = BPscores.Properties.VariableNames;
    nbp = length(varnames);
    BPweightscores = BPscores;
    for i=1:nbp
        BPweightscores.(varnames{i}) = BPscores.(varnames{i})*BPweights.(varnames{i});
    end
    BPtotalscore = sum(BPscores{:,:},2,'omitnan')/nbp;
    BPweightedscore = sum(BPweightscores{:,:},2,'omitnan');

    %apply weights to biophysical scores-----------------------------------
    varnames = SLscores.Properties.VariableNames;
    nsl = length(varnames);
    SLweightscores = SLscores;
    for i=1:nsl
        SLweightscores.(varnames{i}) = SLscores.(varnames{i})*SLweights.(varnames{i});
    end
    SLtotalscore = sum(SLscores{:,:},2,'omitnan')/nsl;
    SLweightedscore = sum(SLweightscores{:,:},2,'omitnan');

    %plot of scores--------------------------------------------------------
    hfs = thumbnailplot(BioPhysData,SocialData,BPscores,SLscores,YrPopData,options,true);
    plotitle = findobj(hfs.Children,'Tag','Thumb');
    mnBPscore = mean(BPtotalscore,'omitnan');
    mnSLscore = mean(SLtotalscore,'omitnan');
    plotitle.String = sprintf('Metric scores, means: R_B_P=%.2f; R_S=%.2f\nCircles represent population size',...
                               mnBPscore,mnSLscore);
    

    hfw = thumbnailplot(BioPhysData,SocialData,BPweightscores,SLweightscores,YrPopData,options,false);
    plotitle = findobj(hfw.Children,'Tag','Thumb');  
    mnBPwtscore = mean(BPweightedscore,'omitnan');
    mnSLwtscore = mean(SLweightedscore,'omitnan');
    plotitle.String = sprintf('Weighted metric scores, means: R_B_P=%.2f; R_S=%.2f\nCircles represent population size',...
                               mnBPwtscore,mnSLwtscore);

    weightedscores_plot(BPweightedscore,SLweightedscore,YrPopData,1000);    
    hfs = weightedscores_plot(BPtotalscore,SLtotalscore,YrPopData,1000);
    % target_marker(mnBPscore,mnSLscore,'SizeData',100,'Alpha',0.2,'MarkerEdgeColor','r','MarkerFaceColor','r');
    plotaxes = findobj(hfs.Children,'Tag','WSplot');
    plotaxes.XLabel.String = 'BioPhysical score';
    plotaxes.YLabel.String = 'Social score';
    plotaxes.Title.String = 'Metric scores by country';
end

%%
function hf = weightedscores_plot(BPvalues,SLvalues,YrPopData,mfactor)
    %create plot to reproduce variant of Figure 2 in paper using weighted
    %scores
    blue = mcolor('dark blue');   
    grey = mcolor('light grey');
    %plot figure
    hf = figure('Tag','PlotFig');
    ax = axes(hf);
    markersz = mfactor*YrPopData.Population/max(YrPopData.Population);
    scatter(ax,BPvalues,SLvalues,markersz,'filled',...
        'MarkerEdgeColor',blue,'MarkerFaceAlpha',0.1)
    xlabel('Weighted BioPhysical score')
    ylabel('Weighted social score')
    title('Weighted metric scores by country')
    subtitle('Circles scaled by population size')    
    hold on
    N = 20;
    dist = 0.2:0.2:1;
    for j=1:length(dist)
        alp = 0:pi/2/N:pi/2;
        x = dist(j)*cos(alp);
        y = dist(j)*sin(alp);
        plot(ax,x,y,'Color',grey,'Tag','plotobj')
    end
    hold off
    axis equal
    ax.XLim = [0,1];
    ax.YLim = [0,1];    
    ax.Tag = 'WSplot';
    
    %add index point
    % mnBPscore = mean(BPvalues,'omitnan');
    % mnSLscore = mean(SLvalues,'omitnan');
    % target_marker(mnBPscore,mnSLscore,'SizeData',100,'Alpha',0.2,...
    %                           'MarkerEdgeColor','r','MarkerFaceColor','r');

    pwm = pop_weighted_mean(BPvalues,SLvalues,YrPopData.Population);
    target_marker(pwm.x,pwm.z,'SizeData',150,'Alpha',0.2,'Tag','plotobj',...
                              'MarkerEdgeColor','r','MarkerFaceColor','r');
    text(ax,0.8,0.9,sprintf('R-I = %.2f\nB-I = %.2f\nN = %.0f',pwm.RI,pwm.BI,pwm.Nrec));
end

%%
function pwm = pop_weighted_mean(BPvalues,SLvalues,popdata)
    %estimate the population weighted mean value
    idincl = ~(isnan(BPvalues) | isnan(SLvalues));
    totpop = sum(popdata(idincl));
    x = sum((BPvalues(idincl).*popdata(idincl))./totpop);
    z = sum((SLvalues(idincl).*popdata(idincl))./totpop);
    RI = sqrt(x.^2+z.^2);
    BI = (atan(z./x)*180/pi-45)/45;  %report as (angle from 45 degree line)/45 deg
    Nrec = sum(idincl);
    pwm = table(x,z,RI,BI,Nrec);
end

%%
function hf = thumbnailplot(BioPhysData,SocialData,BPscaled,SLscaled,YrPopData,y0,islim)
    %Thumbnail plots of all variables from the biophysical and social
    %datasets
    BPvartxt = BioPhysData.Properties.VariableDescriptions;
    nBPvar = width(BioPhysData);
    markersz = 100*YrPopData.Population/max(YrPopData.Population);

    %plot the biophysical datasets
    blue = mcolor('dark blue');
    hf = figure('Tag','PlotFig');
    ax = axes(hf);
    subplot(6,3,1,ax);
    s = gobjects(1,nBPvar);
    for i=1:nBPvar
        s(i) = subplot(6,3,i);
        xvar = BioPhysData.(i);
        yvar = BPscaled.(i);
        scatter(s(i),xvar,yvar,markersz,'filled',...
                            'MarkerEdgeColor',blue,'MarkerFaceAlpha',0.1)
        title(BPvartxt{i})
        if any(i==[1,4,7,10,13,16])
            ylabel('Score')
        end
        if any(i==16:18)
            xlabel('Variable')
        end     
        if islim
            s(i).YLim = [0,1];
        end
    end

    %plot the social datasets
    SLvartxt = SocialData.Properties.VariableDescriptions;
    nSLvar = width(SocialData);
    for i=1:nSLvar
        j = nBPvar+i;
        s(j) = subplot(6,3,j);
        xvar = SocialData.(i);
        yvar = SLscaled.(i);
        scatter(s(j),xvar,yvar,markersz,'filled',...
                            'MarkerEdgeColor',blue,'MarkerFaceAlpha',0.1)
        title(SLvartxt{i})
        if any(j==[1,4,7,10,13,16])
            ylabel('Score')
        end
        if any(j==16:18)
            xlabel('Metric values')
        end    
        if islim
            s(j).YLim = [0,1];
        end
    end
    if ~isempty(y0.BP) && ~isnan(y0.BP(1))
        titletxt = sprintf('Metric scores (y_B_P= %.1f/%.1f; y_S_L= %.1f/%.1f)\nCircles represent population size',...
                                     y0.BP(1),y0.BP(2),y0.SL(1), y0.SL(2));
    elseif ~isempty(y0.BP) && isnan(y0.BP(1))
        titletxt = sprintf('Metric scores using Linear slope over data range\nCircles represent population size');   
    else
        titletxt = sprintf('Metric scores using Logistic curve\nCircles represent population size');       
    end
    sgtitle(titletxt,'Tag','Thumb')       
end

%%
function varweights = setVariableWeights(var)
    %prompt user to define weights for all variables in data set
    %NB must sum to zero
    varlist = var.Properties.VariableDescriptions;
    varname = var.Properties.VariableNames;
    nvar = length(varlist);
    defaults = cellstr(num2str(ones(nvar,1)*10/nvar));
    promptxt = 'Set weights for each variable (must sum to 10)';
    varlist{1} = sprintf('%s\n%s',promptxt,varlist{1});
    ok = 0;
    while ok==0
        answers = inputdlg(varlist,'Weights',1,defaults);
        if isempty(answers), varweights = []; return; end
        sumvar = 0;
        for i=1:nvar
            varweights.(varname{i}) = str2double(answers{i})/10;
            sumvar = sumvar+varweights.(varname{i});
        end
        if ismembertol(sumvar,1,0.01), ok=1; end
    end
end

%%
function country_selection(BioPhysData,SocialData,YrPopData)
    %prompt user to select country and plot scores and weighted scores on
    %radar plots
    %get the user to define weights for the variables
%     BPweights = setVariableWeights(BioPhysData);
%     SLweights = setVariableWeights(SocialData);
%     if isempty(BPweights) || isempty(SLweights), return; end %user cancelled

    BPweights = table(.2,.1,.1,.1,.2,.15,.15,'VariableNames',BioPhysData.Properties.VariableNames);
    SLweights = table(0.05,.15,.15,.1,.1,.1,.1,0.05,0.05,.1,0.05,'VariableNames',SocialData.Properties.VariableNames);
    %apply weights to biophysical scores-----------------------------------
    [BPscores,SLscores,~] = scaled_data(BioPhysData,SocialData,YrPopData,false);
    varnames = BPscores.Properties.VariableNames;
    nbp = length(varnames);
    BPweightscores = BPscores;
    for i=1:nbp
        BPweightscores.(varnames{i}) = BPscores.(varnames{i})*BPweights.(varnames{i});
    end

    %apply weights to biophysical scores-----------------------------------
    varnames = SLscores.Properties.VariableNames;
    nsl = length(varnames);
    SLweightscores = SLscores;
    for i=1:nsl
        SLweightscores.(varnames{i}) = SLscores.(varnames{i})*SLweights.(varnames{i});
    end

    Country = BioPhysData.Properties.RowNames;

    ok = 1;
    while ok>0
        selection = listdlg("PromptString",'Select a country:',...
                               'SelectionMode','single','ListString',Country);    
        if isempty(selection), ok=0; continue; end
    
        var.BPscores = BPscores(selection,:);
        var.BPweights = BPweightscores(selection,:);
        var.SLscores = SLscores(selection,:);
        var.SLweights = SLweightscores(selection,:);      
        country_score_plot(var,Country{selection});        
    end
end

%%
function country_score_plot(var,country)
    %plot figure
    hf = figure('Tag','PlotFig');
    ax = axes(hf); 
    blue = mcolor('dark blue'); green = mcolor('green');
 
    scoredata = var.BPscores{1,:};    
    scoredata(scoredata<0.01 | isnan(scoredata)) = 0.01;
    axlabels = var.BPweights.Properties.VariableDescriptions;
    nrec = length(scoredata);  
    dataset = var.BPweights{1,:}*nrec;
    dataset(dataset<0.01 | isnan(dataset)) = 0.01;
    axlim = max(1,max(dataset));    
    dataset = [dataset;scoredata];

    s1 = subplot(1,2,1,ax);
    spider_plot(dataset,'AxesHandle',s1,...
                    'AxesLabels',axlabels,...
                    'AxesLimits',repmat([0;axlim],1,nrec),...
                    'AxesOffset',0,...
                    'AxesPrecision', 2,...
                    'AxesDisplay', 'one',... %limit axes labels to one axis
                    'AxesLabelsEdge','none',...%border for axis labels off
                    'Color',[blue;green],...
                    'Marker','none',...
                    'FillOption','on',...
                    'FillTransparency',0.5);
    hl1 = legend('Weighted','Unweighted','Location','southwest');
    hl1.Position(1:2) = [0.4,0.2];
    title(sprintf('%s - BioPhysical',country))
    clear dataset scoredata axlabels axlim nrec

    scoredata = var.SLscores{1,:};
    scoredata(scoredata<0.01 | isnan(scoredata)) = 0.01;
    axlabels = var.SLweights.Properties.VariableDescriptions;
    nrec = length(scoredata);  
    dataset = var.SLweights{1,:}*nrec;   
    dataset(dataset<0.01 | isnan(dataset)) = 0.01;
    axlim = max(1,max(dataset));
    dataset = [dataset;scoredata];

    s2 = subplot(1,2,2);
    spider_plot(dataset,'AxesHandle',s2,...
                    'AxesLabels',axlabels,...
                    'AxesLimits',repmat([0;axlim],1,nrec),...
                    'AxesOffset',0,...
                    'AxesPrecision', 2,...
                    'AxesDisplay', 'one',... %limit axes labels to one axis
                    'AxesLabelsEdge','none',...%border for axis labels off
                    'Color',[blue;green],...
                    'Marker','none',...
                    'FillOption','on',...
                    'FillTransparency',0.5);
    %legend('Weighted','Unweighted','Location','southwest')
    title(sprintf('%s - Social',country))
    sgtitle(sprintf('Weighted and Unweighted scores\n(Weighted scores are multiplied by number of variables)'))    
end

%--------------------------------------------------------------------------
% Radial-Index timeseries plots
%% ------------------------------------------------------------------------
function index_ts(BioPhysTS,SocialTS,PopData)
    %compute and plot the time varying index
    yearvar = fieldnames(BioPhysTS);
    nyear = length(yearvar);
    temp = split(yearvar,'r');
    yeardata = str2double(temp(:,2));

    %get the user to define weights for the variables
    BPdata1 = subsampleTable(BioPhysTS,1,true);
    SLdata1 = subsampleTable(SocialTS,1,false);

%     BPweights = setVariableWeights(BPdata1);
%     SLweights = setVariableWeights(SLdata1);
%     if isempty(BPweights) || isempty(SLweights), return; end %user cancelled
    BPweights = table(.2,.1,.1,.1,.2,.15,.15,'VariableNames',BPdata1.Properties.VariableNames);
    SLweights = table(0.05,.15,.15,.1,.1,.1,.1,0.05,0.05,.1,0.05,'VariableNames',SLdata1.Properties.VariableNames);
    %select score scaling to use
    options = selectScoreScale();

    %BPpws = zeros(nyear,1); SLpws = BPpws; RI = BPpws; phi = BPpws; Nrec = BPpws;   
    pwm = [];
    for i=1:nyear
        %extract data for a year and format to match O'Neill dataset
        anBPdata = subsampleTable(BioPhysTS,i,true);
        anSLdata = subsampleTable(SocialTS,i,false);

        anPopData = PopData.(yearvar{i});
        [BioPhysData,SocialData,anPopData] = sortPop(anBPdata,anSLdata,anPopData);
        BioPhysData.Unused = [];

        %apply weights to biophysical scores-----------------------------------
        [BPscores,SLscores,~] = scaled_data(BioPhysData,SocialData,anPopData,false,options);
        varnames = BPscores.Properties.VariableNames;
        nbp = length(varnames);
        BPweightscores = BPscores;
        for ii=1:nbp
            BPweightscores.(varnames{ii}) = BPscores.(varnames{ii})*BPweights.(varnames{ii});
        end
        BPweightedscore = sum(BPweightscores{:,:},2,'omitnan');

        %apply weights to biophysical scores-----------------------------------
        varnames = SLscores.Properties.VariableNames;
        nsl = length(varnames);
        SLweightscores = SLscores;
        for jj=1:nsl
            SLweightscores.(varnames{jj}) = SLscores.(varnames{jj})*SLweights.(varnames{jj});
        end
        SLweightedscore = sum(SLweightscores{:,:},2,'omitnan');

        pwm = [pwm;pop_weighted_mean(BPweightedscore,SLweightedscore,anPopData.Population)]; %#ok<AGROW> 
        if i==1
            first.yrBP = BPweightedscore;
            first.yrSL = SLweightedscore;
            first.yrPop = anPopData;
        elseif i==nyear
            last.yrBP = BPweightedscore;
            last.yrSL = SLweightedscore;
            last.yrPop = anPopData;
        end
        clear anBPdata anSLdata anPopData BPscores BPweightscores SLscores SLweightscores
    end
    radial_index_ts_plot(yeardata,pwm)
    weightedscores_ts_plot(yeardata,first,last,pwm)
end

%%
function radial_index_ts_plot(yr,pwm)
    %plot a time series of RI and phi
    hf = figure('Tag','PlotFig');
    ax = axes(hf); 
    blue = mcolor('dark blue'); orange = mcolor('orange');
    yyaxis left
    plot(ax,yr,pwm.RI,'Color',blue);
    xlabel('Year')
    ylabel('Radial-Index')
    yyaxis right
    plot(ax,yr,pwm.BI,'Color',orange);
    ylabel('Bias-Index')
    legend({'Radial-Index','Bias-Index'},'Location','Best')
    title(sprintf('Global indices for years %d-%d',yr(1),yr(end)))
end

%%
function weightedscores_ts_plot(yr,first,last,pwm)
    %Weighted scores plot for first and last year with trace of 
    %radial-index over time
    mfactor = 1000;  %factor to scale population markers
    orange = mcolor('orange');
    hf = weightedscores_plot(last.yrBP,last.yrSL,last.yrPop,mfactor);
    ax = findobj(hf.Children,'Type','axes');
    hold on
    markersz = mfactor*last.yrPop.Population/max(last.yrPop.Population);
    scatter(ax,first.yrBP,first.yrSL,markersz,'filled',...
                            'MarkerEdgeColor',orange,'MarkerFaceAlpha',0.05);
    plot(ax,pwm.x,pwm.z,'r','LineWidth',1);
    plot(ax,pwm.x(1),pwm.z(1),'xr','MarkerSize',10,'Tag','plotobj')
    hold off
    hl = findobj(ax.Children,'Tag','plotobj');
    for i=1:length(hl)
        hl(i).Annotation.LegendInformation.IconDisplayStyle = 'off';  
    end
    txt1 = sprintf('Data for %d',yr(1));
    txt2 = sprintf('Data for %d',yr(end));
    legend({txt1,txt2,'R-I (t) trace'},'Location','north')
    title(sprintf('Initial and Final weighted scores and Radial-Index for years %d-%d',...
                             yr(1),yr(end)))
end
%--------------------------------------------------------------------------
% Utlities
%% ------------------------------------------------------------------------
function options = selectScoreScale()
%select the type of scaling and define the limits if Linear segments
    options.type = questdlg('Select scaling function','Scaling','Logistic','Linear segments','Linear','Logistic');

    if strcmp(options.type,'Linear segments')
        answers = inputdlg({'BioPhysical threshold at x=0.8 and 1','Social threshold at x=0.8 and 1'},...
                                        'Thresholds',1,{'0.8 0.2','0.5 0.8'});
        if isempty(answers), return; end
        options.BP = str2num(answers{1});  %#ok<ST2NM> %score value for variable-threshold
        options.SL = str2num(answers{2});  %#ok<ST2NM> 
    elseif strcmp(options.type,'Linear')
        options.BP = NaN; options.SL = NaN;
    else
        options.BP = []; options.SL = [];
    end
end
%%
function [BPexPop,SLexPop,PopData] = sortPop(BioPhysData,SocialData,YrPopData)
    %combine the O'Neil et al, 2018, data and population data from 
    %https://ourworldindata.org/grapher/population
        %Add population data to the data sets
    BPdata = addvars(BioPhysData,BioPhysData.Properties.RowNames,...
                                          'NewVariableName',{'Country'});
    BPexPop = innerjoin(BPdata,YrPopData,...
                    'LeftKeys',8,'RightKeys',1,'RightVariables',4);
    BPexPop.Properties.RowNames = BPexPop.Country;
    BPexPop.Country = [];

    SLdata = addvars(SocialData,SocialData.Properties.RowNames,...
                                          'NewVariableName',{'Country'});
    SLexPop = innerjoin(SLdata,YrPopData,...
                    'LeftKeys',12,'RightKeys',1,'RightVariables',4);
    SLexPop.Properties.RowNames = SLexPop.Country;
    SLexPop.Country = [];
    SLexPop.Population = [];

    PopData = BPexPop;
    PopData(:,1:7) = [];  %Remove all but population data
    BPexPop.Population = [];
end
%%
function var = subsampleTable(tsvar,idx,isBP)
    %subsample the Fanning timeseries data to match O'Neill format
    yearvar = fieldnames(tsvar);
    blank = NaN(height(tsvar.(yearvar{1})),1);
    var = tsvar.(yearvar{idx});
    var.Properties.RowNames = var.Country;
    var(:,1:3) = [];                    %remove location and year
    if isBP
        %var = addvars(var,blank,'After',3,'NewVariableName','Unused');
        var = addvars(var,blank,'NewVariableName','Unused');
    end
end
  