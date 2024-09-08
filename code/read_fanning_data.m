%script to read Fanning et al data from spreadsheet and save to mat file 
%as a struct of years with tables of country data for each year.
% source files in ../Work/Research/Resilience

    userprompt = 'Select Excel file>';
    [fname, path]=uigetfile('*.xlsx',userprompt,'MultiSelect','off');
    if fname==0, return; end
    filename = [path,fname];

    varnames = {'Country','ID','Year','CO2','Phosphorus','Nitrogen','LandSystem','EcoFoot','MatFoot'};
    opts = spreadsheetImportOptions('NumVariables',9,...
                                'Sheet','BioPhysical_Historical',...
                                'VariableNames',varnames,...
                                'VariableDescriptionsRange','A1:I1',...
                                'VariableTypes',[{'char','char'},repmat({'double'},1,7)],...
                                'DataRange', 'A2:I3553'); 
    AllBioPhysData = readtable(filename,opts);

    % mnyr = min(AllPopData.Year,[],'omitnan');
    % mxyr = max(AllBioPhysData.Year,[],'omitnan');
    mnyr = 1992; mxyr = 2015;
    
    range = mnyr:1:mxyr;
    for i=1:length(range)
        yri = range(i);
        idx = find(AllBioPhysData.Year==yri);
        BioPhysData.(['yr',num2str(yri)]) =  AllBioPhysData(idx,:);
    end

    varnames = {'Country','ID','Year','LifeSat','HealthyLife','Nutrition','Sanitation','Income',...
                'Energy','Educatiom','SocialSupport','Democracy','Equality','Employment'};
    opts = spreadsheetImportOptions('NumVariables',14,...
                                'Sheet','Social_Historical',...
                                'VariableNames',varnames,...
                                'VariableDescriptionsRange','A1:N1',...
                                'VariableTypes',[{'char','char'},repmat({'double'},1,12)],...
                                'DataRange', 'A2:N3553'); 
    AllSocialData = readtable(filename,opts);

    range = mnyr:1:mxyr;
    for i=1:length(range)
        yri = range(i);
        idx = find(AllSocialData.Year==yri);
        SocialData.(['yr',num2str(yri)]) =  AllSocialData(idx,:);
    end

    save([path,'FanningData.mat'],'BioPhysData','SocialData','-mat')

    clearvars