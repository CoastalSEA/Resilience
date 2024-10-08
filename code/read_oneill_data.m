%script to read O'Neill et al data from spreadsheet and save to mat file 
%as a table
%source files in ../Work/Research/Resilience
    userprompt = 'Select Excel file>';
    [fname, path]=uigetfile('*.xlsx',userprompt,'MultiSelect','off');
    if fname==0, return; end
    filename = [path,fname];

    varnames = {'CO2','Phosphorus','Nitrogen','Water','eHANPP','EcoFoot','MatFoot'};
    opts = spreadsheetImportOptions('NumVariables',7,...
                                'Sheet','BioPhysical',...
                                'VariableNames',varnames,...
                                'VariableDescriptionsRange','B1:H1',...
                                'VariableTypes',repmat({'double'},1,7),...
                                'RowNamesRange','A2:A152',...
                                'DataRange', 'B2:H152'); 

    BioPhysData = readtable(filename,opts);

    varnames = {'LifeSat','HealthyLife','Nutrition','Sanitation','Income',...
                'Energy','Educatiom','SocialSupport','Democracy','Equality','Employment'};
    opts = spreadsheetImportOptions('NumVariables',11,...
                                'Sheet','Social',...
                                'VariableNames',varnames,...
                                'VariableDescriptionsRange','B1:L1',...
                                'VariableTypes',repmat({'double'},1,11),...
                                'RowNamesRange','A2:A152',...
                                'DataRange', 'B2:L152'); 

    SocialData = readtable(filename,opts);

    save([path,'ONeillData.mat'],'BioPhysData','SocialData','-mat')

    clearvars