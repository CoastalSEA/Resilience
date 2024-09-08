%script to read population data from spreadsheet and save to mat file 
%as a struct of years with tables of country data for each year.
% source data is from 
%   Gapminder, 2022, Gapminder - Population v7 and other sources – with 
%   major processing by Our World in Data, Our World in Data, 
%   https://ourworldindata.org/grapher/population.
% source files in ../Work/Research/Resilience
    userprompt = 'Select Excel file>';
    [fname, path]=uigetfile('*.xlsx',userprompt,'MultiSelect','off');
    if fname==0, return; end
    filename = [path,fname];

    varnames = {'Country','Code','Year','Population'};
    opts = spreadsheetImportOptions('NumVariables',4,...
                                'Sheet','population',...
                                'VariableNames',varnames,...
                                'VariableDescriptionsRange','A1:D1',...
                                'VariableTypes',{'char','char','double','double'},...
                                'DataRange', 'A2:D58253'); 


    AllPopData = readtable(filename,opts);

    mnyr = 1900;
    mxyr = max(AllPopData.Year,[],'omitnan');
    range = mnyr:1:mxyr;
    for i=1:length(range)
        yri = range(i);
        idx = find(AllPopData.Year==yri);
        PopData.(['yr',num2str(yri)]) =  AllPopData(idx,:);
    end
    outname = sprintf('PopData%d-%d.mat',mnyr,mxyr);
    save([path,outname],'PopData','-mat')

    clearvars


