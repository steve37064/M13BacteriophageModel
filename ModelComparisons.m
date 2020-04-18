
set(0,'DefaultFigureWindowStyle','docked'); 
clc;
clear all; 
close all; 
timepoints = linspace(0,1*60*60,1000)';


%==========================================================================
%Orginal Model: Current Gold Standard
%==========================================================================
OrginalModel='Model';
OrginalModelInfo = getValues(OrginalModel);
OrginalModelObservables = SimulateModel(OrginalModel,timepoints);


%==========================================================================
%Generated Models: 
%==========================================================================
GeneratedModelFolder = "./GeneratedModels/MatlabVersions/";
addpath(genpath(GeneratedModelFolder))
GeneratedModels = dir(GeneratedModelFolder);
GeneratedModels = {GeneratedModels.name};
GeneratedModels = {GeneratedModels{3:end}};

DataStore = {};
ModelInfoStore = {};

for ModelTest = 1:length(GeneratedModels)
    ModelToRun = split(GeneratedModels{ModelTest},".");
    ModelToRun = ModelToRun{1};
    DataStore{ModelTest} = SimulateModel(ModelToRun,timepoints);
    ModelInfoStore{ModelTest} =  getValues(GeneratedModelFolder+ModelToRun);
end  

%==========================================================================
%Plotting Models: 
%==========================================================================
timepointsMinutes = timepoints./60;
ListOfAnalytesToPlot = ["Phage","P1","P2","P3","P4","P5","P6","P7","P8"...
                        ,"P9","P10","P11","P2P10","As","PI","PE","PF"];
                    
Samples = zeros(length(GeneratedModels)+1,length(ListOfAnalytesToPlot));
ColAnalyte = 0;
for Analyte = ListOfAnalytesToPlot
    ColAnalyte = ColAnalyte + 1; 
    OrginalAnalyte = GetSimulatedData(Analyte,OrginalModelObservables,OrginalModelInfo);
    Samples(1,ColAnalyte) = OrginalAnalyte(end);
    figure 
    plot(timepointsMinutes,OrginalAnalyte,"--","DisplayName","Orginal","LineWidth",4)
    hold on 
    for i = 1:length(GeneratedModels)
        NextAnalytesToPlot = GetSimulatedData(Analyte,DataStore{i},ModelInfoStore{i});
        Samples(i+1,ColAnalyte) = NextAnalytesToPlot(end);
        ModelToRun = split(GeneratedModels{i},".");
        ModelToRun = ModelToRun{1};
        plot(timepointsMinutes,NextAnalytesToPlot,"DisplayName",ModelToRun,"LineWidth",4)
    end 
    %legend('Interpreter', 'none',"FontSize",20)
    grid
    xlabel("Time [Minutes]","FontSize",20)
    ylabel(Analyte,"FontSize",20)
    %legend
    SaveName = "Graphs/04_18_2020___1___RemoveSwaps/"+Analyte+".png";
    saveas(gca,SaveName)
end 

rowNames = {'Orginal' GeneratedModels{:}};
Table2Store = array2table(Samples,'VariableNames',ListOfAnalytesToPlot,"RowNames",rowNames);
writetable(Table2Store,"FinalValues.csv",'WriteRowNames',true)




function [Data2Plot] = GetSimulatedData(Analyte,Observables,ModelInfo)
    ObservLoc = ismember(ModelInfo("observable_labels"),Analyte); 
    if sum(ObservLoc) ~= 1 
        disp("Something is wrong!") 
        disp(Analyte+" not found!!!!")
    else 
        Data2Plot = Observables(:,ObservLoc); 
    end 
    
end 

%Runnig the model where we can choose which model to simulate
function [observables_out] = SimulateModel(RunModel,timepoints)
    Model2Run = str2func(RunModel);
    [~, ~, ~, observables_out] = Model2Run(timepoints,[],[],1);
end 

%Function to get the required information to run the model in Matlab 
function [ModelInfo] = getValues(model)
    % Get parameter names from param_labels string definition
    [~, result]= system(sprintf('grep "param_labels = { " %s.m', model));
    eval(result);
    % Get default parameter values from model 
    [~, result]= system(sprintf('grep "parameters = \\[ " %s.m', model));
    eval(result);
    [~, result]= system(sprintf('grep "observable_labels = { " %s.m', model));
    eval(result);
    ModelInfo = containers.Map(["param_labels", "parameters", "observable_labels"],...
                                {param_labels, parameters,  observable_labels}  );
end 


