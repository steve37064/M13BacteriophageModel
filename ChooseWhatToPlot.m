
close all 
GeneratedModels
ListOfAnalytesToPlot = ["Phage","P1","P2","P3","P4","P5","P6","P7","P8"...
                        ,"P9","P10","P11"];
                    
PlotList = ["RemakeOrginal.m","Swap_3_5.m","Swap_4_9.m"]; 

ShortList = ismember(GeneratedModels,PlotList);
DataStoreShortList = DataStore(ShortList);
ModelInfoStoreShortList = ModelInfoStore(ShortList);
GeneratedModelsShortList = GeneratedModels(ShortList);
%==========================================================================
%Plotting Models: 
%==========================================================================
timepointsMinutes = timepoints./60;
ListOfAnalytesToPlot = ["Phage","P1","P2","P3","P4","P5","P6","P7","P8"...
                        ,"P9","P10","P11"];

ColAnalyte = 0;
for Analyte = ListOfAnalytesToPlot
    ColAnalyte = ColAnalyte + 1; 
    OrginalAnalyte = GetSimulatedData(Analyte,OrginalModelObservables,OrginalModelInfo);
    Samples(1,ColAnalyte) = OrginalAnalyte(end);
    figure 
    plot(timepointsMinutes,OrginalAnalyte,"--","DisplayName","Orginal","LineWidth",4)
    hold on 
    for i = 1:length(GeneratedModelsShortList)
        NextAnalytesToPlot = GetSimulatedData(Analyte,DataStoreShortList{i},ModelInfoStoreShortList{i});
        Samples(i+1,ColAnalyte) = NextAnalytesToPlot(end);
        ModelToRun = split(GeneratedModelsShortList{i},".");
        ModelToRun = ModelToRun{1};
        plot(timepointsMinutes,NextAnalytesToPlot,"DisplayName",ModelToRun,"LineWidth",4)
    end 
    legend('Interpreter', 'none',"FontSize",20)
    grid
    xlabel("Time [Minutes]","FontSize",20)
    ylabel(Analyte,"FontSize",20)
    %SaveName = "Graphs/04_18_2020___1___RemoveSwaps/"+Analyte+".png";
    %saveas(gca,SaveName)
end 


function [Data2Plot] = GetSimulatedData(Analyte,Observables,ModelInfo)
    ObservLoc = ismember(ModelInfo("observable_labels"),Analyte); 
    if sum(ObservLoc) ~= 1 
        disp("Something is wrong!") 
        disp(Analyte+" not found!!!!")
    else 
        Data2Plot = Observables(:,ObservLoc); 
    end 
    
end 