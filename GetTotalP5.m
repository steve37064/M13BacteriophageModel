PlotListList = {};


close all 
Exlcude_1_2_10_11 = {'RemakeOrginal.m','Swap_3_4.m','Swap_3_6.m','Swap_3_8.m','Swap_4_8.m','Swap_6_8.m','Swap_8_9.m','Swap_4_5.m','Swap_5_9.m','Swap_3_5.m','Swap_3_9.m','Swap_4_9.m','Swap_6_9.m','Swap_5_6.m','Swap_4_6.m','Swap_5_8.m'};
TitleNames = {};
for i = 2:length(Exlcude_1_2_10_11)
    Orginal = string(Exlcude_1_2_10_11{1});
    Next = string(Exlcude_1_2_10_11{i});
    PlotListList{i-1} = [Orginal(1),Next(1) ];
    TitleNames{i-1} = "Comparision To Orginal";
end

for ClusterNum = 1:length(PlotListList)
    PlotList = PlotListList{ClusterNum};

    ShortList = ismember(GeneratedModels,PlotList);
    DataStoreShortList = DataStore(ShortList);
    ModelInfoStoreShortList = ModelInfoStore(ShortList);
    GeneratedModelsShortList = GeneratedModels(ShortList);
    
    P5_Orginal = GetSimulatedData("P5",DataStoreShortList{1},ModelInfoStoreShortList{1});
    P5_Swap = GetSimulatedData("P5",DataStoreShortList{2},ModelInfoStoreShortList{2});
    
    P5DNA_Orginal = GetSimulatedData("P5DNA",DataStoreShortList{1},ModelInfoStoreShortList{1});
    P5DNA_Swap = GetSimulatedData("P5DNA",DataStoreShortList{2},ModelInfoStoreShortList{2});
    
    Assembly_Site_Orginal = GetSimulatedData("PE",DataStoreShortList{1},ModelInfoStoreShortList{1});
    Assembly_Site_Swap    = GetSimulatedData("PE",DataStoreShortList{2},ModelInfoStoreShortList{2});

    
    TotalP5_Orginal = P5_Orginal + 1600*( P5DNA_Orginal + Assembly_Site_Orginal);
    TotalP5_Swap = P5_Swap       + 1600*( P5DNA_Swap    + Assembly_Site_Swap );
    
    figure 
    chosenLinestyle = {'-.','-','--',':'};
    LinstyleToPlot = chosenLinestyle{rem(1,4)+1};
    ModelToRun = split(GeneratedModelsShortList{1},".");
    ModelToRun = ModelToRun{1};
    plot(timepointsMinutes,TotalP5_Orginal,LinstyleToPlot,"DisplayName",ModelToRun,"LineWidth",2)
    hold on 
    LinstyleToPlot = chosenLinestyle{rem(2,4)+1};
    ModelToRun = split(GeneratedModelsShortList{2},".");
    ModelToRun = ModelToRun{1};
    plot(timepointsMinutes,TotalP5_Swap,LinstyleToPlot,"DisplayName",ModelToRun,"LineWidth",2)
    legend()
    %legend('Location','NorthEastOutside','Interpreter', 'none',"FontSize",15)
    legend('Interpreter', 'none',"FontSize",15)
    grid
    xlabel("Time [Minutes]","FontSize",20)
    ylabel("Total P5 (Bound + Unbound)","FontSize",20)
    
    FileName = split(PlotList(2),".");

    saveas(gca,"TotalP5/" + FileName(1)+"__Total_P5.png")
    
end 
%GeneratedModelsShortList = GeneratedModels(ShortList)
%NextAnalytesToPlot = GetSimulatedData(Analyte,DataStoreShortList{i},ModelInfoStoreShortList{i});
%Samples(i+1,ColAnalyte) = NextAnalytesToPlot(end);
%ModelToRun = split(GeneratedModelsShortList{i},".");
%ModelToRun = ModelToRun{1};

function [Data2Plot] = GetSimulatedData(Analyte,Observables,ModelInfo)
    ObservLoc = ismember(ModelInfo("observable_labels"),Analyte); 
    if sum(ObservLoc) ~= 1 
        disp("Something is wrong!") 
        disp(Analyte+" not found!!!!")
    else 
        Data2Plot = Observables(:,ObservLoc); 
    end 
    
end 
