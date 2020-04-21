MainOutputFolder = "../../../Box/PhD/Spring_2020_Classes/Cell_and_Systems/CellSystemsProject/Graphs/";
CurrentTime = clock;
DateDir = CurrentTime(1)+"_"+CurrentTime(2)+"_"+CurrentTime(3)+"___";
TimeDir = CurrentTime(4)+"_"+CurrentTime(5)+"_"+floor(CurrentTime(6))+"___"; 
PersonalLabel = "ClusterByPhage_Factor10"; 
%PersonalLabel = "Exclude_2_10_1_11_UpdateLegend_Log10";

OutputDirectory = MainOutputFolder+DateDir+TimeDir+PersonalLabel;

Log10Plot = false; 
PlotOrginalOnEveryPlot = false;
Status = mkdir(OutputDirectory);

%Clustering 1 
PlotListList = {
    ["Orginal","RemakeOrginal.m","Swap_10_11.m","Swap_1_10.m","Swap_1_11.m","Swap_1_3.m","Swap_1_6.m","Swap_1_8.m","Swap_2_10.m","Swap_3_10.m","Swap_3_11.m","Swap_3_4.m","Swap_3_6.m","Swap_3_8.m","Swap_4_8.m","Swap_6_10.m","Swap_6_11.m","Swap_6_8.m","Swap_8_9.m"],...
    ["Swap_1_2.m","Swap_2_11.m","Swap_2_3.m","Swap_2_5.m","Swap_2_6.m","Swap_4_5.m","Swap_5_10.m","Swap_5_6.m","Swap_5_9.m"],...
    ["Swap_1_9.m","Swap_3_9.m","Swap_4_9.m","Swap_6_9.m"],...
    ["Swap_1_4.m","Swap_2_4.m","Swap_2_8.m","Swap_4_10.m","Swap_4_11.m","Swap_4_6.m","Swap_5_8.m","Swap_8_10.m"],...
    ["Swap_8_11.m","Swap_9_11.m"],...
    ["Swap_1_5.m","Swap_2_9.m","Swap_3_5.m","Swap_9_10.m"],...
    ["Swap_5_11.m"]...
};

%Clustering 2 
% excluding: [2, 10, 1, 11]
%PlotListList = {
%    ["Orginal","RemakeOrginal.m","Swap_3_4.m","Swap_3_6.m","Swap_3_8.m","Swap_4_8.m","Swap_6_8.m","Swap_8_9.m"], 
%    ["Swap_4_5.m","Swap_4_6.m","Swap_5_8.m","Swap_5_9.m"],
%    ["Swap_3_5.m","Swap_3_9.m","Swap_4_9.m","Swap_6_9.m"],
%    ["Swap_5_6.m"]
%    };


PlotListList = {
 ["Swap_1_2.m","Swap_2_11.m","Swap_9_10.m","Swap_1_5.m","Swap_8_10.m","Swap_5_10.m"], ... %No Phage 
 ["Swap_8_11.m","Swap_1_8.m","Swap_3_8.m","Swap_5_11.m"], ... %0<Phage<10 
 ["Swap_6_8.m","Swap_3_5.m","Swap_3_9.m","Swap_4_9.m","Swap_5_6.m","Swap_4_8.m","Swap_1_9.m","Swap_2_6.m","Swap_4_5.m","Swap_2_8.m","Swap_2_3.m","Swap_2_5.m"], ... %10<Phage<100 
 ["Swap_6_9.m","Swap_5_8.m","Swap_5_9.m","Swap_9_11.m","Swap_2_9.m","Swap_3_4.m","Swap_10_11.m","Swap_3_10.m","Swap_1_6.m","Swap_6_11.m","Swap_1_11.m","Swap_8_9.m","Swap_3_6.m","RemakeOrginal.m","Swap_1_3.m","Swap_3_11.m","Swap_4_11.m","Swap_4_10.m","Swap_1_4.m","Swap_4_6.m","Swap_6_10.m","Orginal","Swap_2_10.m","Swap_1_10.m","Swap_2_4.m"]... %100<Phage<1000 
}
TitleNames = ["No Phage Produced","0<Page<10","10<Page<100","100<Page<1000"]; 



         
for ClusterNum = 1:length(PlotListList)
    close all 
    PlotList = PlotListList{ClusterNum};

    ShortList = ismember(GeneratedModels,PlotList);
    DataStoreShortList = DataStore(ShortList);
    ModelInfoStoreShortList = ModelInfoStore(ShortList);
    GeneratedModelsShortList = GeneratedModels(ShortList);
    %==========================================================================
    %Plotting Models: 
    %==========================================================================
    timepointsMinutes = timepoints./60;
    ListOfAnalytesToPlot = ["Phage","P1","P2","P3","P4","P5","P6","P7","P8"...
                            ,"P9","P10","P11","P2P10","As","PI","PE","PF"];

    ColAnalyte = 0;
    for Analyte = ListOfAnalytesToPlot
        ColAnalyte = ColAnalyte + 1; 
        OrginalAnalyte = GetSimulatedData(Analyte,OrginalModelObservables,OrginalModelInfo);
        Samples(1,ColAnalyte) = OrginalAnalyte(end);
        figure 
        if PlotOrginalOnEveryPlot
            if Log10Plot
                plot(timepointsMinutes,log10(abs(OrginalAnalyte)),"--","DisplayName","Orginal","LineWidth",4)
            else 
                plot(timepointsMinutes,OrginalAnalyte,"--","DisplayName","Orginal","LineWidth",4)
            end 
        end 
        hold on 
        for i = 1:length(GeneratedModelsShortList)
            NextAnalytesToPlot = GetSimulatedData(Analyte,DataStoreShortList{i},ModelInfoStoreShortList{i});
            Samples(i+1,ColAnalyte) = NextAnalytesToPlot(end);
            ModelToRun = split(GeneratedModelsShortList{i},".");
            ModelToRun = ModelToRun{1};
            if Log10Plot
                plot(timepointsMinutes,log10(abs(NextAnalytesToPlot)),"DisplayName",ModelToRun,"LineWidth",4)
            else 
                plot(timepointsMinutes,NextAnalytesToPlot,"DisplayName",ModelToRun,"LineWidth",4)
            end 
        end 
        legend('Location','NorthEastOutside','Interpreter', 'none',"FontSize",20)
        grid
        xlabel("Time [Minutes]","FontSize",20)
        if Log10Plot
            ylabel("Log_{10}("+Analyte+")","FontSize",20)
            ylim([0 inf])
        else 
            ylabel(Analyte,"FontSize",20)
            ylim([0 inf])
        end 
        if exist("TitleNames")
            title(TitleNames(ClusterNum),"FontSize",20)
        end 
        SAVDIR = OutputDirectory + "/";
        Specific = "Cluster_"+num2str(ClusterNum-1) + "___" + Analyte + ".png";
        SaveName = SAVDIR+Specific;
        saveas(gca,SaveName)
    end 
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