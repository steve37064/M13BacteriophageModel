MainOutputFolder = "../../../Box/PhD/Spring_2020_Classes/Cell_and_Systems/CellSystemsProject/Graphs/";
%MainOutputFolder = "TestingGraphs/";
CurrentTime = clock;
DateDir = CurrentTime(1)+"_"+CurrentTime(2)+"_"+CurrentTime(3)+"___";
TimeDir = CurrentTime(4)+"_"+CurrentTime(5)+"_"+floor(CurrentTime(6))+"___"; 
PersonalLabel = "CompareEachToOrginal_Exclude_1_2_10_11"; 
PersonalLabel = "ClusterComparison_Exclude_1_2_10_11";
PersonalLabel = "CompareToOrginalRemade";
%PersonalLabel = "DirecComparisionForReport";

OutputDirectory = MainOutputFolder+DateDir+TimeDir+PersonalLabel;

Log10Plot = false; 
PlotOrginalOnEveryPlot = false;
Status = mkdir(OutputDirectory);

%Clustering 1 
PlotListList = {
    % cluster number: 0
       ["Orginal","Remake_Original.m","Swap_10_11.m","Swap_1_10.m","Swap_1_11.m","Swap_1_3.m","Swap_1_6.m","Swap_2_10.m","Swap_3_10.m","Swap_3_11.m","Swap_3_4.m","Swap_3_6.m","Swap_6_10.m","Swap_6_11.m","Swap_6_8.m","Swap_8_9.m"],
% cluster number: 1
       ["Swap_1_2.m","Swap_2_11.m","Swap_2_3.m","Swap_2_6.m","Swap_5_10.m"],
% cluster number: 2
       ["Swap_1_9.m","Swap_3_9.m","Swap_4_9.m","Swap_6_9.m","Swap_9_11.m"],
% cluster number: 3
       ["Swap_2_9.m","Swap_3_5.m","Swap_5_11.m","Swap_9_10.m"],
% cluster number: 4
       ["Swap_1_4.m","Swap_1_5.m","Swap_2_4.m","Swap_4_10.m","Swap_4_11.m","Swap_4_5.m","Swap_4_6.m","Swap_5_8.m","Swap_5_9.m"],
% cluster number: 5
       ["Swap_2_8.m","Swap_8_10.m"],
% cluster number: 6
       ["Swap_2_5.m","Swap_5_6.m"],
% cluster number: 7
       ["Swap_1_8.m","Swap_3_8.m","Swap_4_8.m","Swap_8_11.m"],
};

if false  
PlotListList = { ["Remake_Original.m"] }; 
TitleNames = ["Comparison to Original"];
end 

%Clustering 2 
% excluding: [2, 10, 1, 11]
if false 
PlotListList = {
% cluster number: 0
        ["Orginal","Remake_Original.m","Swap_3_4.m","Swap_3_6.m","Swap_3_8.m","Swap_4_8.m","Swap_6_8.m","Swap_8_9.m"],
% cluster number: 1
        ["Swap_4_5.m","Swap_5_9.m"],
% cluster number: 2
        ["Swap_3_5.m","Swap_3_9.m","Swap_4_9.m","Swap_6_9.m"],
% cluster number: 3
        ["Swap_5_6.m"],
% cluster number: 4
        ["Swap_4_6.m","Swap_5_8.m"]
    };
end 

%PlotListList = {
% ["Swap_1_2.m","Swap_2_11.m","Swap_9_10.m","Swap_1_5.m","Swap_8_10.m","Swap_5_10.m"], ... %No Phage 
% ["Swap_8_11.m","Swap_1_8.m","Swap_3_8.m","Swap_5_11.m"], ... %0<Phage<10 
% ["Swap_6_8.m","Swap_3_5.m","Swap_3_9.m","Swap_4_9.m","Swap_5_6.m","Swap_4_8.m","Swap_1_9.m","Swap_2_6.m","Swap_4_5.m","Swap_2_8.m","Swap_2_3.m","Swap_2_5.m"], ... %10<Phage<100 
% ["Swap_6_9.m","Swap_5_8.m","Swap_5_9.m","Swap_9_11.m","Swap_2_9.m","Swap_3_4.m","Swap_10_11.m","Swap_3_10.m","Swap_1_6.m","Swap_6_11.m","Swap_1_11.m","Swap_8_9.m","Swap_3_6.m","RemakeOrginal.m","Swap_1_3.m","Swap_3_11.m","Swap_4_11.m","Swap_4_10.m","Swap_1_4.m","Swap_4_6.m","Swap_6_10.m","Orginal","Swap_2_10.m","Swap_1_10.m","Swap_2_4.m"]... %100<Phage<1000 
%};
%TitleNames = ["No Phage Produced","0<Page<10","10<Page<100","100<Page<1000"]; 

if false  
    PlotListList = {};
    TitleNames = {};
    Exlcude_1_2_10_11 = {'Remake_Original.m','Swap_3_4.m','Swap_3_6.m','Swap_3_8.m','Swap_4_8.m','Swap_6_8.m','Swap_8_9.m','Swap_4_5.m','Swap_5_9.m','Swap_3_5.m','Swap_3_9.m','Swap_4_9.m','Swap_6_9.m','Swap_5_6.m','Swap_4_6.m','Swap_5_8.m'};
    %Exlcude_1_2_10_11 = {'RemakeOrginal.m','Swap_5_6.m'};
    TitleNames = {};
    for i = 2:length(Exlcude_1_2_10_11)
        Orginal = string(Exlcude_1_2_10_11{1});
        Next = string(Exlcude_1_2_10_11{i});
        PlotListList{i-1} = [Orginal(1),Next(1) ];
        TitleNames{i-1} = "Comparision To Original";
    end
end 

         
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
                            ,"P9","P10","P11","P2P10","As","PI","PE","PF",...
                             "ssDNA","ssPDNA", "RF1","RF2","RF2DP3","P5DNA" ];
    ListOfAnalytesToPlot = [ "DP3", "RNAP", "R", "ssDNA", "ssPDNA", "RF1", "RF2", "RF2DP3", "P5DNA", "DA", "DB", "DH", "DZ", "DW", "EA", "EB", "EH", "EZ", "EW", "ELA", "ELB", "ELH", "ELZ", "ELW", "A", "B", "D", "E", "F", "G", "H", "W", "Y", "Z", "RBS2", "RBS10", "RBS5", "RBS9", "RBS8", "RBS1", "RBS3", "RBS4", "RBS6", "RBS11", "RBS1R", "RBS2R", "RBS3R", "RBS4R", "RBS5R", "RBS6R", "RBS8R", "RBS9R", "RBS10R", "RBS11R", "PD1", "PD2", "PD3", "PD4", "PD5", "PD6", "PD8", "PD9", "PD10", "PD11", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P2P10", "As", "PI", "PE", "PF", "Phage"];

    ColAnalyte = 0;
    for Analyte = ListOfAnalytesToPlot
        ColAnalyte = ColAnalyte + 1; 
        OrginalAnalyte = GetSimulatedData(Analyte,OrginalModelObservables,OrginalModelInfo);
        Samples(1,ColAnalyte) = OrginalAnalyte(end);
        figure 
        ax = axes;
        ax.ColorOrder = linspecer(length(PlotList));
        %ax.FontSize = 40; 
        chosenLinestyle = {'-.','-','--',':'};
        if PlotOrginalOnEveryPlot
            if Log10Plot
                plot(timepointsMinutes,log10(abs(OrginalAnalyte)),"--","DisplayName","Original","LineWidth",4)
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
            LinstyleToPlot = chosenLinestyle{rem(i,4)+1};
            lineWidth=2;
            if Log10Plot
                plot(timepointsMinutes,log10(abs(NextAnalytesToPlot)),LinstyleToPlot,"DisplayName",ModelToRun,"LineWidth",lineWidth)
            else 
                plot(timepointsMinutes,NextAnalytesToPlot,LinstyleToPlot,"DisplayName",ModelToRun,"LineWidth",lineWidth)
            end 
        end 
        %legend('Location','NorthEastOutside','Interpreter', 'none',"FontSize",15)
        otherFontSizes = 20;
        legend('Interpreter', 'none',"FontSize",15)
        grid
        xlabel("Time [Minutes]","FontSize",otherFontSizes)
        if Log10Plot
            ylabel("Log_{10}("+Analyte+")","FontSize",otherFontSizes)
            ylim([0 inf])
        else 
            ylabel(Analyte,"FontSize",otherFontSizes)
            ylim([0 inf])
        end 
        if exist("TitleNames")
            %title(TitleNames(ClusterNum),"FontSize",otherFontSizes)
        else
            title("Cluster Number: " + num2str(ClusterNum-1),"FontSize",otherFontSizes)
        end 
        SAVDIR = OutputDirectory + "/";
        Specific = "Cluster_"+num2str(ClusterNum-1) + "___" + Analyte + ".png";
        %Specific = PlotList{2} + "___" + Analyte + ".png";
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