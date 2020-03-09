clc;  clear; 
set(0,'DefaultFigureWindowStyle','docked'); 
close all; 


%==========================================================================
%            Choose What you want to plot [Default All Species]
%==========================================================================
%Setting the value to 1 will plot the group, otherwise it will not be
%plotted 
PlotHostResources = 1; 
PlotPromoterSites = 1; 
PlotViralDNA = 1; 
PlotOccupiedPromoterSites = 1;
PlotActivelyTranscribed = 1; 
PlotFreeRBSSites = 1 ;
PlotmRNA = 1;
PlotViralRegulatoryProteins = 1; 
%==========================================================================
%        Running the Model and Plotting [Using Default Parameters]
%==========================================================================
%Get parameter and obserable label list and names
global param_labels observable_labels observables_out timepoints
[param_labels,observable_labels] = GetParamNameAndLabels("Model.m");

LengthToSimulation = 1*60*60; %In Seconds 
timepoints = linspace(0,LengthToSimulation,1000)';
[err, timepoints, species_out, observables_out] = Model(timepoints); 
close all 
ViralDNA = ["ssDNA","ssPDNA","RF1","RF2"];
PlotGroup("Viral DNA",ViralDNA,PlotViralDNA)

PlotGroup("Host Resources","DP3",PlotHostResources)
PlotGroup("Host Resources","RNAP",PlotHostResources)
PlotGroup("Host Resources","R",PlotHostResources)

FreePromoterSites = ["DA","DB","DH","DZ","DW"];
PlotGroup("Free Promoter Sites",FreePromoterSites,PlotPromoterSites)

OccupiedPromoterSites = ["EA"];
PlotGroup("Occupied Promoter Sites",OccupiedPromoterSites,PlotOccupiedPromoterSites)

ActivelyTranscribedPromoter = ["ELA"];
PlotGroup("Actively Transcribed mRNA",ActivelyTranscribedPromoter,PlotActivelyTranscribed)

FreeRBSSites = ["RBS2","RBS5","RBS9","RBS8"];
PlotGroup("Free Ribosome Binding Sites",FreeRBSSites,PlotFreeRBSSites)

mRNA = ["A","B","D","E","F","G","H"];
PlotGroup("mRNA Strands",mRNA,PlotmRNA)

mRNA = ["Z","Y","W"];
PlotGroup("mRNA Strands",mRNA,PlotmRNA)

regulatoryProteins = ["P2","P10","P5","P2P10"];
PlotGroup("Viral Proteins",regulatoryProteins,PlotViralRegulatoryProteins)

regulatoryProteins = ["P8","P7","P9","P3","P6"];
PlotGroup("Viral Proteins",regulatoryProteins,PlotViralRegulatoryProteins)

regulatoryProteins = ["P1","P4","P11"];
PlotGroup("Viral Proteins",regulatoryProteins,PlotViralRegulatoryProteins)

%==========================================================================
%Plotting function that accepts a list of species names and plots the list
%on a single graph. 
%==========================================================================
function PlotGroup(GroupType,SpeciesList,ShowPlot)
    global observable_labels observables_out timepoints
    if ShowPlot
        %plottingIndex = FindIndexToPlot(SpeciesList)
        [tf,idx] = ismember(SpeciesList,observable_labels);
        if all(tf)==1 
            figure()
            hold on 
            for i = 1:length(idx)
                IndexToPlot = idx(i);
                plot(timepoints,observables_out(:,IndexToPlot),'DisplayName',observable_labels{IndexToPlot},"linewidth",4);
            end
            legend
            xlabel("Time [Seconds]",'FontSize',20)
            ylabel("Number of Molecules",'FontSize',20) 
            title(GroupType,'FontSize',20)
            grid()
            hold off 
            
        else 
            disp("One of your species in the list is not contained in the obserables. Check the "+GroupType +" list to ensure no mistakes are made");
        end 
    else 
        disp(GroupType+" plot was supressed");
    end 
    
end


%==========================================================================
%Function to scan the bngl file and get the order in which the observables 
%and parameters are defined. This was done to make this file robust against
%changes to the BNGL model. 
%==========================================================================
function [param_labels,observable_labels] = GetParamNameAndLabels(mFileModelLoc)

    modelFile = fopen(mFileModelLoc);
    SearchTerms = ["param_labels = " "observable_labels = "];

    FoundLines = cell(1,length(SearchTerms)); 
    n = 1;
    while true  
        TestLine = fgetl(modelFile);
        if contains(TestLine,SearchTerms(n))
            FoundLines{n} = TestLine;
            n =n+1;
        end 
        if TestLine == -1 | n == 3
            break 
        end 
    end 
    clear TestLine
    fclose(modelFile);

    param_labels      = FindNames(FoundLines{1});
    observable_labels = FindNames(FoundLines{2});



    function output = FindNames(StringInput)
        StringInput = replace(StringInput," ","");
        StringInput = replace(StringInput,"'","");
        StringInput = strsplit(StringInput,["{","}"]);
        StringInput = StringInput{2};
        StringInput = strsplit(StringInput,",");
        output = StringInput;
    end 
end 