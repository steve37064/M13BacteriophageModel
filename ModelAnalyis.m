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
%==========================================================================
%        Running the Model and Plotting [Using Default Parameters]
%==========================================================================
%Get parameter and obserable label list and names
global param_labels observable_labels
[param_labels,observable_labels] = GetParamNameAndLabels("Model.m");

LengthToSimulation = 1*60*5; %In Seconds 
timepoints = linspace(0,LengthToSimulation,1000)';
[err, timepoints, species_out, observables_out] = Model(timepoints); 

ViralDNA = ["ssDNA","ssPDNA","RF1"];
PlotGroup("Viral DNA",ViralDNA,timepoints,observables_out,PlotViralDNA)

HostResource = ["DP3"];
PlotGroup("Host Resources",HostResource,timepoints,observables_out,PlotHostResources)

FreePromoterSites = ["DA","DB","DH","DZ","DW"];
PlotGroup("Free Promoter Sites",FreePromoterSites,timepoints,observables_out,PlotPromoterSites)

OccupiedPromoterSites = ["EA"];
PlotGroup("Occupied Promoter Sites",OccupiedPromoterSites,timepoints,observables_out,PlotOccupiedPromoterSites)

ActivelyTranscribedPromoter = ["ELA"];
PlotGroup("Actively Transcribed mRNA",ActivelyTranscribedPromoter,timepoints,observables_out,PlotActivelyTranscribed)

FreeRBSSites = ["RBS2","RBS5","RBS7","RBS9","RBS8"];
PlotGroup("Free Ribosome Binding Sites",FreeRBSSites,timepoints,observables_out,PlotFreeRBSSites)

mRNA = ["A","D"];
PlotGroup("mRNA Strands",mRNA,timepoints,observables_out,PlotmRNA)


%==========================================================================
%Plotting function that accepts a list of species names and plots the list
%on a single graph. 
%==========================================================================
function PlotGroup(GroupType,SpeciesList,Timepoints,ObservablesOut,ShowPlot)
    global observable_labels
    if ShowPlot
        %plottingIndex = FindIndexToPlot(SpeciesList)
        [tf,idx] = ismember(SpeciesList,observable_labels);
        if all(tf)==1 
            figure()
            hold on 
            for i = 1:length(idx)
                IndexToPlot = idx(i);
                plot(Timepoints,ObservablesOut(:,IndexToPlot),'DisplayName',observable_labels{IndexToPlot},"linewidth",4);
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