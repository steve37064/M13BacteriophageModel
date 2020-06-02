clc;  clear; 
set(0,'DefaultFigureWindowStyle','docked'); 
close all; 

%Path to RuleBenderDistribution 
%/Applications/RuleBender.app/Contents/eclipse/BioNetGen-2.5.0/include/ 



num_of_samples = 50; 
sigma = 1; 
ModelName = "Model";
Save_Directory = "TempSaveDirectory";

num_of_mutations_list = [1];

for num_of_mutations = num_of_mutations_list
disp("Working on simulation with " + num2str(num_of_mutations) + " number of mutations")
%==========================================================================
%                          Generate Model Mutations 
%==========================================================================
[param_labels,parameters,observable_labels] = GetModelInformation(ModelName);


Parameter_Mutations = repmat(parameters,num_of_samples,1) ; 
LocationOfMutations = zeros(num_of_samples,length(parameters));

for i_samples = 1:num_of_samples 
    [~,mutation] = datasample(param_labels,num_of_mutations);
    LocationOfMutations(i_samples,mutation) = 1;
    for i_mutation = mutation
        Parameter2Mutate = log10(Parameter_Mutations(i_samples,i_mutation)); 
        MutatedParameter = 10^(normrnd(Parameter2Mutate,sigma));
        Parameter_Mutations(i_samples,mutation) = MutatedParameter ;
    end 
end 

DataObservation = zeros(num_of_samples,length(observable_labels));
FailedSimultions = [];
timepoints = linspace(0,60*60,1000); 
for i_SimulateModel = 1:num_of_samples
    [~, ~, ~, observables_out] = Model( timepoints', [], Parameter_Mutations(i_SimulateModel,:), 1 );
    try
        DataObservation(i_SimulateModel,:) = observables_out(end,:);
    catch 
        FailedSimultions = [FailedSimultions i_SimulateModel]; 
    end 
    %disp(i_SimulateModel)
end 

[~, ~, ~, orginal_observables] = Model( timepoints', [], Parameter_Mutations(i_SimulateModel,:), 1 );
orginal_observables = orginal_observables(end,:);


DataObservation(FailedSimultions,:) = []; 
Parameter_Mutations(FailedSimultions,:) = []; 

%Remove simulation with negative values
BadSimulations = sum(DataObservation<-1,2)~=0;
DataObservation(BadSimulations,:) = [];
Parameter_Mutations(BadSimulations,:) = [];

%==========================================================================
%                          PCA Analysis 
%==========================================================================
G = (orginal_observables - DataObservation + 1)./(1 + orginal_observables);
G = DataObservation./(max(DataObservation)+1);

[coefficents,score,~] = pca(G);

PC1 = score(:,1);
PC2 = score(:,2);

Phage = DataObservation(:,end);
figure 
pcaFigure = scatter(PC1,PC2,[],Phage,'filled');
colorbar
xlabel("Principal Component 1","FontSize",20) 
ylabel("Principal Component 2","FontSize",20)
title("Number of Samples: " + num2str(num_of_samples) + ... 
      " Number of Mutations: " + num2str(num_of_mutations),"FontSize",15)
drawnow

VectorLengths = sqrt(PC1.^2 + PC2.^2);
[values,idx] = max(VectorLengths);


%==========================================================================
%                     Select Point for Analysis 
%==========================================================================
ReRun = 1; 
selectionNumber = 0; 
disp("Click on the simulations you wish to study further and click enter")
disp("If you click multiple times before pressing enter, only the last click"+ ...
     "will be used to record")
disp("Press enter twice to exit and save")
while ReRun
    selectionNumber = selectionNumber+1;
    [Coordinates] = ginput;

    if isempty(Coordinates)
        ReRun = 0; 
        saveas(pcaFigure,Save_Directory+"/PCA____Samples_"+num2str(num_of_samples)+"_Mutations_"+num2str(num_of_mutations)+".png");
    else 
        %Default is the last click: 
        x = Coordinates(end,1);
        y = Coordinates(end,2); 
        
        %Euclidean distance: 
        [~,IDX] = min(sqrt((PC1-x).^2 + (PC2-y).^2));
        
        hold on 
        scatter(PC1(IDX),PC2(IDX),80)
        text(PC1(IDX),PC2(IDX),"\leftarrow "+num2str(selectionNumber))
    end 
end 

end 





function [param_labels,parameters,observable_labels] = GetModelInformation(model)
    % Get parameter names from param_labels string definition
    [status, result]= system(sprintf('grep "param_labels = { " %s.m', model));
    eval(result);
    % Get default parameter values from model 
    [status, result]= system(sprintf('grep "parameters = \\[ " %s.m', model));
    eval(result);
    [status, result]= system(sprintf('grep "observable_labels = { " %s.m', model));
    eval(result);
end 