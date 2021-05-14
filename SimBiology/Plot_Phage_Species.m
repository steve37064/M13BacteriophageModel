function [] = Plot_Phage_Species(WhichPlot,t,x)

%------------------------Documentation-------------------------------------
%This file is used to plot the designed model species as a function of time
%The Input is a vector containing zeros and ones. If it is equal to 1 then
%the group of species will be plotted. If zero than it will not be plotted.

%The Input is as follows: 
%WhichPlot = [ 
%              All DNA species will be Plotted 
%              Free RNAP promoter sites  
%              RNAP promoter sites bound to a RNA polymerase 
%              Elongating RNAP away from promoter site 
%              Polycistronic mRNA species 
%              Free Ribosome Binding sites 
%              Ribosome binding site with ribosome attached 
%              Elongating Ribosomes away from ribosome binding site 
%              Phage Proteins 
%              New Phage  
%              Assembly Site formation species 
%              E.Coli Proteins ]
%t = simulation time vector (In seconds) 
%x = Matrix of phage species as a function of time  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


xaxis = ['Time [Seconds]'];

%Plotting All DNA Species 
if WhichPlot(1) == 1 
    Species_Name = {'ssDNA';'ssPDNA';'RF1';'RF2';'RF2DP3';'P5DNA'};
    for n = 1:6 
        figure 
        plot(t,x(:,n))
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
       
    end 
end 

%Plotting all free RNAP promoter sites 
if WhichPlot(2) == 1 
    Species_Name = {'DA';'DB';'DH';'DZ';'DW'}; 
    for n = 1:5 
        figure 
        plot(t,x(:,n+6)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 

%Plotting all promoter sites occpied by an RNAP enzyme 
if WhichPlot(3) == 1 
    Species_Name = {'EA';'EB';'EH';'EZ';'EW'};
    for n = 1:5 
        figure 
        plot(t,x(:,n+11)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 

%Plotting all RNAP enzymes that are away from promoter site but still
%elongating on the mRNA 
if WhichPlot(4) == 1 
    Species_Name = {'ELA';'ELB';'ELH';'ELZ';'ELW'};
    for n = 1:5 
        figure 
        plot(t,x(:,n+16)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 

%Plotting All mRNA species 
if WhichPlot(5) == 1 
 Species_Name = {'A';'B';'C';'D';'E';'F';'G';'H';'Z';'Y';'W'};
    for n = 1:11 
        figure 
        plot(t,x(:,n+21)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
       
    end 
end 

%Plotting all free ribosome binding sites 
if WhichPlot(6) == 1 
 Species_Name = {'RBS2';'RBS10';'RBS5';'RBS9';'RBS8';'RBS3';'RBS6';...
                 'RBS1';'RBS11';'RBS4'};
    for n = 1:10 
        figure 
        plot(t,x(:,n+32)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 
%Plottign all ribsosome binidng site with a ribsosome attached to it 
if WhichPlot(7) == 1 
 Species_Name = {'RBS2R';'RBS10R';'RBS5R';'RBS9R';'RBS8R';'RBS3R';...
                 'RBS6R';'RBS1R';'RBS11R';'RBS4R'};
    for n = 1:10 
        figure 
        plot(t,x(:,n+42)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 

%Plotting all ribosome still making protein and bound to mRNA 
if WhichPlot(8) == 1 
 Species_Name = {'PD2';'PD10';'PD5';'PD8';'PD3';'PD6';'PD1';'PD11';'PD4'};
    for n = 1:9 
        figure 
        plot(t,x(:,n+52)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 

%Plotting ALL phage proteins 
if WhichPlot(9) == 1 
 Species_Name = {'P2';'P10';'P5';'P7';'P9';'P8';'P3';'P6';'P1';'P11';...
                 'P4';'P2P10'};
    for n = 1:12 
        figure 
        plot(t,x(:,n+61)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 
%Plotting New Phage 
if WhichPlot(10) == 1 
    Species_Name = {'Phage'};
    for n = 1:1 
        figure 
        plot(t,x(:,n+73)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 

%Plotting Assembly Site Proteins 
if WhichPlot(11) == 1 
    Species_Name = {'Free Assembly Site';'Phage Initiation';...
                    'Phage Elongation'; 'Phage Termination' };
    for n = 1:4 
        figure 
        plot(t,x(:,n+77)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
        
    end 
end 


%Plotting E.Coli Proteins 
if WhichPlot(12) == 1 
    Species_Name = {'DNA Polymerse';'Ribosomes';'RNA Polymerase'};
    for n = 1:3 
        figure 
        plot(t,x(:,n+74)) 
        xlabel(xaxis) 
        ylabel(Species_Name(n)) 
        grid 
       
    end 
    

end 




