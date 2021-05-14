%--------------------------------------------------------------------------
%Model Run------------------Version 3----------------------February 27,2014
%--------------------------------------------------------------------------
clc; clear all 
close all 
%--------------------------------------------------------------------------
%-------------------------What to Plot-------------------------------------
%--------------------------------------------------------------------------
Plot_What = zeros(1,12); 
Plot_What(1)  = 0;    %All DNA species will be plotted 
Plot_What(2)  = 0;    %Free RNAP promoter sites  
Plot_What(3)  = 0;    %RNAP promoter sites bound to a RNA polymerase 
Plot_What(4)  = 0;    %Elongating RNAP away from promoter site 
Plot_What(5)  = 0;    %Polycistronic mRNA species 
Plot_What(6)  = 0;    %Free Ribosome Binding sites 
Plot_What(7)  = 0;    %Ribosome binding site with ribosome attached 
Plot_What(8)  = 0;    %Elongating Ribosomes away from ribosome binding site 
Plot_What(9)  = 0;    %Phage Proteins 
Plot_What(10) = 1;    %New Phage  
Plot_What(11) = 0;    %Assembly Site formation species 
Plot_What(12) = 0;    %E.Coli Proteins 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%----------------------------Inputs----------------------------------------
%--------------------------------------------------------------------------

%------------------------------General-------------------------------------
%-----Defining Which Model to Use .sproj file
Load_Model = 'Final_Model.sbproj';
%-----The length of Simulation Time in seconds 
Simulation_Time = 3600*1; 
%-----Paramters That May be Varied 
P = zeros(1,66);     %Empty Vector that will contain all constant variables 
%-----Paramters That need to be converted to rate constants 
PC = zeros(1,16);    %Empty Vector  
%--------------------------------------------------------------------------

%------------------------------Replication---------------------------------
P(1) = 1;   %Forward Reaction Rate for P2 nicking DNA (1/(molecule*second))
P(2) = P(1)*1000;      %Reverse Reaction Rate for P2 nicking DNA (1/second) 
P(3) = 1;             %P2/P10 foward binding constant (1/(molecule*second))
P(4) = 9.5454e-16;    %P2/P10 reverse binding constant (1/(second))
P(5) = .001;          %Rate of P5 Sequestering ssDNA
P(6) = 1100;     %Number of P5 molcules needed to start S.S. DNA Production
P(7) = 0.1;           %Rate of DNA polymerase Binding 
PC(1) = 800;          %Rate of DNA polymerase elongation (Nucleotides/S)
%--------------------------------------------------------------------------

%------------------------------Transcrption--------------------------------
PC(2) = 0.01;           %RNA Polymerase Binding (1/(molecule*second))
PC(3) = 67;          %Rate 0f RNAP elongation (Nucleotides/S)
%--------------------------------------------------------------------------

%------------------------------Translation---------------------------------
PC(4) = 20;          %Ribosome Polymerase Binding (1/(molecule*second))
PC(5) = 14*3;        %Ribosome Elongation Rate (nucleotides/second)
P(61) = 0.00015;     %Efficency of P5 Inhibition to P2 
P(62) = 0.000026;    %Efficency of P5 Inhibition to P10 
P(63) = 0.01;        %Efficency of P5 Inhibition to P5 
P(64) = 4;           %Efficency of P5 Inhibition to P3 
P(65) = 80;           %Efficency of P5 Inhibition to P1 
%Note-An efficency factor less than 1, will increase the inhibition of P5
%on translation. An efficency factor greater than one will decrease the
%effect P5 has on inhibition. 
%--------------------------------------------------------------------------

%------------------------------mRNA Degragation----------------------------
PC(6)  = 3.4/6;          %Half-Life of mRNA A in minutes 
PC(7)  = 2;            %Half-Life of mRNA B in minutes 
PC(8)  = 2;            %Half-Life of mRNA C in minutes 
PC(9)  = 2;            %Half-Life of mRNA D in minutes 
PC(10) = 2.5;          %Half-Life of mRNA E in minutes 
PC(11) = 6;            %Half-Life of mRNA F in minutes 
PC(12) = 8;            %Half-Life of mRNA G in minutes 
PC(13) = 10;           %Half-Life of mRNA H in minutes 
PC(14) = 6;            %Half-Life of mRNA Z in minutes 
PC(15) = 6;            %Half-Life of mRNA Y in minutes 
PC(16) = 6;            %Half-Life of mRNA W in minutes 
%--------------------------------------------------------------------------

%------------------------Phage Production----------------------------------

P(55) = 0.00123;      %Rate of new phage elongation (1/s) 
P(56) = 0.1;           %Rate of assembly site formation (1/(Second*Molecule)) 
P(57) = 0.1;           %Rate of P7/P9 forming a complex with an assembly site (1/(Second*Molecule^2)) 
P(58) = 0.1;           %Rate of P5DNA binding to assembly site containt p7 and P9 (1/(Second*Molecule)) 
P(59) = 0.1;          %Rate of Phage Termination (1/(Second*Molecule^2)) 
P(60) = 3000;         %Minium number of P8 molecules needed to have phage producing at full satuatino (Molecule) 
%--------------------------------------------------------------------------

%-------------------Hill Coefficent Values---------------------------------
P(66) = 40;  %This one makes sure no negative P8 is being made, that is why it is so high   
P(67) = 1;   %This is the one that dictates the P5 inhibition reactions 
P(68) = 1;   %This the the hill coefficent constant for the rate of P5DNA formation
%--------------------------------------------------------------------------

%-----------------Initial Amount of E.Coli Proteins------------------------
Initial_Amount = [ 
    3      %DNA Polymerase  
    7880   %Ribosomes 
    1280   %RNA Polymerases  [1500-11400 Total]
    ]; 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------Loading Model and Parameters--------------------
%--------------------------------------------------------------------------
%-----Loading the Model 
sbioloadproject(Load_Model) 
%-----Defining the Legnth of the Simulation Time 
m1.Configset.StopTime = Simulation_Time;
%-----Calculation of the constants that are under the PC array 
[P(8:54)] = Load_Parameters(PC); 
%-----Loading the Parameters  
for n = 1:length(P) 
    m1.parameters(n).value = P(n); 
end 
%----Loading the Intial Amount of E.Coli Species 
for n = 75:77
    m1.species(n).InitialAmount = Initial_Amount(n-74);
end 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------Simulating The Model----------------------------
%--------------------------------------------------------------------------

[t,x] = sbiosimulate(m1); 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------Plotting Species--------------------------------
%--------------------------------------------------------------------------
t = t./60; 

Plot_Phage_Species(Plot_What,t,x)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------