function [Parameters] = Load_Parameters(Master_Variables) 

%This function will recalculate all parameters from a set of Master
%variables. The Master_Variables will be an arrow consisting of the
%following values: 
%Master_Variables = [DNA Polymerase Elongation Rate(Nuc/s)     1
%                    RNA Polymerase Binding                    2
%                    RNA Polymerase Elongation Rate(Nuc/s)     3
%                    Ribosome Binding                          4
%                    Ribosome Elongation Rate(Nuc/s)           5
%                    mRNA A Degragation Half Life (minutes)    6
%                    mRNA B Degragation Half Life (minutes)    7
%                    mRNA C Degragation Half Life (minutes)    8
%                    mRNA D Degragation Half Life (minutes)    9
%                    mRNA E Degragation Half Life (minutes)    10
%                    mRNA F Degragation Half Life (minutes)    11
%                    mRNA G Degragation Half Life (minutes)    12
%                    mRNA H Degragation Half Life (minutes)    13
%                    mRNA Z Degragation Half Life (minutes)    14
%                    mRNA Y Degragation Half Life (minutes)    15
%                    mRNA W Degragation Half Life (minutes)    16


%--------------------------------------------------------------------------
%-------------------Constants that Can Changes-----------------------------
%--------------------------------------------------------------------------

%--------------------------------Replication-------------------------------
Length_Of_Genome = 6400;            %Nucleotides 
%--------------------------------------------------------------------------

%-----------------------------Transcription--------------------------------
PBSA = 0.6;      %Promoter Binding Stregnth of mRNA A Range 0.6-0.9 
PBSB = 1.00;     %Promoter Binding Stregnth of mRNA B Range 1.00
PBSH = 0.65;     %Promoter Binding Stregnth of mRNA H Range 0.5-0.8
PBSZ = 0.02;     %Promoter Binding Stregnth of mRNA Z Range 0.15 - 0.30
PBSW = 0.225;    %Promoter Binding Stregnth of mRNA W Range 0.15 - 0.30
RNA_Length = 50; %Number of basepairs RNA polymerase needs to travel to 
                 %clear the promoter site 
LA = 1900;  %Nucleotide Length of Transcript A 
LB = 1000;  %Nucleotide Length of Transcript B 
LH = 320;   %Nucleotide Length of Transcript H 
LZ = 1600;  %Nucleotide Length of Transcript Z 
LW = 1270;  %Nucleotide Length of Transcript W 
%--------------------------------------------------------------------------

%------------------------------Translation---------------------------------
RBSP2 = 0.54;          %Ribosome Binding Site Strength of P2 
RBSP10 = 0.054;        %Ribosome Binding Site Strength of P10
RBSP5 = 1.00;          %Ribosome Binding Site Strength of P5
RBSP9 = 0.08;          %Ribosome Binding Site Strength of P9
RBSP8 = 0.22;          %Ribosome Binding Site Strength of P8
RBSP3 = 0.40;          %Ribosome Binding Site Strength of P3
RBSP6 = 0.04;          %Ribosome Binding Site Strength of P6
RBSP1 = 0.08;          %Ribosome Binding Site Strength of P1
RBSP11 = 0.20;         %Ribosome Binding Site Strength of P11
RBSP4 = 0.23;          %Ribosome Binding Site Strength of P4
RSR = 173;             %Ribosome Spacing Requirment (Nucleotides) 
LOP2  = 1200;          %Nucleotide Length of P2 Coding Region  
LOP10 =  350;          %Nucleotide Length of P10 Coding Region
LOP5  =  270;          %Nucleotide Length of P5 Coding Region
LOP8  =  220;          %Nucleotide Length of P8 Coding Region
LOP3  = 1300;          %Nucleotide Length of P3 Coding Region
LOP6  =  340;          %Nucleotide Length of P6 Coding Region
LOP1  = 1040;          %Nucleotide Length of P1 Coding Region
LOP11 = 310;           %Nucleotide Length of P11 Coding Region
LOP4  = 1270;          %Nucleotide Length of P4 Coding Region
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%----------------------------Calculations----------------------------------
k_Replication = Master_Variables(1)/Length_Of_Genome;   %1/s 
Ka = Master_Variables(2)*PBSA;                          %1/(s*molecule)      
Kb = Master_Variables(2)*PBSB;                          %1/(s*molecule)
Kh = Master_Variables(2)*PBSH;                          %1/(s*molecule)
Kz = Master_Variables(2)*PBSZ;                          %1/(s*molecule)
Kw = Master_Variables(2)*PBSW;                          %1/(s*molecule)
RNAP_Sat = Master_Variables(3)/RNA_Length;              %1/s
RNAP_A = Master_Variables(3)/(LA-RNA_Length);           %1/s
RNAP_B = Master_Variables(3)/(LB-RNA_Length);           %1/s
RNAP_H = Master_Variables(3)/(LH-RNA_Length);           %1/s
RNAP_Z = Master_Variables(3)/(LZ-RNA_Length);           %1/s
RNAP_W = Master_Variables(3)/(LW-RNA_Length);           %1/s
KP2  = Master_Variables(4)*RBSP2;                       %1/(s*molecule)
KP10 = Master_Variables(4)*RBSP10;                      %1/(s*molecule)                       
KP5  = Master_Variables(4)*RBSP5;                       %1/(s*molecule)
KP9  = Master_Variables(4)*RBSP9;                       %1/(s*molecule)
KP8  = Master_Variables(4)*RBSP8;                       %1/(s*molecule)
KP3  = Master_Variables(4)*RBSP3;                       %1/(s*molecule)
KP6  = Master_Variables(4)*RBSP6;                       %1/(s*molecule)
KP1  = Master_Variables(4)*RBSP1;                       %1/(s*molecule)
KP11 = Master_Variables(4)*RBSP11;                      %1/(s*molecule)
KP4  = Master_Variables(4)*RBSP4;                       %1/(s*molecule)
Rib_Sat = Master_Variables(5)/RSR;                      %1/s 
KEL2  = Master_Variables(5)/(LOP2-RSR);                 %1/s
KEL10 = Master_Variables(5)/(LOP10-RSR);                %1/s
KEL5  = Master_Variables(5)/(LOP5-RSR);                 %1/s
KEL8  = Master_Variables(5)/(LOP8-RSR);                 %1/s
KEL3  = Master_Variables(5)/(LOP3-RSR);                 %1/s
KEL6  = Master_Variables(5)/(LOP6-RSR);                 %1/s
KEL1  = Master_Variables(5)/(LOP1-RSR);                 %1/s
KEL11 = Master_Variables(5)/(LOP11-RSR);                %1/s
KEL4  = Master_Variables(5)/(LOP4-RSR);                 %1/s
Ads = log(2)/(Master_Variables(6)*60);                  %1/s 
Bds = log(2)/(Master_Variables(7)*60);                  %1/s 
Cds = log(2)/(Master_Variables(8)*60);                  %1/s 
Dds = 0.7*log(2)/(Master_Variables(9)*60);              %1/s 
Dds_Alt = 0.3*log(2)/(Master_Variables(9)*60);          %1/s 
Eds = 0.7*log(2)/(Master_Variables(10)*60);             %1/s 
Eds_Alt = 0.3*log(2)/(Master_Variables(10)*60);         %1/s 
Fds = 0.7*log(2)/(Master_Variables(11)*60);             %1/s 
Fds_Alt = 0.3*log(2)/(Master_Variables(11)*60);         %1/s 
Gds = 0.7*log(2)/(Master_Variables(12)*60);             %1/s 
Gds_Alt = 0.3*log(2)/(Master_Variables(12)*60);         %1/s 
Hds = log(2)/(Master_Variables(13)*60);                 %1/s 
Zds = log(2)/(Master_Variables(14)*60);                 %1/s
Yds = log(2)/(Master_Variables(15)*60);                 %1/s
YWs = log(2)/(Master_Variables(16)*60);                 %1/s

%--------------------------------------------------------------------------
Parameters = [
k_Replication 
Ka    
Kb 
Kh 
Kz 
Kw 
RNAP_Sat 
RNAP_A 
RNAP_B 
RNAP_H 
RNAP_Z 
RNAP_W 
KP2  
KP10                     
KP5  
KP9 
KP8  
KP3  
KP6 
KP1 
KP11 
KP4  
Rib_Sat 
KEL2  
KEL10 
KEL5  
KEL8  
KEL3 
KEL6  
KEL1  
KEL11 
KEL4  
Ads 
Bds 
Cds 
Dds 
Dds_Alt
Eds 
Eds_Alt
Fds 
Fds_Alt
Gds 
Gds_Alt
Hds 
Zds 
Yds 
YWs ]; 



