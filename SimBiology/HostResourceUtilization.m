

%--------------------------------------------------------------------------
%Calculating Total Polmerases Elongating mRNA A 
PEA = x(:,12) + x(:,17); %Polymerases Elongating A 
PEB = x(:,13) + x(:,18); %Polymerases Elongating B  
PEH = x(:,14) + x(:,19); %Polymerases Elongating H 
PEZ = x(:,15) + x(:,20); %Polymerases Elongating Z 
PEW = x(:,16) + x(:,21); %Polymerases Elongating W 

Total = PEA + PEB + PEH + PEZ + PEW + x(:,77);
figure
plot(t,PEA./Total.*100) 
hold on 
plot(t,PEB./Total.*100,'--') 
plot(t,PEH./Total.*100,'-.') 
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Total RNAP Elongating an mRNA') 
legend('mRNA A','mRNA B','mRNA H') 

figure 
plot(t,PEZ./Total.*100) 
hold on 
plot(t,PEW./Total.*100,'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Total RNAP Elongating an mRNA') 
legend('mRNA Z','mRNA W') 

figure 
plot(t,x(:,77)./Total.*100)
xlabel('Time = Minutes') 
ylabel('Percent of Free RNAP') 
grid 

figure 
plot(t,x(:,7)./(x(:,7) + x(:,12)).*100)
hold on 
plot(t,x(:,8)./(x(:,8) + x(:,13)).*100,'--')
plot(t,x(:,9)./(x(:,9) + x(:,14)).*100,'-.')
legend('DA','DB','DH') 
grid 
ylim([0 75])
xlabel('Time = Minutes') 
ylabel('Percent of Free Promoter Sites')

figure
plot(t,x(:,10)./(x(:,10) + x(:,15)).*100)
hold on 
plot(t,x(:,11)./(x(:,11) + x(:,16)).*100,'--')
legend('DZ','DW') 
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Free Promoter Sites')




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Calculating the percentages of ribosmes actively translating various 
%proteins
RP2  = x(:,43) + x(:,53); 
RP10 = x(:,44) + x(:,54); 
RP5  = x(:,45) + x(:,55); 
RP9  = x(:,46); 
RP8  = x(:,47) + x(:,56); 
RP3  = x(:,48) + x(:,57); 
RP6  = x(:,49) + x(:,58); 
RP1  = x(:,50) + x(:,59); 
RP11 = x(:,51) + x(:,60); 
RP4  = x(:,52) + x(:,61); 

Total = RP2+RP10+RP5+RP9+RP8+RP3+RP6+RP1+RP11+RP4+x(:,76)

figure
plot(t,RP2./Total.*100) 
hold on 
plot(t,RP10./Total.*100,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating a Protein ') 
legend('P2','P10')
grid 

figure
plot(t,RP5./Total.*100) 
hold on 
plot(t,RP3./Total.*100,'--')
plot(t,RP1./Total.*100,'-.')

xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating a Protein ') 
legend('P5','P3','P1')
grid 

figure
plot(t,RP8./Total.*100) 
hold on 
plot(t,RP9./Total.*100,'--')
plot(t,RP4./Total.*100,'-.')

xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating a Protein ')
legend('P8','P9','P4')
grid 

figure
hold on 
plot(t,RP11./Total.*100)
plot(t,RP6./Total.*100,'-.') 
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating a Protein ') 
legend('P11','P6')
grid 

figure
hold on 
plot(t,x(:,76)./Total.*100) 
xlabel('Time = Minutes')
ylabel('Percent of Free Ribosomes') 
grid 
%--------------------------------------------------------------------------




































