%This File will be use to custom generate plots specified for the papers
%needs. 

%--------------------------------------------------------------------------
figure 
plot(t,x(:,76))
hold on 
plot(t,x(:,77),'--')
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
grid 
legend('Number of Free Ribosomes','Number of Free RNA Polymerase 3') 
figure 
plot(t,x(:,75))
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
grid 
legend('Number of Free DNA Polymerase 3') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,78)) 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
grid 
legend('Free Assembly Sites (A_S)') 
figure 
plot(t,x(:,79)) 
hold on 
plot(t,x(:,80),'--')
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
grid 
legend('Phage Initiation (P_I) ','Phage Elongaton (P_E)') 
figure 
plot(t,x(:,81)) 
hold on 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
grid 
legend('Assembly Termination (P_F)') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,51) + x(:,60)) 
hold on 
plot(t,x(:,52) + x(:,61),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Actively Translating Ribosomes') 
legend('P11 Translation','P4 Translation') 
figure 
plot(t,x(:,47) + x(:,46)) 
hold on 
plot(t,x(:,49) + x(:,58),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Actively Translating Ribosomes') 
legend('P8 Translation','P6 Translation') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,48) + x(:,57)) 
hold on 
plot(t,x(:,50) + x(:,59),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Actively Translating Ribosomes') 
legend('P3 Translation','P1 Translation') 
figure 
plot(t,x(:,45) + x(:,55)) 
hold on 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Actively Translating Ribosomes') 
legend('P5 Translation') 
figure 
plot(t,x(:,43) + x(:,53)) 
hold on 
plot(t,x(:,44) + x(:,54),'--')
grid 
xlabel('Time = Minutes') 
ylabel('Number of Actively Translating Ribosomes') 
legend('P2 Translation','P10 Translation') 
%--------------------------------------------------------------------------
figure 
plot(t,x(:,36)./(x(:,36) + x(:,46)))
hold on 
plot(t,x(:,37)./(x(:,37) + x(:,47)),'--')
plot(t,x(:,39)./(x(:,39) + x(:,49)),'-.')
legend('RBS9','RBS8', 'RBS6')
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Free Ribosome Binding Sites ')
figure 
plot(t,x(:,41)./(x(:,41) + x(:,51)))
hold on 
plot(t,x(:,42)./(x(:,42) + x(:,52)),'--')
legend('RBS11','RBS4')
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Free Ribosome Binding Sites ')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,33)./(x(:,33) + x(:,43)))
hold on 
plot(t,x(:,34)./(x(:,34) + x(:,44)),'--')
plot(t,x(:,35)./(x(:,35) + x(:,45)),'-.')
legend('RBS 2','RBS 10', 'RBS5')
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Free Ribosome Binding Sites ')
figure 
plot(t,x(:,38)./(x(:,38) + x(:,48)))
hold on 
plot(t,x(:,40)./(x(:,40) + x(:,50)),'--')
legend('RBS 3','RBS 1')
grid 
xlabel('Time = Minutes') 
ylabel('Percent of Free Ribosome Binding Sites ')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,(sum(x(:,12:21)')')./(x(:,3)+x(:,4)+x(:,5)))
hold on 
grid 
xlabel('Time = Minutes') 
ylabel('Elongating RNAP Per dsDNA')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,(x(:,15) + x(:,20))./(x(:,3)+x(:,4)))
hold on 
plot(t,(x(:,16) + x(:,21))./(x(:,3)+x(:,4)),'--')
legend('Z','W') 
grid 
xlabel('Time = Minutes') 
ylabel('Elongating RNAP Per dsDNA')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,(x(:,12) + x(:,17))./(x(:,3)+x(:,4)))
hold on 
plot(t,(x(:,13) + x(:,18))./(x(:,3)+x(:,4)),'--')
plot(t,(x(:,14) + x(:,19))./(x(:,3)+x(:,4)),'-.')
legend('A','B','H') 
grid 
xlabel('Time = Minutes') 
ylabel('Elongating RNAP Per dsDNA')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,(x(:,12) + x(:,17))./(x(:,3)+x(:,4)))
hold on 
plot(t,(x(:,13) + x(:,18))./(x(:,3)+x(:,4)),'--')
plot(t,(x(:,14) + x(:,19))./(x(:,3)+x(:,4)),'-.')
plot(t,(x(:,15) + x(:,20))./(x(:,3)+x(:,4)),'r')
plot(t,(x(:,16) + x(:,21))./(x(:,3)+x(:,4)),'--r')
legend('A','B','H','Z','W') 
grid 
xlabel('Time = Minutes') 
ylabel('Elongating RNAP Per dsDNA')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,32)) 
hold on 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('mRNA W') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,30)) 
hold on 
plot(t,x(:,31),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('mRNA Z','mRNA Y') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,28)) 
hold on 
plot(t,x(:,29),'--')
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('mRNA G','mRNA H') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,25)) 
hold on 
plot(t,x(:,26),'--')
plot(t,x(:,27),'-.')
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('mRNA D','mRNA E','mRNA F') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,22)) 
hold on 
plot(t,x(:,23),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('mRNA A','mRNA B') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,(x(:,74))) 
hold on 
plot(t,(x(:,80)),'--')
grid 
xlabel('Time = Minutes') 
ylabel('Number of Particles')
legend('Phage','Elongating Phage') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,log10(x(:,70) + 14.*(x(:,78) + x(:,79) + x(:,80) + x(:,81)))) 
hold on 
plot(t,log10(x(:,71)+ 14.*(x(:,78) + x(:,79) + x(:,80) + x(:,81))),'--') 
plot(t,log10(x(:,72)+ 14.*(x(:,78) + x(:,79) + x(:,80) + x(:,81))),'-.')
grid 
xlabel('Time = Minutes') 
ylabel('log_1_0(Number of Molecules)')
ylim([ 0 6])
legend('Total P1','Total P11','Total P4') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,x(:,67)) 
hold on 
plot(t,x(:,67) + 2700.* x(:,74),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('Total P8 in the Cell','Total P8 Synthesized') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,68)) 
hold on 
plot(t,x(:,69),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('Free P3','Free P6') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,68) + 5.*x(:,74)) 
hold on 
plot(t,x(:,69)+ 5.*x(:,74),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('Total P3 Produced','Total P6 Produced') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
plot(t,log10(x(:,65)))
hold on 
plot(t,log10(x(:,66)),'--') 
grid 
xlabel('Time = Minutes') 
ylabel('log_1_0(Number of Molecules)')
legend('Free P7','Free P9') 
%--------------------------------------------------------------------------
% --------------------------------------------------------------------------
figure 
plot(t,x(:,62))
hold on 
plot(t,x(:,63),'--') 
plot(t,x(:,73),'-.')
grid 
xlabel('Time = Minutes') 
ylabel('Number of Molecules')
legend('Free P2','Free P10','P2/P10 Complex') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,log10(x(:,64)))
hold on 
plot(t,log10(1600.*(x(:,6)+ x(:,80) + x(:,81))),'r')
plot(t,log10(x(:,64) + 1600.*(x(:,6)+ x(:,80) + x(:,81))),'k')
xlabel('Time = Minutes') 
ylabel('Log_1_0(Number of Molecules)') 
legend('Free P5','Bound P5','Total P5')  
ylim([0 6])
grid 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,6)+ x(:,80) + x(:,81))
hold on 
plot(t,x(:,79) + x(:,80) + x(:,81),'--') 
xlabel('Time = Minutes') 
ylabel('Number of Molecules') 
grid 
legend('Total P5DNA','Total Assembly Sites') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,3)) 
hold on 
plot(t,x(:,4),'--') 
plot(t,x(:,5),'-.') 
grid
xlabel('Time = Minutes') 
ylabel('Number of Molecules') 
legend('RF1 DNA', 'RF2 DNA','RF2/DP3 Complex') 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure 
plot(t,x(:,1)) 
hold on 
plot(t,x(:,2),'--') 
grid
xlabel('Time = Minutes') 
ylabel('Number of Molecules') 
legend('Free ssDNA Strands','ssDNA/DP3 Complex') 
%--------------------------------------------------------------------------

