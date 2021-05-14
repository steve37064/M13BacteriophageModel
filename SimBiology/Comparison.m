%Comparison of Cellular Resources Usage With and Withou P5 Inhibtion 

%These Graphs Are Useful when comparing how cellular resources shift when
%P5 does not inhibit P2 Translation. 

%--------------------------------------------------------------------------
%Calculating the number  of ribosmes actively translating various 
%proteins when P5 is turned on for all proteins 
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
Total = RP2+RP10+RP5+RP9+RP8+RP3+RP6+RP1+RP11+RP4+x(:,76);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Calculating the number  of ribosmes actively translating various 
%proteins when P5 is turned off for a single protein
RP2_1  = x1(:,43) + x1(:,53); 
RP10_1 = x1(:,44) + x1(:,54); 
RP5_1  = x1(:,45) + x1(:,55); 
RP9_1  = x1(:,46); 
RP8_1  = x1(:,47) + x1(:,56); 
RP3_1  = x1(:,48) + x1(:,57); 
RP6_1  = x1(:,49) + x1(:,58); 
RP1_1  = x1(:,50) + x1(:,59); 
RP11_1 = x1(:,51) + x1(:,60); 
RP4_1  = x1(:,52) + x1(:,61); 
Total_1 = RP2_1+RP10_1+RP5_1+RP9_1+RP8_1+RP3_1+RP6_1+RP1_1+RP11_1+RP4_1 ...
          +x1(:,76);
%--------------------------------------------------------------------------


figure 
plot(t,RP1./Total) 
hold on 
plot(t1,RP1_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P1') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 


figure 
plot(t,RP2./Total) 
hold on 
plot(t1,RP2_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P2') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 

figure 
plot(t,RP3./Total) 
hold on 
plot(t1,RP3_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P3') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 

figure 
plot(t,RP4./Total) 
hold on 
plot(t1,RP4_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P4') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 

figure 
plot(t,RP5./Total) 
hold on 
plot(t1,RP5_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P5') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 

figure 
plot(t,RP6./Total) 
hold on 
plot(t1,RP6_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P6') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 




figure 
plot(t,RP8./Total) 
hold on 
plot(t1,RP8_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P8') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 

figure 
plot(t,RP9./Total) 
hold on 
plot(t1,RP9_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P9') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 

figure 
plot(t,RP10./Total) 
hold on 
plot(t1,RP10_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P10') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 


figure 
plot(t,RP11./Total) 
hold on 
plot(t1,RP11_1./Total_1,'--')
xlabel('Time = Minutes')
ylabel('Percent of Total Ribosomes Translating P11') 
legend('Inhibition Turned On','Inhibition Turned Off')
grid 