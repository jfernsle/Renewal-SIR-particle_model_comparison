%% Simulation and model statisitcs
%This code combines multiple particle simulation runs to find averages and
%standard deviations.  Each simulation is modeled using Renewal model and
%SIR models. The average of each curve is plotted, along with a shaded area
%showing the standard deviation.
%Example particle simulation data sets are included in the Github
%repository. The code uses option 'D' to read gi(tau) from generation
%intervals measured from the first 15 days of particle simulations and fits
%the distribution to a linear function. Options A-C are therefore commented
%out of the code. If using a different type of simulation, uncomment the
%appropriate section.

clc
clear all
close all

%Here we read in multiple particle simulation runs from .csv files
SIRdir = dir('S_I_R_irate*.csv');
nfile = length(SIRdir);
gdir = dir('gi_function*.csv');

%Reads data from each file into a structure SIRdata, where .file contains
%the data in columns t, S, I, R, i and .g contains gi from generation
%intervals measured from the simulation. It is assumed that the number of 
%SIR and g_function files is the same.

for i=1:nfile   
    SIRdata(i).file=readmatrix(SIRdir(i).name,'NumHeaderLines',1);
    SIRdata(i).g = readmatrix(gdir(i).name,'NumHeaderLines',1);
end
%% Main loop to run models on particle simulation data

for q = 1:length(SIRdata)
    dt = SIRdata(q).file(2,1)-SIRdata(q).file(1,1);     %Find dt from data
    %Read data N_S(t), N_I(t), N_R(t), irate saved from previous particle model simulation
    t = SIRdata(q).file(:,1)';  %time (days)
    N_S = SIRdata(q).file(:,2)';  %Number of susceptible particles
    N_I = SIRdata(q).file(:,3)';  %Number of infected particles
    N_R = SIRdata(q).file(:,4)';  %Number of recovered particles
    i_rate = SIRdata(q).file(:,5)';
    
    N0 = 4000;  %Population of the simulation must be specified

        
    %****************Plot particle simulation data***********************
    figure(2);
    hold on
    plot(t,N_S/N0,'b.','MarkerSize',14) %S(t) 
    plot(t,N_I/N0,'r.','MarkerSize',14) %I(t)
    plot(t,N_R/N0,'g.','MarkerSize',14) %R(t)
    legend('S(t)','I(t)','R(t)')
    xlabel('Time (days)')
    ylabel('Population fraction')

    %Display the final infected fraction and compute R0 from final infected
    %fraction using Eq. 3.
    final_sick = 1-N_S(length(N_S))/N0;
    fprintf('Final fraction infected = %.3f \n',final_sick);
    R0_meas = -log(1-final_sick)/final_sick;
    fprintf('Measured R0 = %.3f \n',R0_meas);

    %Plot particle simulation infection rate
    figure(3);
    xlabel('Time (days)')
    ylabel('i(t)')
    hold on;
    plot(t,smooth(i_rate, 100)/dt,'r','LineWidth',2)     %infection rate

    %***********************end plot********************************

    %i_init is the early infection data used to find parameters and
    %initialize models
    t_i = 15;    %the duration of intial data (days)
    i_init = i_rate(1:t_i/dt+1)/dt;

   %Define infectious probability function by picking from following sections
%A, B, C, or D:
%% *A Fig. 2 and 3: Constant infectiousness gi(tau) = H(tau_m - tau)
% tau_m = 10; %Maximum infectious time (days)
% gi = ones(1,round(tau_m/dt)+1)/tau_m;   %infectious probability is constant
% P_I = ones(1,round(tau_m/dt)+1);        %Infectors are infected for entire time tau_m
%% *B Fig. 4 and 5: Linear gi(t) = 2/tau_m(1-tau/tau_m)H(tau_m - tau), used for random infectious period

% tau_m = 10; %Maximum infectious time (days)
% gi = zeros(1,tau_m/dt+1);
% P_I = zeros(1,tau_m/dt+1);
% for i=1:tau_m/dt+1
%     gi(i) = 2/tau_m*(1-(i-1)*dt/tau_m); %rand and linear, gi(tau) = 2/tau_m*(1-tau/tau_m)
%     P_I(i) = 1-((i-1)*dt/tau_m);  %random case, P_I(tau) = 1-tau/tau_m
% end
% P_I = ones(1,round(tau_m/dt)+1);    %Use for linear. Comment out for random

%% *C Exponential gi(tau). Does not appear in paper, but is equivalent to SIR results, as described in Discussion.
%Exponential gi(t) = gamma*exp(-gamma*t)
% tau_m = 20;   % the max recovery time should be set such that gi(tau) = gamma*exp(-gamma*tau_m) ~ 0.
% gi = zeros(1,tau_m/dt+1);
% P_I = zeros(1,tau_m/dt+1);
% i_init = ones(1,length(gi));
% i_init(1) = 0.00001;
% gamma = 1/3.0;    
% for i=1:length(gi)
%     i_init(i) = i_init(1)*exp((R0-1)*gamma*(i-1)*dt);
%     gi(i) = gamma*exp(-gamma*(i-1)*dt);
%     P_I(i) = exp(-gamma*(i-1)*dt);
% end
%% *D Use generation-interval distribution for gi(t) from simulation 'contact tracing'
%P_I(tau) must also be defined, if I(t) and R(t) are calculated. This
%depends on the simulation that was run.

    tau_m = 10; %Maximum infectious time (days)
    gi = zeros(1,tau_m/dt+1);

    %Find the best-fit gi(t) function from contact tracing
    %Linear gi(t) assumed
    coeff = polyfit([0:dt:(length(gi)-1)*dt], SIRdata(q).g/ sum(SIRdata(q).g*dt),1);
    for i=1:tau_m/dt+1
        gi(i) = coeff(1)*(i-1)*dt+coeff(2);
        if gi(i) < 0
            gi(i) = 0;
        end
    end

    gi_flip = flip(gi); %gi is applied backwards in time, so it should be flipped
    P_I = ones(1,round(tau_m/dt)+1);    %Use for linear or any situation where recovery time tau_r is the same for all particles.


%% *Contact tracing: use the generation-interval distribution and Eq. 1 to find R0_calc
%R0 can be calculated from knowledge of the generation-interval distribution
%and the early infection rate. However, especially early in an epidemic 
%when i(t) is small, there is a large uncertainty in that calculation, in 
%addition to uncertainty in the generation-interval distribution. Therefore,
%we apply Eq. 1 to find R0 or a range of times from tau_m (max length of
%the generation interval) to t_i, the initial data interval. 

R0_calc = zeros(1,round(t_i/dt));   %array to store values of R0 calculated from Eq.1. 
for i=round(tau_m/dt)+1:round(t_i/dt)    %Loop to find multiple values of R0_calc
    N_S_init = N0-sum(i_rate(1:i))*dt;   %Solve for N_S(t) = integral(i(t),{0,t}) from Eq.1
    R0_calc(i) = i_rate(i)/(trapz(i_rate(i-round(tau_m/dt):i).*flip(gi)*dt))/N_S_init*N0; %Solve for R0_calc from Eq. 1
end

R0 = sum(R0_calc)/((t_i-tau_m)/dt); %Average all calculated values of R0
fprintf('Calcuated R0 = %.3f \n',R0);

%% Main Renewal Model algorithm
%Here we use Eq. 1 and Eq. 2 to predict the epidemic after time ti.

% Generate the sequence of infections
N_S_RM = zeros(1,length(i_rate));    %The RM prediction of the number Susceptible, N_S
N_I_RM = zeros(1,length(i_rate)); %The RM prediction of the number Infected, N_I
N_R_RM = zeros(1,length(i_rate)); %The RM prediction of the number Recovered, N_R
i_RM = zeros(1,length(i_rate));  %The RM prediction of the infection rate i(t)
Reff_RM = zeros(1,length(i_rate));    %The RM prediction of the effective reproduction number Reff = R0*S(t)

n_init = length(i_init);    %number of values in initial infection rate i_init

%Set initial conditions at time t_i
i_RM(1:n_init) = i_init;    %Initial infection rate data spanning time 0 to ti is used to start RM model
I_init = sum(i_init)*dt;    %Initial Infected I_S(ti) = integral(i(t)*dt,{0,ti})
N_S_RM(n_init+1) = N0-I_init;   %Initial Susceptible N_S(ti) = N0-I(ti)
Reff_RM(n_init+1) = R0*N_S_RM(n_init+1)/N0; %Effective reproduction number Reff(ti) = R0*S(ti)

%Main loop to compute the Renewal model N_S(t), N_I(t), N_R(t) and i(t) using initial infection rate data i_init
for j=n_init+2:length(i_rate)
    i_RM(j) = Reff_RM(j-1)*trapz(i_RM(j-round(tau_m/dt)-1:j-1).*flip(gi))*dt; %Eq. 1, first equation
    N_S_RM(j) = N_S_RM(j-1)-i_RM(j)*dt; %Number of Susceptible. Eq. 1, second equation
    N_I_RM(j) = trapz(dt*i_RM(j-round(tau_m/dt)-1:j-1).*flip(P_I));   %Number of Infected. Eq. 2    
    N_R_RM(j) = N0-N_S_RM(j)-N_I_RM(j); %Recovered particles. Computed to conserve N0ulation
    Reff_RM(j) = R0*N_S_RM(j)/N0;   %Effective reproduction number Reff(t) = R0*S(t)
end

% Set up SIR to start the same as Renewal model for comparison
i_init_SIR = i_RM(length(i_init)+2);
N_S_init = N_S_RM(length(i_init)+2);
%% Plot RM model predictions
figure(2)
hold on

%Fill in initial data for t < t_i with particle simulation data
N_S_RM(1:round(t_i)/dt+2) = N_S(1:round(t_i)/dt+2);
N_I_RM(1:round(t_i)/dt+2) = N_I(1:round(t_i)/dt+2);
N_R_RM(1:round(t_i)/dt+2) = N_R(1:round(t_i)/dt+2);
i_RM(1:round(t_i)/dt+2) = smooth(i_rate(1:round(t_i)/dt+2),100)/dt;

%Plot Renewal model S, I, R population fractions. 
index_RM = round(t_i/dt)+3:length(t); %only plot array indices that come from RM prediction, not the initial data
%index_RM = 1:length(t);
plot(t(index_RM),N_S_RM(index_RM)/N0,'k--','LineWidth',1.5)
plot(t(index_RM),N_I_RM(index_RM)/N0,'k--','LineWidth',1.5)
plot(t(index_RM),N_R_RM(index_RM)/N0,'k--','LineWidth',1.5)

%Plot Renewal model infection rate i(t).
figure(3)
hold on
plot(t(index_RM),i_RM(index_RM),'k--','LineWidth',1.5)
xlabel('t (days)')
ylabel('infection rate (/day)')

%% SIR model 
% This script computes SIR model solutions, primarily as a comparison with
% the particle simulation and the Renewal model. SIR requires initial i(t)
% exponential growth, R0, and gamma, where gamma is interpreted as the 
% 'recovery rate'. This routine requires i_rate from the particle
% simulation for the initial exponential fit to i(t).

    %Exponential fit of seed data to find lambda, the exponential constant
    i_fit = i_rate/dt; %Convert i_rate in the infection rate for time step dt
    fit_coeff = fit(t(1:round(t_i/dt))',i_fit(1:round(t_i/dt))','exp1');    %Select on the initial data before t_i
    fit_num = coeffvalues(fit_coeff);
    lambda = fit_num(2);    %The exponential fit constant from i(t) = i(0)*exp(lambda*t)
    gamma = lambda/(R0-1);  %Use the non-parametric approach where R0 is determined from initial i(t) data using RM model
    
    %Find time shift to plot all data with similar initial growth
    ishift = 40;
    tshift = log(ishift/fit_num(1))/lambda;
    SIRdata(q).tshift = floor(tshift/dt)*dt;    %store tshift in increments of dt

%% Main SIR algorithm
%Initial conditions for SIR are set to match the RM model infection rate at
%t = ti for comparison purposes. See calculation for N_S_init and i_init
%above. An alternate is to use the infection rate at t = ti determined
%from the exponential fit to initial infection rate data.

beta = R0*gamma;    %SIR transmission coefficient
I_init = 1/(gamma*R0*N_S_init/N0)*i_init_SIR;   %Determine initial I(t) from the initial i(t) according to Eq. 4
R_init = (N0-I_init-N_S_init);  %Determine initial R(t) from population conservation

%Initialize SIR arrays
S_SIR = ones(1,length(i_rate));
I_SIR = zeros(1,length(i_rate));
R_SIR = zeros(1,length(i_rate));
i_SIR = zeros(1,length(i_rate));

%Set initial t = ti population fractions
S_SIR(length(i_init)+2)= N_S_init/N0;
I_SIR(length(i_init)+2)= I_init/N0;
R_SIR(length(i_init)+2)= R_init/N0;

%Main loop to compute SIR model
for i=length(i_init)+3:length(i_rate) 
        S_SIR(i) = S_SIR(i-1) - beta*S_SIR(i-1)*I_SIR(i-1)*dt;  %Susceptible population fraction. Eq. 4, first equation
        I_SIR(i) = I_SIR(i-1) + beta*S_SIR(i-1)*I_SIR(i-1)*dt - gamma*I_SIR(i-1)*dt; % Infected population fraction. Eq. 4, second equation
        R_SIR(i) = R_SIR(i-1) + gamma*I_SIR(i-1)*dt;  % Recovered population fraction.
        i_SIR(i) = beta.*S_SIR(i).*I_SIR(i)*N0;
end


%% Plot SIR model predictions
figure(2)
hold on
%Set initial values of SIR (before prediction) to be the same as RM for
%comparison
    S_SIR(1:t_i/dt+2) = N_S_RM(1:t_i/dt+2)/N0;
    I_SIR(1:t_i/dt+2) = N_I_RM(1:t_i/dt+2)/N0;
    R_SIR(1:t_i/dt+2) = N_R_RM(1:t_i/dt+2)/N0;
    i_SIR(1:t_i/dt+2) = i_RM(1:t_i/dt+2);

    index_RM = round(t_i/dt)+3:length(t); %only plot array indices that come from SIR prediction, not the initial data

    plot(t(index_RM),S_SIR(index_RM),':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
    plot(t(index_RM),I_SIR(index_RM),':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
    plot(t(index_RM),R_SIR(index_RM),':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
    legend('S','I','R')

    %Plot infection rate, comparable to i(t)
    figure(3)
    hold on
    plot(t(index_RM),beta.*S_SIR(index_RM).*I_SIR(index_RM)*N0,':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
    
    %Save General Infectiousness and SIR SIR model data
    SIRdata(q).RM = [t; N_S_RM; N_I_RM; N_R_RM; i_RM];
    SIRdata(q).SIR = [t; S_SIR; I_SIR; R_SIR; i_SIR];
end     %Main loop cycling through files

%% Process particle simulation and analytic model data

%Find the minimum time value where i(t) = ishift occurs.
tshiftmin = floor(SIRdata(1).tshift/dt)*dt;

for q=2:length(SIRdata)
    tshiftmin = min(floor(SIRdata(q).tshift/dt)*dt,tshiftmin);
end

tshiftmin = floor(tshiftmin/dt)*dt;

%Find the minimum length of the shifted data
lengthmin = length(SIRdata(1).RM(1,:))-(SIRdata(1).tshift-tshiftmin)/dt;

for q=2:length(SIRdata)
    lengthq = length(SIRdata(q).RM(1,:))-(SIRdata(q).tshift-tshiftmin)/dt;
    lengthmin = min(lengthq,lengthmin);
end
lengthmin = floor(lengthmin);

%Initialize arrays containing timeshifted, truncated data
N_S_part = zeros(length(SIRdata),lengthmin);
N_I_part = zeros(length(SIRdata),lengthmin);
N_R_part = zeros(length(SIRdata),lengthmin);
i_part = zeros(length(SIRdata),lengthmin);

N_S_RM = zeros(length(SIRdata),lengthmin);
N_I_RM = zeros(length(SIRdata),lengthmin);
N_R_RM = zeros(length(SIRdata),lengthmin);
i_RM = zeros(length(SIRdata),lengthmin);

S_SIR = zeros(length(SIRdata),lengthmin);
I_SIR = zeros(length(SIRdata),lengthmin);
R_SIR = zeros(length(SIRdata),lengthmin);
i_SIR = zeros(length(SIRdata),lengthmin);

t = 0:dt:(lengthmin-1)*dt;

gi = zeros(length(SIRdata),length(SIRdata(1).g));
t_gi = 0:dt:(length(gi(1,:))-1)*dt;

for q=1:length(SIRdata)
    
    tshift = floor((SIRdata(q).tshift-tshiftmin)/dt)*dt;   %Find how much to shift data
    
    N_S_part(q,:) = SIRdata(q).file(tshift/dt+1:tshift/dt+lengthmin,2)';
    N_I_part(q,:) = SIRdata(q).file(tshift/dt+1:tshift/dt+lengthmin,3)';
    N_R_part(q,:) = SIRdata(q).file(tshift/dt+1:tshift/dt+lengthmin,4)';
    i_part(q,:) = smooth(SIRdata(q).file(tshift/dt+1:tshift/dt+lengthmin,5),100)';

    N_S_RM(q,:) = SIRdata(q).RM(2,tshift/dt+1:tshift/dt+lengthmin);
    N_I_RM(q,:) = SIRdata(q).RM(3,tshift/dt+1:tshift/dt+lengthmin);
    N_R_RM(q,:) = SIRdata(q).RM(4,tshift/dt+1:tshift/dt+lengthmin);
    i_RM(q,:) = smooth(SIRdata(q).RM(5,tshift/dt+1:tshift/dt+lengthmin),100);

    S_SIR(q,:) = SIRdata(q).SIR(2,tshift/dt+1:tshift/dt+lengthmin);
    I_SIR(q,:) = SIRdata(q).SIR(3,tshift/dt+1:tshift/dt+lengthmin);
    R_SIR(q,:) = SIRdata(q).SIR(4,tshift/dt+1:tshift/dt+lengthmin);
    i_SIR(q,:) = smooth(SIRdata(q).SIR(5,tshift/dt+1:tshift/dt+lengthmin),100);

    gi(q,:) = smooth(SIRdata(q).g,100)/sum(SIRdata(q).g*dt);

    figure(2)
    hold on
    plot(t,N_S_part(q,:)/N0,'b.','MarkerSize',14) 
    plot(t,N_I_part(q,:)/N0,'r.','MarkerSize',14);
    plot(t,N_R_part(q,:)/N0,'g.','MarkerSize',14)

    xlabel('t (days)')
    ylabel('Population fraction')

    
    plot(t,S_SIR(q,:),'k--','LineWidth',1.5)
    plot(t,I_SIR(q,:),'k--','LineWidth',1.5)
    plot(t,R_SIR(q,:),'k--','LineWidth',1.5)


    figure(3)
    hold on

    plot(t,smooth(i_part(q,:),100)/dt,'r','LineWidth',2)
    plot(t,smooth(i_RM(q,:),100),'k--','LineWidth',1.5)
    xlabel('t (days)')
    ylabel('infection rate (/day)')

end

figure(2)
legend('S(t)','I(t)','R(t)')

%Compute mean and standard deviation
N_S_part_mean = mean(N_S_part);
N_I_part_mean = mean(N_I_part);
N_R_part_mean = mean(N_R_part);
i_part_mean = mean(i_part);


N_S_part_std = std(N_S_part);
N_I_part_std = std(N_I_part);
N_R_part_std = std(N_R_part);
i_part_std = std(i_part);

N_S_RM_mean = mean(N_S_RM);
N_I_RM_mean = mean(N_I_RM);
N_R_RM_mean = mean(N_R_RM);
i_RM_mean = mean(i_RM);

N_S_RM_std = std(N_S_RM);
N_I_RM_std = std(N_I_RM);
N_R_RM_std = std(N_R_RM);
i_RM_std = std(i_RM);

S_SIR_mean = mean(S_SIR);
I_SIR_mean = mean(I_SIR);
R_SIR_mean = mean(R_SIR);
i_SIR_mean = mean(i_SIR);

S_SIR_std = std(S_SIR);
I_SIR_std = std(I_SIR);
R_SIR_std = std(R_SIR);
i_SIR_std = std(i_SIR);

gi_mean = mean(gi);
gi_std = std(gi);

%Generate fill area +- standard deviation for RM
Scurve1 = (N_S_RM_mean+N_S_RM_std)/N0;
Scurve2 = (N_S_RM_mean-N_S_RM_std)/N0;
tfill = [t,fliplr(t)];
SinBetween = [Scurve1,fliplr(Scurve2)];

Icurve1 = (N_I_RM_mean+N_I_RM_std)/N0;
Icurve2 = (N_I_RM_mean-N_I_RM_std)/N0;
IinBetween = [Icurve1,fliplr(Icurve2)];

Rcurve1 = (N_R_RM_mean+N_R_RM_std)/N0;
Rcurve2 = (N_R_RM_mean-N_R_RM_std)/N0;
RinBetween = [Rcurve1,fliplr(Rcurve2)];

icurve1 = (i_RM_mean+i_RM_std);
icurve2 = (i_RM_mean-i_RM_std);
iinBetween = [icurve1,fliplr(icurve2)];

%Generate fill area +- standard deviation for SIR
SSIRcurve1 = (S_SIR_mean+S_SIR_std);
SSIRcurve2 = (S_SIR_mean-S_SIR_std);
tfill = [t,fliplr(t)];
SSIRinBetween = [SSIRcurve1,fliplr(SSIRcurve2)];

ISIRcurve1 = (I_SIR_mean+I_SIR_std);
ISIRcurve2 = (I_SIR_mean-I_SIR_std);
ISIRinBetween = [ISIRcurve1,fliplr(ISIRcurve2)];

RSIRcurve1 = (R_SIR_mean+R_SIR_std);
RSIRcurve2 = (R_SIR_mean-R_SIR_std);
RSIRinBetween = [RSIRcurve1,fliplr(RSIRcurve2)];

iSIRcurve1 = (i_SIR_mean+i_SIR_std);
iSIRcurve2 = (i_SIR_mean-i_SIR_std);
iSIRinBetween = [iSIRcurve1,fliplr(iSIRcurve2)];

gicurve1 = (gi_mean+gi_std);
gicurve2 = (gi_mean-gi_std);
giinBetween = [gicurve1,fliplr(gicurve2)];
tgifill = [t_gi,fliplr(t_gi)];

%Plot mean particle simulation population fractions and mean SIR/RM models +- stand dev.
figure(4)
hold on
fill(tfill,SinBetween,[.2 .2 1],'EdgeColor',[.2 .2 1])
fill(tfill,IinBetween,'r','EdgeColor','r')
fill(tfill,RinBetween,'g','EdgeColor','g')

fill(tfill,SSIRinBetween,[.5 .5 .9],'EdgeColor',[.5 .5 .9])
fill(tfill,SinBetween,[.2 .2 1],'EdgeColor',[.2 .2 1])
plot(t,N_S_part_mean/N0,'k','LineWidth',1.5) 

fill(tfill,ISIRinBetween,[.9 .5 .5],'EdgeColor',[.9 .5 .5])
fill(tfill,IinBetween,'r','EdgeColor','r')
plot(t,N_I_part_mean/N0,'k','LineWidth',1.5)

fill(tfill,RSIRinBetween,[.5 .9 .5],'EdgeColor',[.5 .9 .5])
fill(tfill,RinBetween,'g','EdgeColor','g')
plot(t,N_R_part_mean/N0,'k','LineWidth',1.5)
xlabel('t (days)')  
ylabel('S,I,R (population fraction)')
legend('S','I','R')

ylim([0 1])

figure(5)
hold on
%Plot the infection rate +- the standard devation for SIR and RM models
fill(tfill,iSIRinBetween,[.5 .5 .9],'EdgeColor',[.5 .5 .9])
fill(tfill,iinBetween,[.2 .2 1],'EdgeColor',[.2 .2 1])
plot(t,i_part_mean/dt,'r','LineWidth',1.5) 
xlabel('t (days)')
ylabel('infection rate (/day)')
ylim([0 250])
legend('SIR','RM','Particle')

figure(6)
hold on
fill(tgifill,giinBetween,[.5 .5 .9],'EdgeColor',[.5 .5 .9])
plot(t_gi,-0.02*t_gi+0.2,'r','LineWidth',1.5) %Plot the function p_i(tau) used in the particle simulation
xlabel('\tau (days)')
ylabel('serial interval (/day)')
ylim([0 .3])
legend('serial interval','g_i')
