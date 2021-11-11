%% RM and SIR model code for "Active particle infection dynamics as a test of analytic epidemic models"
%How to use this script: When running RM and SIR models, choices must be
%made regarding the amount of initial particle code data used, whether or
%not to use 'contact tracing', and choice of parameters. Please go through
%each section to determine whether or not to run it. Sections that require
%a choice are marked with a "*" in their section title.

%RM and SIR models use "initial data" from particle code to determine
%parameters and predict the rest of the simulated epidemic. If using
%"i_rate" data, you must run the particle simulation to generate i_rate
%data first. Otherwise, i_rate can be read from a file
N0 = 4000;  %Use population from simulation
R0 = 2;     %If R0 is assumed known from particel simulation, input here. Otherwise compute using R0_calc below
dt = 0.02;  %Timestep used in particle simulation

%% *Read particle simulation data from a file.
%If reading simulation data from file, run the following section, otherwise
%skip this section.

dt = m(1,2)-m(1,1);     %Find dt from data

%Read data S(t), I(t), R(t) saved from previous particle model simulation
m = readmatrix('S_I_R_irate1.csv');
t = m(:,1)';
N_S = m(:,2)';
N_I = m(:,3)';
N_R = m(:,4)';
i_rate = m(:,5)'/dt;

%Plot simulation from file. Commment out following lines if using
%simulation
figure(2)
hold on
plot(t,N_S/N0,'b.','MarkerSize',14)
plot(t,N_I/N0,'r.','MarkerSize',14)
plot(t,N_R/N0,'g.','MarkerSize',14)
legend('S','I','R')
xlabel('t (days)')  
ylabel('S,I,R (population fraction)')

figure(3)
hold on
plot(t,smooth(i_rate, 100),'r','LineWidth',2)    %Plot smoothed infection rate per day

%% Pick infectious probability function gi(tau) and the fraction of infectors function P_I(tau)

%Select the initial data from particle code to run RM and SIR models.
t_i = 15;    %the duration of intial data (days) 
i_init = i_rate(1:t_i/dt+1);  %Select infection rate data from R0_calc

%Define infectious probability function by picking from following sections
%A, B, C, or D:
%% *A Fig. 2 and 3: Constant infectiousness gi(tau) = H(tau_m - tau)
tau_m = 10; %Maximum infectious time (days)
gi = ones(1,round(tau_m/dt)+1)/tau_m;   %infectious probability is constant
P_I = ones(1,round(tau_m/dt)+1);        %Infectors are infected for entire time tau_m
%% *B Fig. 4 and 5: Linear gi(t) = 2/tau_m(1-tau/tau_m)H(tau_m - tau), used for random infectious period

tau_m = 10; %Maximum infectious time (days)
gi = zeros(1,tau_m/dt+1);
P_I = zeros(1,tau_m/dt+1);
for i=1:tau_m/dt+1
    gi(i) = 2/tau_m*(1-(i-1)*dt/tau_m); %rand and linear, gi(tau) = 2/tau_m*(1-tau/tau_m)
    P_I(i) = 1-((i-1)*dt/tau_m);  %random case, P_I(tau) = 1-tau/tau_m
end
P_I = ones(1,round(tau_m/dt)+1);    %Use for linear. Comment out for random

%% *C Exponential gi(tau). Does not appear in paper, but is equivalent to SIR results, as described in Discussion.
%Exponential gi(t) = gamma*exp(-gamma*t)
tau_m = 20;   % the max recovery time should be set such that gi(tau) = gamma*exp(-gamma*tau_m) ~ 0.
gi = zeros(1,tau_m/dt+1);
P_I = zeros(1,tau_m/dt+1);
i_init = ones(1,length(gi));
i_init(1) = 0.00001;
gamma = 1/3.0;    
for i=1:length(gi)
    i_init(i) = i_init(1)*exp((R0-1)*gamma*(i-1)*dt);
    gi(i) = gamma*exp(-gamma*(i-1)*dt);
    P_I(i) = exp(-gamma*(i-1)*dt);
end
%% *D Use generation-interval distribution for gi(t) from simulation 'contact tracing'
gi = gi_fit;    %Use best linear fit to particle simulation generation-interval distribution

%P_I(tau) must also be defined, if I(t) and R(t) are calculated. This
%depends on the simulation that was run.

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

t = linspace(0,(length(i_rate)-1)*dt,length(i_rate)-1);  % build time array

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
    N_R_RM(j) = N0-N_S_RM(j)-N_I_RM(j); %Recovered particles. Computed to conserve population
    Reff_RM(j) = R0*N_S_RM(j)/N0;   %Effective reproduction number Reff(t) = R0*S(t)
end

% Set up SIR to start the same as Renewal model
i_init_SIR = i_RM(length(i_init)+2);
N_S_init = N_S_RM(length(i_init)+2);
%% Plot RM model predictions
figure(2)
hold on
%Plot Renewal model S, I, R population fractions. 
index_RM = round(t_i/dt)+3:length(t); %only plot array indices that come from RM prediction, not the initial data
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

t = linspace(0,(length(i_rate)-1)*dt,length(i_rate)-1);  % build time array if it does not already exist

%Exponential fit of seed data to find lambda, the exponential constant
i_fit = i_rate/dt; %Convert i_rate in the infection rate for time step dt
fit_coeff = fit(t(1:round(t_i/dt))',i_fit(1:round(t_i/dt))','exp1');    %Select on the initial data before t_i
fit_num = coeffvalues(fit_coeff);
lambda = fit_num(2);    %The exponential fit constant from i(t) = i(0)*exp(lambda*t)

%% *Pick SIR parameters from following choices
%The SIR model uses either the A) parametric or the B) non-parametric approach to
%find parameters R0 and gamma to run the model, along with initial i(t)
%data.

%% *A) Parametric approach (Ma 2020)  
%The parametric approach has the modeler choose a value of gamma, the
%recovery rate from knowledge of the generation-interval distribution.
%Common choices are gamma = 1/<tau_r> (using the average generation
%interval, or gamma = 1/tau_m (using the maximum recovery time.

gamma = 1/10;            %Guess value of gamma
R0 = lambda/gamma+1;     %Solve for R0 based on gamma
%% *B) Non-parametric approach (Ma 2020) 
%The non-parametric approach has the modeler find R0 using the Renewal
%model equation (Eq. 1) from initial i(t) data. Run "Contact tracing"
%routine above to calculate R0.  Or simply assume R0 from a particle
%simulation.
gamma = lambda/(R0-1);

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

%Set initial t = ti population fractions
S_SIR(length(i_init)+2)= N_S_init/N0;
I_SIR(length(i_init)+2)= I_init/N0;
R_SIR(length(i_init)+2)= R_init/N0;

%Main loop to compute SIR model
for i=length(i_init)+3:length(i_rate) 
        S_SIR(i) = S_SIR(i-1) - beta*S_SIR(i-1)*I_SIR(i-1)*dt;  %Susceptible population fraction. Eq. 4, first equation
        I_SIR(i) = I_SIR(i-1) + beta*S_SIR(i-1)*I_SIR(i-1)*dt - gamma*I_SIR(i-1)*dt; % Infected population fraction. Eq. 4, second equation
        R_SIR(i) = R_SIR(i-1) + gamma*I_SIR(i-1)*dt;  % Recovered population fraction.
end

%% Plot SIR model predictions
figure(2)
hold on
index_RM = round(t_i/dt)+3:length(t); %only plot array indices that come from SIR prediction, not the initial data
plot(t(index_RM),S_SIR(index_RM),':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
plot(t(index_RM),I_SIR(index_RM),':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
plot(t(index_RM),R_SIR(index_RM),':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
legend('S','I','R')

%Plot infection rate, comparable to i(t)
figure(3)
hold on
plot(t(index_RM),beta.*S_SIR(index_RM).*I_SIR(index_RM)*N0,':','LineWidth',1.5,'Color',[0.30,0.75,0.93])
