%% Active particle code for "Active particle infection dynamics as a test of analytic epidemic models"
% This code sets up a square region with periodic boundary conditions that
% contains particles in the Susceptible, Infected, and Recovered classes. A
% selected number of particles start in the Infected class and move in
% randomly chosen directions at fixed velocity through a sea of stationary
% Susceptible class particles. When an Infected class particle comes within
% one infection radius of a Susceptible particle a random number p1 is
% generated and if p1 < pi, the probability of infection per encounter,
% then an infection is passed and the Susceptible particle becomes Infected
% and starts moving in its own random direction. The probability of
% infection per encounter pi can be a function of the time after infection
% tau or a function of the particle itself. This dependence is determined
% in the initialization of the system. 
% The simulation may be displayed as an animation by un-commenting line
% xxx, but this slows down the code runtime. Particles are displayed as
% blue Susceptible, red Infected, and green Recovered circles. At the end
% of the simulation, the Susceptible, Infected, and Recovered fractions of
% the population are plotted, along with the smoothed infection rate
% (infections/day). 
% The simulation records the time between when a particle is infected and
% when it passes on an infection to measure the generation interval
% distribution of infections. This distribution is recorded only for the
% initial infection data and is comparable to serial intervals recorded by
% contact tracing a real epidemic. The generation interval distribution is
% plotted at the end of the simulation, gi(tau). Because of the small
% number of recorded infections, the generation interval is subject to
% uncertainty and will not represent the true generation interval
% distribution.

%Clear all variables and figures
clc
close all
clear all

%Main simulation variables (see Table 1)
N0 = 1000; %Number of particles
L = 8; %Width of the simulation (arbitrary length unit)
v0 = 3.0;   %Infected particle velocity (arbitrary length unit/day)
duration = 10000;   %Maximum duration of a simulation (days)
dt = .02;   %Time step (days)
nsteps = round(duration / dt) + 1; % number of time steps in simulation
t = 0; %Time elapsed in simulation (days)
ri = 0.1;   %Radius of infection (arbitrary length unit)
R0_expect = 2; %Expected basic reproduction number


% Choose from the following options the probability of infection per encounter and the infectious time distribution
%Infected particles may have a time-dependent probability of passing an
%infection pi and a distribution in infectious times tau_r in our
%simulation. Figures in our paper were generated with different assumptions
%that can be chosen below.

%% Fig. 2 and 3: Constant infectiousness. All particles recover at time tau_r = tau_m (days)
tau_m = 10; %Recovery time for all particles (days)
tau_r = ones(1,N0)*tau_m;   %Array to store recovery time for each particle
tau_ave = tau_m;    %Average recovery time is the same as tau_m

Nc = 2*ri*v0*tau_ave*N0/L^2;    %Average number of close encounters per infector
pi = R0_expect/Nc*ones(1,tau_m/dt+1); %Probability per encounter of passing infection as a function of tau time after infection


%% Fig. 4 and 5: "Linear" pi(tau), probability of infectiousness drops linearly with time after infection
tau_m = 10; %Recovery time for all particles (days)
tau_r = ones(1,N0)*tau_m;   %Array to store recovery time for each particle
tau_ave = 10;   %Average recovery time is the same as tau_m

Nc = 2*ri*v0*tau_ave*N0/L^2;    %Average number of close encounters per infector
pi = zeros(1,tau_m/dt+1);   %Probability per encounter of passing infection as a function of time after infection
p = R0_expect/Nc; 

for i=1:length(pi)
    tau = (i-1)*dt;
    pi(i) = 2*p*(1-tau/tau_m);
end

%% Fig. 4: Random distribution of recovery times, each particle with constant infectiousness
tau_m = 10;
tau_r = rand(1,N0)*tau_m;   %Generate random recovery times from 0 to tau_m for each particle
tau_ave = tau_m/2;   %Average recovery time is 1/2 tau_m

Nc = 2*ri*v0*tau_ave*N0/L^2;    %Average number of close encounters per infector
pi = R0_expect/Nc*ones(1,tau_m/dt+1);   %Probability per encounter of passing infection as a function of time after infection

%% Not found in paper: COVID serial interval distribution from Nishiura et al., assuming each particles' pi(tau) has the same shape as the serial interval distribution 
tau_m = 20; %Recovery time for all particles (days)
tau_r = ones(1,N0)*tau_m;   %Array to store recovery time for each particle
tau_ave = 20;   %Average recovery time is the same as tau_m

%Generate lognormal distribution following Nishiura et al.
mean_si = 4.7;
std_si = 2.9;
mu_si = log(mean_si^2/sqrt(std_si^2+mean_si^2));
sigma_si = sqrt(log(std_si^2/mean_si^2 + 1));
pd = makedist('Lognormal','mu',mu_si,'sigma',sigma_si);
Nc = 2*ri*v0*tau_ave*N0/L^2;    %Average number of close encounters per infector
gi_t = pdf(pd,gen_time);
pi = R0_expect/Nc*pdf(pd,0:dt:tau_m);

%% Main simulation section
%Here we set up the simulation to run until no infected particles remain.

%Eq. 5: Solve for final fraction Recovered using the expected R0
fun1 = @(x)log(1-x)+R0_expect*x;   %Generate the transcendental equation 5
R_final = fzero(fun1,[0.0001 0.9999]);   %Solve Eq. 5 for 1-S_final and print to Command Window
fprintf('Prediction of final fraction infected = %.3f \n',R_final);

%Generation interval distribution function, gi(tau), similar to contact
%tracing
gi = zeros(1,length(pi));   %Initialize generation interval distribution function
ti = 15;    %Record generation intervals up to ti (days)

%Record the time each particle becomes infected. This is used for
%calculating generation intervals, generating pi, and finding recovery
%time.
t_infected = zeros(1,N0);%the t value that a particle was infected on


%The following series of arrays temporarily stores the states of each
%particle as the simulation advances.
Infected = zeros(1,N0); %The particles currently infected
Recovered = zeros(1,N0); %The particles currently recovered

Color = zeros(N0,3); % Array representing RGB colors of each particle in the simulation
Color(:,3) = 1; %Initialize all particles as blue Susceptible particles

I_init = 5; %Number of initial infected particles
patientzero = randi(N0,[1 I_init]); %Chooses random indices as first infected particles
Infected(patientzero) = 1;  %Sets these particles as currently infected
t_infected(patientzero) = dt/10;   %Sets the time these particles become infected as right after t = 0;
Color(patientzero,1) = 1;   %Infected particles are red
Color(patientzero,3) = 0;   %not blue

%Records particles who have ever been infected (currently or previously)
Infected_prev = Infected; 

%The following arrays store particle numbers at each time step dt.
N_S = zeros(1,nsteps);  %Susceptible particles as a function of time in the simulation
N_I = zeros(1,nsteps); %Infected particles as a function of time in the simulation
N_R = zeros(1,nsteps);  %Recovered particles as a function of time in the simulation
%time = linspace(0,duration,round(nsteps)); %array for keeping track
%of timesteps for final plots
N_I_total = zeros(1,nsteps);%Total number of infections as a function of time

prox_detect = zeros(N0,N0); %Records which two particles are with one infection radius ri of each other

%Fern edit
%i_rateDaily = zeros(1,duration); %Infection rate i(t) in infections per day
i_rate = zeros(1,nsteps); %Infection rate i(t) (/day)

%diameter = 0.02; % diameter = particle size     
%r = diameter/2;
x = zeros(1,N0);
y = zeros(1,N0);
v_x = zeros(1,N0);
v_y = zeros(1,N0);
%a = zeros(1,N0);
%b = zeros(1,N0);

for z = 1:N0         % Initialising positions and velocities of each particle
    x(z) = rand()*L;
    y(z) = rand()*L;
    v_x(z) = rand()*2*v0-v0;
    %v_y(z) = rand()*2*v0-v0;
    v_y(z) = (randi([0 1])*2-1)*sqrt(v0^2 - v_x(z)^2);
end


%% Plot particle positions
%This section should be run if plotting particle positions is desired,
%otherwise it should be commented out.  Plotting particles significantly
%slows running of the code.

res = get(0,'screensize');
figure('rend','painters','pos',[0.2*res(3) 0.2*res(4) 0.6*res(3) 0.5*res(4)])
fig_res = get(gcf, 'Position');

s = scatter(x,y,10,Color,'filled'); %Plots a scatter. Positions are updated through handle.
hold on
plot([L,L],[0,L],'color',[0 0 0]);
plot([0,L],[L,L],'color',[0 0 0]);
plot([0,L],[0,0],'color',[0 0 0]);
plot([0,0],[0,L],'color',[0 0 0]);
axis([0 L 0 L]);
axis square
%% Main particle code loop

dt_count = 1; %Initialize time counter
t = dt; 
while t < duration && sum(Infected) > 0 %Continue to run the loop until there are no Infected particles 
    for i = 1:N0   %Loop to move infected particles
         if Infected(i) == 1
           x(i) = x(i) + dt*v_x(i); % Updates the x coordinates of each particle based on velocity
           y(i) = y(i)+dt*v_y(i); % Updates the y coordinates of each particle based on velocity 
        end
        if x(i) >= L    % Periodic boundary condition
           x(i) = L - x(i);
        elseif x(i) <= 0
           x(i) = L - x(i);
        end      
        if y(i) >= L    % Periodic boundary conditio
           y(i) = L - y(i);
        elseif y(i) <= 0     
           y(i) = L - y(i);
        end
    end   %loop moving infected particles
    
    for j=1:N0      %Loop through particles to find proximity within ri
        for k=j+1:N0    %of other particles
            if (x(j)-x(k))^2 + (y(j)-y(k))^2 < ri^2 % Find if two particles are within ri distance  
                if Recovered(j) == 0 && Recovered(k) ==0 && Infected(j) ~= Infected(k) && prox_detect(j,k) == 0 %if particles are susceptible and infected and have not previously encountered each other
                    if Infected(j) ==1
                        Prob = pi(floor((t-t_infected(j))/dt)+1);   %Find pi(tau) for the particle
                    else
                        Prob = pi(floor((t-t_infected(k))/dt)+1);   %Find pi(tau) for the particle
                    end
                    if rand()<=Prob     %if p1 <= pi an infection will be passed
                        Infected_prev(j) = 1;
                        Infected_prev(k) = 1;
                        i_rate(dt_count) = i_rate(dt_count)+1;    % record the infection in i(t)
                    
                        if t_infected(j) == 0   %Find whether particle has already been infected
                           t_infected(j) = t;   %if not, then record the time of infection

                           %Generation interval record
                           if t_infected(k) < ti    %Record generation interval for particle k if earlier than ti
                              gi(floor((t-t_infected(k))/dt)+1) = gi(floor((t-t_infected(k))/dt)+1) + 1;
                           end 
                        else
                           t_infected(k) = t;   %Record time of infection if particle k just became infected
                            
                           if t_infected(j) < ti    %Record generation interval for particle j if earlier than ti
                               gi(floor((t-t_infected(j))/dt)+1) = gi(floor((t-t_infected(j))/dt)+1) + 1;
                           end
                        end
                        Color(j,3) = 0; %Set both j and k particles infected
                        Color(j,1) = 1;
                        Color(k,3) = 0;
                        Color(k,1) = 1;
                end
                
            prox_detect(j,k)=1; %records proximity from this encounter
            prox_detect(k,j)=1;
                end    
            else prox_detect(j,k)=0;    %or not proximate
                prox_detect(k,j) = 0;
            end

        end     %kth particle loop
    end    %jth particle loop
    Infected = Infected_prev;   %
    for q=1:N0
        if Infected(q) == 1 && t - t_infected(q) >= tau_r(q) %if a particle is infected and is past their recovery time, it recovers
            Infected(q) = 0;    %Particle is no longer infected
            Recovered(q) = 1;   %It is recovered
            Color(q,1) = 0;     %Set colors for plotting particles
            Color(q,2) = 1;
            Color(q,3) = 0;
        end
    end
    
    t = t + dt;
    
    % Fern edit
    %dt_count = round(t/dt);
    dt_count = dt_count + 1;
    N_I(dt_count) = sum(Infected);
    N_I_total(dt_count) = sum(Infected_prev);
    N_R(dt_count) = sum(Recovered);
    
 %*************Uncomment if plotting positions. Otherwise, 
    s.XData = x; %handle for datapoints
    s.YData = y;
    pause(0.001)
    title(sprintf('Active particle simulator - t = %.3f',t));
    set(s,'cdata',Color);
%*******************
    
end     %Main simulation loop

t = linspace(0,t,dt_count); %Generate time array
i_rate = i_rate(1:dt_count)/dt; %measure i(t) in infections/day
N_S = N0 - N_I_total(1:dt_count);   %Calculate Susceptible from those who have not been infected
N_I = N_I(1:dt_count);  %Infected particles vs. time
N_R = N_R(1:dt_count);  %Recovered particles vs. time


%% Plot S, I, R vs. t
figure(2);
hold on

%Plot S, I, R vs. time
plot(t,N_S/N0,'b.','MarkerSize',14) 
plot(t,N_I/N0,'r.','MarkerSize',14);
plot(t,N_R/N0,'g.','MarkerSize',14)
legend('S','I','R')
xlabel('t (days)')
ylabel('S,I,R (population fraction)')

final_sick = 1-N_S(length(N_S))/N0;  %Find the final sick fraction
fprintf('Measured final fraction infected = %.3f \n',final_sick);

R0_meas = -log(1-final_sick)/final_sick;    %The effective R0 for the simulation can be computed using Eq. 5
fprintf('Measured effective R0 = %.3f \n',R0_meas);

%% Infection Rate Plot
% The array i_rate contains the number of infections at each time step. To
% convert this into the infection rate, we dived by the time step dt.
figure(3);
xlabel('t (days)')
ylabel('infection rate (/day)')
hold on;

plot(t,smooth(i_rate, 100),'r','LineWidth',2)    %Plot smoothed infection rate per day
%% Generation interval distribution gi(t) plot. Comparable to contact tracing

figure(5)
hold on
tau = [0:dt:(length(pi)-1)*dt]; %Time tau after infection
plot(tau, smooth(gi,100) / sum(gi*dt), 'r-','LineWidth',1.5)   %Plot smoothed generation interval distribution
plot(tau,pi/(sum(pi)*dt),'b-','LineWidth',1.5) %Plot normalized pi(tau)

%Find best fit line to the generation-interval distribution.  Use only for linear and constant gi(t)
coeff = polyfit(tau, gi/sum(gi*dt),1); %Linear fit coefficients to generation-interval distribution gi(tau)
gi_fit = tau*coeff(1)+coeff(2); %Define the best-fit linear gi(tau)
gi_fit(gi_fit<0)=0; %probability of infectiousness cannot be zero
plot(tau,gi_fit,'k-','LineWidth',1.5)

legend('g_i(\tau)','p_i(\tau)','linear fit g_i(\tau)')
xlabel('time after infection (days)')
ylabel('probability of infectiousness')


%% Record simulation data to a spreadsheet

headers = ["time","N_S", "    N_I    ", "    N_R    ", "i(t)"];
Data = [t' N_S' N_I' N_R' i_rate'/dt];
writematrix(headers,'S_I_R_irate1.csv')
writematrix(Data, 'S_I_R_irate1.csv','WriteMode','append')

header2 = ["gi(t) function"];
writematrix(header2, 'gi1.csv')
writematrix(gi', 'gi1.csv','WriteMode','append')
