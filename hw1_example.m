%% Example code for BGGN 260 Homework 1, Fall 2016
% Questions to anp054@eng.ucsd.edu, juw087@eng.ucsd.edu, or gert@ucsd.edu

%% Setup constants
V_t = 26.7; % thermal voltage, in mV
C_m = 1;    % membrane capacitance, in uF/cm^2

% Free ion concentrations of Na+, K+, and Cl- 
C_Na_o = 145;  % in mM
C_Na_i = 10;
C_K_o  = 5;
C_K_i  = 140;
C_Cl_o = 110;
C_Cl_i = 4;

% Conductances
g_Na = 12; % in mS/cm^2
g_K  = 4;
g_Cl = 0.1;

% Nernst reversal potentials
E_Na = V_t * log(C_Na_o / C_Na_i); % in mV
E_K  = V_t * log(C_K_o  / C_K_i );
E_Cl = V_t * log(C_Cl_i / C_Cl_o); % reversed since Cl is negatively charged

% Permeabilities
P_Na = g_Na * V_t / E_Na * (1 / C_Na_i - 1 / C_Na_o); % in mS/(cm^2*mM)
P_K  = g_K  * V_t / E_K  * (1 / C_K_i  - 1 / C_K_o );
P_Cl =-g_Cl * V_t / E_Cl * (1 / C_Cl_i - 1 / C_Cl_o);


%% Part 1: Linear and nonlinear current approximations %%
V_m = -150:150; % in mV

J_Na_nonlin = (P_Na .* V_m .*(C_Na_i - C_Na_o.* exp(-V_m./V_t)))./(1-exp(- V_m./V_t));  % in uA/cm^2
J_K_nonlin  = (P_K .* V_m .*(C_K_i - C_K_o.* exp(-V_m./V_t)))./(1-exp(- V_m./V_t));
J_Cl_nonlin = (P_Cl .* V_m .*(C_Cl_i - C_Cl_o.* exp(V_m./V_t)))./(1-exp( V_m./V_t));

% Calculate current densities for all V_m
J_Prime_Na_lin = g_Na * (V_m - E_Na);  % in uA/cm^2
J_Prime_K_lin  = g_K  * (V_m - E_K );
J_Prime_Cl_lin = g_Cl * (V_m - E_Cl);

figure();
%plot(V_m, 0*V_m, 'k:');
hold on
% recommended style: linear uses dashes, GHK uses solid lines, same color for the same ions
plot(V_m, J_Na_nonlin, 'r--')
plot(V_m, J_Prime_Na_lin, 'r')
xlabel('V_m (mV)') % V_m in mV as per neuroscience tradition
ylabel('J (\mu A/cm^2)')
title('Problem 1 - Na');
legend('J_{Na} (GHK)', 'J_{Na} (linear)', 'Location', 'SouthEast')
hold off
figure;
hold on;
plot(V_m, J_K_nonlin, 'g--')
plot(V_m, J_Prime_K_lin, 'g')
plot(V_m,(P_K * V_m *C_K_i));
plot(V_m,(P_K * V_m *C_K_o));
xlabel('V_m (mV)') % V_m in mV as per neuroscience tradition
ylabel('J (\mu A/cm^2)')
title('Problem 1 - K');
legend('J_{K} (GHK)', 'J_{K} (linear)','J_{K} (Vm >> Vt)', 'J_{K} (Vm << -Vt)', 'Location', 'SouthEast')
hold off;
figure;
hold on;
plot(V_m, J_Cl_nonlin, 'b--')
plot(V_m, J_Prime_Cl_lin, 'b')
legend('J_{Cl} (GHK)', 'J_{Cl} (linear)', 'Location', 'SouthEast')
xlabel('V_m (mV)') % V_m in mV as per neuroscience tradition
ylabel('J (\mu A/cm^2)')
title('Problem 1 - Cl');
hold off;

%legend('hide')



% figure();
% %plot(V_m, 0*V_m, 'k:');
% hold on
% recommended style: linear uses dashes, GHK uses solid lines, same color for the same ions
% plot(V_m, J_Prime_Na_lin, 'r')
% plot(V_m, J_Prime_K_lin, 'g')
% plot(V_m, J_Prime_Cl_lin, 'b')
% xlabel('V_m (mV)') % V_m in mV as per neuroscience tradition
% ylabel('J (\mu A/cm^2)')
% title('Problem 1');
% %legend('hide')
% legend('J_{Na} (Non-linear)', 'J_K (Non-linear)', 'J_{Cl} (Non-linear)', 'Location', 'SouthEast')
% hold off
%% Part 2: The resting potential %%
% Derive these equations. You will either need to write the equations write in
% Microsoft Word, Latex, or submit a hard copy of hand-written derivations.

% Computing the resting potentials, in mV
V_r_GHK = V_t*log((P_Na*C_Na_o + P_K*C_K_o + P_Cl*C_Cl_i) / (P_Na*C_Na_i + P_K*C_K_i + P_Cl*C_Cl_o));
fprintf("The Resting potential using GHK is %f .",V_r_GHK);
V_r_lin = (g_Na*E_Na + g_K*E_K + g_Cl*E_Cl)/(g_Na + g_K + g_Cl);
fprintf("The Resting potential using Linear approximation is %f .",V_r_lin);
%% Part 3: GHK membrane dynamics %%

%Non-Linear case

% Na+ is active only during first 500 ms
p_Na_t = @(t) P_Na * (t <= 0.5);
% K+ is active during first 250 ms and last 500 ms
p_K_t  = @(t) P_K  * ((t <= 0.25) + (t >= 0.5));
% Cl- is active during first 250 ms and last 500 ms
p_Cl_t = @(t) P_Cl * ((t <= 0.25) + (t >= 0.5));

% Time derivative of V_m, as a function of t and V_m
dVmdt_GHK = @(t,V_m) -1/C_m * (p_Na_t(t).* V_m .* (C_Na_i - C_Na_o .* exp(- V_m / V_t))./(1 - exp(- V_m./V_t))+ p_K_t(t) .* V_m .* (C_K_i - C_K_o .* exp(- V_m / V_t))./(1 - exp(- V_m./V_t)) + p_Cl_t(t) .* V_m .* (C_Cl_i - C_Cl_o .* exp( V_m / V_t))./(1 - exp( V_m./V_t))); 
[t_nonlin, V_m_nonlin] = ode45(dVmdt_GHK, [0 1], V_r_GHK);
figure;

subplot(2,1,1) % top 1/2 panel, V_m
plot(t_nonlin, 0*t_nonlin, 'k:');
hold on
plot(t_nonlin, V_m_nonlin, 'k--') 
ylabel('V_m (mV)')
legend('', 'linear', 'Location', 'SouthEast');
title('Problem 3');
hold off

subplot(6,1,4) % bottom 4/6 panel, g_K
plot(t_nonlin, p_K_t(t_nonlin), 'k');
ylim([-0.1, P_K+0.2]);
ylabel('P_K mS/(cm^2*mM)');

subplot(6,1,5) % bottom 5/6 panel, g_Na
plot(t_nonlin, p_Na_t(t_nonlin), 'k');
ylim([-0.1, P_Na+0.2]);
ylabel('P_{Na} mS/(cm^2*mM)');

subplot(6,1,6) % bottom 6/6 panel, g_Cl
plot(t_nonlin, p_Cl_t(t_nonlin), 'k');
ylim([-0.1, P_Cl+0.2]);
ylabel('P_{Cl} mS/(cm^2*mM)');

xlabel('Time (s)')

% Linear case

% Na+ is active only during first 500 ms
g_Na_t = @(t) g_Na * (t <= 0.5);
% K+ is active during first 250 ms and last 500 ms
g_K_t  = @(t) g_K  * ((t <= 0.25) + (t >= 0.5));
% Cl- is active during first 250 ms and last 500 ms
g_Cl_t = @(t) g_Cl * ((t <= 0.25) + (t >= 0.5));

% Time derivative of V_m, as a function of t and V_m
dVmdt_lin = @(t,V_m) -1/C_m * (g_Na_t(t)*(V_m-E_Na) + g_K_t(t)*(V_m-E_K) + g_Cl_t(t)*(V_m-E_Cl));
[t_lin, V_m_lin] = ode45(dVmdt_lin, [0 1], V_r_lin);

% Plotting
figure;

subplot(2,1,1) % top 1/2 panel, V_m
plot(t_lin, 0*t_lin, 'k:');
hold on
plot(t_lin, V_m_lin, 'k--') 
ylabel('V_m (mV)')
legend('', 'linear', 'Location', 'SouthEast');
title('Problem 3');
hold off

subplot(6,1,4) % bottom 4/6 panel, g_K
plot(t_lin, g_K_t(t_lin), 'k');
ylim([-0.1, g_K+0.2]);
ylabel('g_K (mS/cm^2)');

subplot(6,1,5) % bottom 5/6 panel, g_Na
plot(t_lin, g_Na_t(t_lin), 'k');
ylim([-0.1, g_Na+0.2]);
ylabel('g_{Na} (mS/cm^2)');

subplot(6,1,6) % bottom 6/6 panel, g_Cl
plot(t_lin, g_Cl_t(t_lin), 'k');
ylim([-0.1, g_Cl+0.2]);
ylabel('g_{Cl} (mS/cm^2)');

xlabel('Time (s)')

Dk=1.3*10^-5;
Dna=1.9*10^-5;
r=5*10^-8;
F=9.649 * 10^4;
n=1;
yk = (2.*pi.*Dk.*r.*1.*F.*n);
yNa = (2.*pi.*Dna.*r.*1.*F.*n);

%Problem 5

J_K_Mod = @(n)(P_K.*V_m .*(2.*pi.*Dk.*r.*1.*F.*n).*(C_K_i - (C_K_o .* (exp(-V_m./V_t)))) )./ (((2.*pi.*Dk.*r.*1.*F.*n).*(1-(exp(-V_m./V_t))))+(P_K.*V_m.*(1+(exp(-V_m./V_t)))));%Zk = 1;
J_Na_Mod = @(n)(P_Na.*V_m .*(2.*pi.*Dna.*r.*1.*F.*n).*(C_Na_i - (C_Na_o .* (exp(-V_m./V_t)))) )./ (((2.*pi.*Dna.*r.*1.*F.*n).*(1-(exp(-V_m./V_t))))+(P_Na.*V_m.*(1+(exp(-V_m./V_t)))));%Zk = 1;

figure();

hold on;

plot(V_m,J_K_Mod(1),'g');
plot(V_m,J_Na_Mod(1),'g');


plot(V_m,J_K_Mod(0.1),'r--');
plot(V_m,J_Na_Mod(0.1),'g--');

plot(V_m,J_K_Mod(10),'r:');
plot(V_m,J_Na_Mod(10),'g:');
plot(V_m,J_K_Mod(5),'r-.');
plot(V_m,J_Na_Mod(5),'g-.');
title('Problem 5 - GHK and GHK modified with n = 1,0.1,10')
xlabel('V_m (mV)') % V_m in mV as per neuroscience tradition
ylabel('Jk (\mu A/cm^2)')
legend('J K Mod(1)', 'J Na Mod(1)','J K Mod(0.1)','J Na Mod(0.1)','J K Mod(10)','J Na Mod(10)', 'J K Mod(5)','J Na Mod(5)', 'Location', 'SouthEast')
hold off










