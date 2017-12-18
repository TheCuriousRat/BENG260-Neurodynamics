%% Example code for BENG 260 Homework 1, Fall 2017
% Questions to vramesh@eng.ucsd.edu or gert@ucsd.edu
clear all;
close all;
%% Morris-Lecar (barnacle muscle fiber) Model
% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2
g_Ca = 1.1; % maximum conducances, in mS/cm^2
g_K  = 2.0;
g_L  = 0.5;
E_Ca = 100; % Nernst reversal potentials, in mV
E_K  = -70;
E_L  = -50;

% Channel gating kinetics
% Functions of membrane voltage
m_infty = @(V) (1 + tanh((V + 1) ./ 15)) ./ 2;
w_infty = @(V) (1 + tanh(V ./ 30)) ./ 2;
tau_w   = @(V) 5 ./ cosh(V ./ 60);  % in ms

% Membrane currents (in uA/cm^2)
I_Ca = @(V)    g_Ca .* m_infty(V) .* (V - E_Ca);
I_K  = @(V, w) g_K  .* w          .* (V - E_K);
I_L  = @(V)    g_L                .* (V - E_L);

% External current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2.3: Morris-Lecar neuron
% vector coding of state variables: X = [V, w]
I = [5, 24.3, 30, 200, 250];
V_lim = [-60 120];
w_lim = [-0.2 1.2];

v = linspace(V_lim(1), V_lim(2), 100);
w = linspace(w_lim(1), w_lim(2), 100);
dVdt = @(V,w,I) (I - I_Ca(V) - I_K(V, w) - I_L(V)) ./ C_m;
dwdt=@(V,w)  (w_infty(V) - w) ./ tau_w(V);
[V, W] = meshgrid(v, w);
[X_null_w, Y_null_w] = getNullcline(dwdt(V, W), v, w);
plot_traj=[];
trajectory=[];
for i=1:5
figure;
hold on;
plot(X_null_w, Y_null_w, 'g--');
[X_null_V_1, Y_null_V_1] = getNullcline(dVdt(V, W, I(i)), v, w);

[xint_1(i), yint_1(i)] = getCrossings(X_null_V_1, Y_null_V_1, X_null_w, Y_null_w);

plot(xint_1(i), yint_1(i),'rx');
J_1 = getJacobian(V, W, dVdt(V, W, I(i)), dwdt(V, W), xint_1(i), yint_1(i));
[V1, D1] = eig(J_1);
disp(eig(J_1));% V and D have real and imaginary components
fprintf("%f %f",D1(1,1),D1(2,2));
% J_2 = getJacobian(V, W, dVdt(V, W, I(2)), dwdt(V, W), xint_2, yint_2);
t = 0:0.001:100;
thetas = 0:45:180;
plot(X_null_V_1,Y_null_V_1);
for j = 1:5
 dXdt1 = @(t1, X1) [dVdt(X1(1), X1(2), I(i)); dwdt(X1(1), X1(2))];
 [t1, X1] = ode23(dXdt1, t, [xint_1(i) + cos(thetas(j) * (3.14159)/180), yint_1(i) + sin(thetas(j) * (3.14159)/180)]);
trajectory=horzcat(trajectory,plot(X1(:,1),X1(:,2)));

end


title('Null clines ');
ylabel('w');
xlabel('V (mV)')
%plot(X1(:,1),X1(:,2));
legend("w null cline","Intersection Point","V null cline")
hold off;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 4: Bonus: Morris-Lecar model with expanded channel gating dynamics  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
% 
% I = [5, 24.3, 35, 200];
% 
% C_m  = 1.0; % membrane capacitance, in uF/cm^2
% g_Na = 120; % maximum conducances, in mS/cm^2
% g_K  = 36;
% g_L  = 0.3;
% E_Na =  115; % Nernst reversal potentials, in mV
% E_K  = -12;
% E_L  = 10.613;
% 
% amV = @(V)(25 - V )./(10 .* (exp((25 - V )./10) - 1));
% bmV = @(V) 4 *exp(-V ./18);
% ahV = @(V) 0.07 *exp(-V ./20);
% bhV = @(V) 1./(exp((30 - V )./10) + 1);
% anV = @(V) (10 - V )./(100 * (exp((10 - V )./10) - 1));
% bnV = @(V) 0.125* exp(-V ./80);
% 
% h1 =@(n1) 0.8 - n1;
% m_inf =@(V) amV(V )./(amV(V ) + bmV(V ));
% % I_Na = @(V,n)     g_Na .*(m_inf(V).^3).*h(n).* (V - E_Na);
% I_Na = @(V,m,h)     g_Na .*(m.^3).*h.* (V - E_Na);
% I_K  = @(V,n)  g_K  .* (n.^4) .* (V - E_K);
% I_L  = @(V)     g_L  .* (V - E_L);
% 
% V1_lim = [-60 60];
% n_lim = [-0.2 1.2];
% v1 = linspace(V1_lim(1), V1_lim(2), 100);
% n = linspace(n_lim(1), n_lim(2), 100);
% dV1dt = @(V,n,I) (I - I_Na(V,m_inf(V),h1(n)) - I_K(V, n) - I_L(V)) ./ C_m;
% dndt=@(V,n) anV(V) .*(1 - n) - bnV(V).* n;
% [V1, N] = meshgrid(v1, n);
% [X_null_n, Y_null_n] = getNullcline(dndt(V1,N), v1, n);
% [X_null_V1, Y_null_V1] = getNullcline(dV1dt(V1, N, I(1)), v1, n);
% [X_null_V2, Y_null_V2] = getNullcline(dV1dt(V1, N, I(2)), v1, n);
% [X_null_V3, Y_null_V3] = getNullcline(dV1dt(V1, N, I(3)), v1, n);
% [X_null_V4, Y_null_V4] = getNullcline(dV1dt(V1, N, I(4)), v1, n);
% figure;
% hold on;
% plot(X_null_n, Y_null_n, 'r--');
% plot(X_null_V1, Y_null_V1);
% plot(X_null_V2, Y_null_V2);
% plot(X_null_V3, Y_null_V3);
% plot(X_null_V4, Y_null_V4);


I = [5, 24.3, 30, 200,250];
C_m  = 1.0; % membrane capacitance, in uF/cm^2
g_Na = 120; % maximum conducances, in mS/cm^2
g_K  = 36;
g_L  = 0.3;
E_Na =  115; % Nernst reversal potentials, in mV
E_K  = -12;
E_L  = 10.613;

V1_lim = [-60 120];
n_lim = [-0.2 1.2];

v1 = linspace(V1_lim(1), V1_lim(2), 100);
n = linspace(n_lim(1), n_lim(2), 100);
amV = @(V)(25 - V )./(10 .* (exp((25 - V )./10) - 1));
bmV = @(V) 4 *exp(-V ./18);
ahV = @(V) 0.07 *exp(-V ./20);
bhV = @(V) 1./(exp((30 - V )./10) + 1);
anV = @(V) (10 - V )./(100 * (exp((10 - V )./10) - 1));
bnV = @(V) 0.125* exp(-V ./80);
h =@(n1) 0.8 - n1;
m_inf =@(V) amV(V )./(amV(V ) + bmV(V ));

I_Na = @(V,n)     g_Na .*(m_inf(V).^3).*h(n).* (V - E_Na);
I_K  = @(V,n)  g_K  .* (n.^4) .* (V - E_K);
I_L  = @(V)     g_L  .* (V - E_L);


dV1dt = @(V,n,I)(I - I_Na(V,n) - I_K(V,n) - I_L(V))./ C_m;
dndt= @(V,n)anV(V) .*(1 - n) - bnV(V).* n;

[V1, N] = meshgrid(v1, n);
[X_null_n, Y_null_n] = getNullcline(dndt(V1, N), v1, n);

%plot(X_null_V1,Y_null_V1);
trajectory_hh=[];
for k=1:5
figure;
hold on;
plot(X_null_n, Y_null_n, 'r--');
%[X_null_V1, Y_null_V1] = getNullcline(dV1dt(V1, N, I(1)), v1, n);
[X_null_V1, Y_null_V1] = getNullcline(dV1dt(V1, N, I(k)), v1, n);

[xint_1_hh(k), yint_1_hh(k)] = getCrossings(X_null_V1, Y_null_V1, X_null_n, Y_null_n);
plot(xint_1_hh(k), yint_1_hh(k),'rx');

J_1hh = getJacobian(V1, N, dV1dt(V1, N, I(k)), dndt(V1, N), xint_1_hh(k), yint_1_hh(k));
[V1hh, D1hh] = eig(J_1hh);
disp(eig(J_1hh));
fprintf("%f %f",D1hh(1,1),D1hh(2,2));% V and D have real and imaginary components
plot(X_null_V1,Y_null_V1);
t = 0:0.001:100;
thetas = 0:45:180;
for l = 1:5
 dXdthh = @(t1, X1hh) [dV1dt(X1hh(1), X1hh(2), I(k)); dndt(X1hh(1), X1hh(2))];
 [t1hh, X1hh] = ode23(dXdthh, t, [xint_1_hh(k) + cos(thetas(l) * (3.14159)/180), yint_1_hh(k) + sin(thetas(l) * (3.14159)/180)]);
 trajectory_hh=horzcat(trajectory_hh,plot(X1hh(:,1),X1hh(:,2)));

end


title('V - n Null Cline ');
ylabel('n K+ activation');
xlabel('V (mV)');
%plot(X1hh(:,1),X1hh(:,2),'g');
legend("n null cline","Intersection Point","V null cline")
hold off;

end

