clear all;
close all;
% Constants
figuresdir = 'D:\Figures\Fig'; 
C_m = 1.0; % membrane capacitance, in uF/cm^2
g_Na = 120.0; % maximum conducances, in mS/cm^2
g_K = 36.0;
g_L = 0.3;
E_Na = 45.0; % Nernst reversal potentials, in mV
E_K = -82.0;
E_L = -59.387;
E_Cl = -80.0; % inhibitory synapse
E_ex = -38.0; % excitatory synapse
% Channel gating kinetics
% Functions of membrane voltage
alpha_m = @(V) 0.1.*(V+45.0)./(1.0 - exp(-(V+45.0) ./ 10.0));
beta_m = @(V) 4.0.*exp(-(V+70.0) ./ 18.0);
alpha_h = @(V) 0.07.*exp(-(V+70.0) ./ 20.0);
beta_h = @(V) 1.0./(1.0 + exp(-(V+40.0) ./ 10.0));
alpha_n = @(V) 0.01.*(V+60.0)./(1.0 - exp(-(V+60.0) ./ 10.0));
beta_n = @(V) 0.125.*exp(-(V+70) ./ 80.0);
% Alpha / Beta constants for inhibitory / excitatory synapses
alpha_r_i = 5.0; % 1/(mM*ms)
beta_r_i = 0.18; % 1/ms
alpha_r_e = 2.4; % 1/(mM*ms)
beta_r_e = 0.56; % 1/ms

I_Na = @(V,m,h) g_Na .* m.^3 .* h .* (V - E_Na);
I_K = @(V, n) g_K .* n.^4 .* (V - E_K);
I_L = @(V) g_L .* (V - E_L);
I_syn_i = @(V,g_GABA,r) r*g_GABA.*(V-E_Cl);
% V=V_post, r is for pre
I_syn_e = @(V,g_Glu,r) r*g_Glu .*(V-E_ex);
% V=V_post, r is for pre

% Membrane voltage differential
dV = @(V, m, h, n, r_i, r_e, I_ext, g_GABA, g_Glu)...
(-I_Na(V, m, h) - I_K(V, n) - I_L(V) - I_syn_i(V, g_GABA, r_i) -I_syn_e(V, g_Glu, r_e) + I_ext) ./ C_m;
% [T] equation for synapses
T_max_i = 1.5; % mM (inhibitory)
T_max_e = 1.0; % mM (excitatory)
K_p = 5.0; % mV
V_p = 7.0; % mV
T_i = @(V_pre) T_max_i ./ (1.0 + exp(-(V_pre - V_p)./K_p));
T_e = @(V_pre) T_max_e ./ (1.0 + exp(-(V_pre - V_p)./K_p));
% Differential gating equations
dm = @(V,m) alpha_m(V).*(1.0-m)-beta_m(V).*m;
dh = @(V,h) alpha_h(V).*(1.0-h)-beta_h(V).*h;
dn = @(V,n) alpha_n(V).*(1.0-n)-beta_n(V).*n;
dr_i = @(V,r) alpha_r_i.*T_i(V).*(1.0-r)-beta_r_i.*r;
% V=V_pre, r is for pre
dr_e = @(V,r) alpha_r_e.*T_e(V).*(1.0-r)-beta_r_e.*r;

% MATLAB
d_single = @(t, x, I_ext, g_GABA, g_Glu) ...
 [ dV(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),I_ext(t),g_GABA,g_Glu); ...
 dm(x(1,:),x(2,:)); ...
 dh(x(1,:),x(3,:)); ...
 dn(x(1,:),x(4,:)); ...
 dr_i(x(1,:),x(5,:));...
 dr_e(x(1,:),x(6,:))];

fig= 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1 - Uncoupled Neurons
% MATLAB
fig= 1;
d = @(t, x, I_exts, g_GABA, g_Glu) ...
 reshape(d_single(t, reshape(x, 6, length(x)/6), I_exts, g_GABA, g_Glu), length(x), 1);
% Shortcut for simulating a whole network
T = linspace(0, 500, 5000); % each division is 0.1 ms
network = @(I_exts, g_GABA, g_Glu) ...
 ode45(d, T, zeros(1, size(g_GABA, 1).*6), [], I_exts, g_GABA, g_Glu);

% One neuron gets 10 uA (starting at 0ms) and the other gets 20 uA
I_exts_ = [10 20]; % uA/cm^2
I_exts = @(t) I_exts_;
% g_GABA and g_Glu matrices are all 0 for no connections
[t,v] = network(I_exts, [0 0; 0 0], [0 0; 0 0]);
v = v';
% MATLAB
V = v(1:6:end,:);
r1=v(5:6:end,:);
r2=v(6:6:end,:);
figure;
plot(t(4001:end), V(:,4001:end));
title('Uncoupled Neurons');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
legend('A - Neuron','B - Neuron');
saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
fig=fig+1;
figure;

subplot(3,1,1);
plot(t(4001:end), V(:,4001:end)); % plot last 100ms
title('Uncoupled Neurons');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
legend('A - Neuron','B - Neuron');
subplot(3,1,2);
plot(t(4001:end), r1(:,4001:end));
title('Uncoupled Neurons');
xlabel('Time (ms)');
ylabel('r_i');
legend('A - Neuron','B - Neuron');
subplot(3,1,3);
plot(t(4001:end), r2(:,4001:end));
title('Uncoupled Neurons');
xlabel('Time (ms)');
ylabel('r_e');
legend('10uA - Neuron','20uA - Neuron');
saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
fig=fig+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2: Two Coupled Neurons -inhibitory

% One neuron gets 10 uA (starting at 0ms) and the other gets 20 uA
I_exts_ = [10 20]; % uA/cm^2
I_exts = @(t) I_exts_;
freq=[];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = [0 0.1 0.2 0.3 0.4 0.5 1.0 1.5 2.0 2.5 3.0 3.5]; % mS/cm^2
m = zeros(2, length(g_GABAs)); % means

for j=1:length(g_GABAs)
    freq_spk=[];
    g_GABA_plot=[];
    for i=1:length(I_exts_)
    [t,v] = network(I_exts, [0 g_GABAs(j); g_GABAs(j) 0], [0 0; 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    subplot(3,1,1);
    plot(t(4001:end), V(:,4001:end));
    title(['2) Coupled Inhibitory Neurons ->  g_GABA -' num2str(g_GABAs(j))]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    legend('A - Neuron 10 uA/cm^2','B - Neuron 20 uA/cm^2');
    subplot(3,1,2)
    plot(t(4001:end), r1(:,4001:end));
    xlabel('Time (ms)');
    ylabel('r_i');
    subplot(3,1,3)
    plot(t(4001:end), r2(:,4001:end));
    xlabel('Time (ms)');
    ylabel('r_e');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
    fig=fig+1;
    m = isi(t(2500:end), V(i,2500:end)); % skip the first 250ms
    freq_spk=horzcat(freq_spk,1000/(m));
    end
    freq=vertcat(freq,freq_spk);
end
figure;
hold on;
plot(g_GABAs',freq(:,1),'r');
plot(g_GABAs',freq(:,2),'b');
plot(g_GABAs',abs(freq(:,1)-freq(:,2)),'g-');
 title('2) Spike Frequency vs Gaba Plot');
    xlabel('GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron 10 uA/cm^2','B - Neuron 20 uA/cm^2','Absolute Freq diff');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
    fig=fig+1;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 3: 
I_exts_ = [10 10.1]; % uA/cm^2
I_exts = @(t) I_exts_ .*(t>=50);
% g_GABA has reciprocal connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = 1; % mS/cm^2
beta_r_i = [0.5 0.4 0.3 0.2 0.1];
p = zeros(1, length(beta_r_i)); % means
for j=1:length(beta_r_i)
    dr_i_3 = @(V,r) alpha_r_i.*T_i(V).*(1.0-r)-beta_r_i(j).*r;
     d_single_3 = @(t, x, I_ext, g_GABA, g_Glu) ...
     [dV(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),I_ext(t),g_GABA,g_Glu); ...
     dm(x(1,:),x(2,:)); ...
     dh(x(1,:),x(3,:)); ...
     dn(x(1,:),x(4,:)); ...
     dr_i_3(x(1,:),x(5,:));...
     dr_e(x(1,:),x(6,:))];
    d_3 = @(t, x, I_exts, g_GABA, g_Glu) ...
    reshape(d_single_3(t, reshape(x, 6, length(x)/6), ...
    I_exts, g_GABA, g_Glu), length(x), 1);
    % Shortcut for simulating a whole network
    T = linspace(0, 500, 5000); % each division is 0.1 ms
    
    network_3 = @(I_exts, g_GABA, g_Glu) ...
    ode45(d_3, T, zeros(1, size(g_GABA, 1).*6), [], I_exts, g_GABA, g_Glu);
    [t,v] = network_3(I_exts, [0 g_GABAs;g_GABAs 0], [0 0; 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
     subplot(3,1,1);
    plot(t(4001:end), V(:,4001:end));
    title(['3-1)Coupled Inhibitory Neurons -> beta r -' num2str(beta_r_i(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron');
    subplot(3,1,2);
    plot(t(4001:end), r1(:,4001:end));
    xlabel('Time (ms)');
ylabel('r_i');
    subplot(3,1,3);
    plot(t(4001:end), r2(:,4001:end));
    xlabel('Time (ms)');
ylabel('r_e');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
    fig=fig+1;
    p(j) = spk_phase(t(2500:end), V(1,2500:end), V(2,2500:end)); % skip first 250ms
end
figure;
plot(beta_r_i,mod(p,2*pi),'r');
title('3-1) Coupled Inhibitory Neurons phase difference');
xlabel('beta r i');
ylabel('Phase Difference');
  saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
fig=fig+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 3: 
I_exts_ = [10 10.1]; % uA/cm^2
I_exts = @(t) I_exts_;
% g_GABA has reciprocal connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = 1; % mS/cm^2
beta_r_i = [0.5:-0.01:0.1];
p = zeros(1, length(beta_r_i)); % means
for j=1:length(beta_r_i)
    dr_i_3 = @(V,r) alpha_r_i.*T_i(V).*(1.0-r)-beta_r_i(j).*r;
     d_single_3 = @(t, x, I_ext, g_GABA, g_Glu) ...
     [dV(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),I_ext(t),g_GABA,g_Glu); ...
     dm(x(1,:),x(2,:)); ...
     dh(x(1,:),x(3,:)); ...
     dn(x(1,:),x(4,:)); ...
     dr_i_3(x(1,:),x(5,:));...
     dr_e(x(1,:),x(6,:))];
    d_3 = @(t, x, I_exts, g_GABA, g_Glu) ...
    reshape(d_single_3(t, reshape(x, 6, length(x)/6), ...
    I_exts, g_GABA, g_Glu), length(x), 1);
    % Shortcut for simulating a whole network
    T = linspace(0, 500, 5000); % each division is 0.1 ms
    
    network_3 = @(I_exts, g_GABA, g_Glu) ...
    ode45(d_3, T, zeros(1, size(g_GABA, 1).*6), [], I_exts, g_GABA, g_Glu);
    [t,v] = network_3(I_exts, [0 g_GABAs;g_GABAs 0], [0 0; 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
%     figure;
%     subplot(3,1,1)
%     plot(t(4001:end), V(:,4001:end));
%     xlabel('Time (ms)');
%     ylabel('Volatage (mV)');
%      title(['3-2) Coupled Inhibitory Neurons -> beta r -' num2str(beta_r_i(j))]);
%     
%     legend('A - Neuron','B - Neuron');
%     subplot(3,1,2)
%     plot(t(4001:end), r1(:,4001:end));
%     xlabel('Time (ms)');
% ylabel('r_i');
%     subplot(3,1,3)
%     plot(t(4001:end), r2(:,4001:end));
%    xlabel('Time (ms)');
% ylabel('r_e');
%     saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
%     fig=fig+1;
    p(j) = spk_phase(t(2500:end), V(1,2500:end), V(2,2500:end)); % skip first 250ms
    end
    
figure;
plot(beta_r_i,mod(p,2*pi),'ro-');
title('3-2) Coupled Inhibitory Neurons phase difference');
xlabel('beta r i');
ylabel('Phase Difference');
  saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 4: 
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GLUs = [0 0.1 0.2 0.3 0.4 0.5]; % mS/cm^2
m_2 = zeros(2, length(g_GLUs)); % means
freq_2=[];
for j=1:length(g_GLUs)
    freq_spk_2=[];
    for i = 1:length(I)
    I_exts_ = [I(i) 0]; % uA/cm^2
    for k =1:length(I_exts_)
    I_exts = @(t) I_exts_ .*(t>=100 && t<=350);
    [t,v] = network(I_exts, [0 0;0 0], [0 g_GLUs(j);0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    subplot(3,1,1);
    plot(t(1:end), V(:,1:end));
    title(['4)Coupled Excitatory Neurons -> I_ext -' num2str(I(i)) ' GLU-' num2str(g_GLUs(j))]);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    legend('A - Neuron','B - Neuron');
    subplot(3,1,2)
    plot(t(1:end), r1(:,1:end));
    xlabel('Time (ms)');
    ylabel('r_i');
    subplot(3,1,3)
    plot(t(1:end), r2(:,1:end));
    xlabel('Time (ms)');
    ylabel('r_e');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
      m_2 = isi(t(2500:3500), V(k,2500:3500));
     freq_spk_2=horzcat(freq_spk_2,1000/(m_2));
    end
    end
    freq_2=vertcat(freq_2,freq_spk_2);
end
figure;
hold on;
plot(g_GLUs',freq_2(:,1),'r');
plot(g_GLUs',freq_2(:,2),'b');
 title('4) Spike Frequency vs GLU Plot 10uA/cm^2');
    xlabel('GLU');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 5: 
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = [0.1 0.2 0.3 0.4 0.5];
g_GLUs = [0.1 0.2 0.3 0.4 0.5]; % mS/cm^2
m_3 = zeros(3, length(g_GLUs)); % means
freq_3=[];
for j=1:length(g_GLUs)
    freq_spk_3=[];
    for i = 1:length(I)
    I_exts_ = [I(i) 0 0]; % uA/cm^2
    for k=1:length(I_exts_)
    I_exts = @(t) I_exts_ .*(t>=100 && t<=350);
    [t,v] = network(I_exts, [0 0 0;0 0 g_GABAs(j);0 0 0],...
        [0 g_GLUs(j) g_GLUs(j);0 0 0;0 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
     subplot(3,1,1);
    plot(t(1:end), V(:,1:end));
    title(['5) Coupled Neurons Circuits Feed forward -> I_ext -'...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j)) ' GABA-' num2str(g_GABAs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron');
    subplot(3,1,2)
    plot(t(1:end), r1(:,1:end));
    xlabel('Time (ms)');
ylabel('r_i');
    subplot(3,1,3)
    plot(t(1:end), r2(:,1:end));
   xlabel('Time (ms)');
ylabel('r_e');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
    m_3 = isi(t(2500:3500), V(k,2500:3500)); % skip the first 250ms
    freq_spk_3=horzcat(freq_spk_3,1000/(m_3));
    end
    end
    freq_3=vertcat(freq_3,freq_spk_3);
end
figure;
hold on;
plot(g_GLUs',freq_3(:,1),'r');
plot(g_GLUs',freq_3(:,2),'b');
plot(g_GLUs',freq_3(:,3),'g');
 title('5) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron');
     saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 6:
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = [0.1 0.2 0.3 0.4 0.5];
g_GLUs = [0.1 0.2 0.3 0.4 0.5]; % mS/cm^2
m_4 = zeros(3, length(g_GLUs)); % means
freq_4=[];
for j=1:length(g_GLUs)
    freq_spk_4=[];
    for i = 1:length(I)
    I_exts_ = [I(i) 0 0]; % uA/cm^2
    for k=1:length(I_exts_)
    I_exts = @(t) I_exts_ .*(t>=100 && t<=350);
    [t,v] = network(I_exts, [0 0 0;0 0 g_GABAs(j);0 0 0], [0 0 g_GLUs(j);0 0 0;0 g_GLUs(j) 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
     subplot(3,1,1);
    plot(t(1:5000), V(:,1:5000));
    title(['6) Coupled Neurons Circuits Feed forward -> I_ext -'...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j)) ' GABA-' num2str(g_GABAs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron');
    subplot(3,1,2)
    plot(t(1:5000), r1(:,1:5000));
    xlabel('Time (ms)');
ylabel('r_i');
    subplot(3,1,3)
    plot(t(1:5000), r2(:,1:5000));
   xlabel('Time (ms)');
ylabel('r_e');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
    m_4 = isi(t(2500:3500), V(k,2500:3500)); % skip the first 250ms
    freq_spk_4=horzcat(freq_spk_4,1000/(m_4));
    end
    end
    freq_4=vertcat(freq_4,freq_spk_4);
end
figure;
hold on;
plot(g_GLUs',freq_4(:,1),'r');
plot(g_GLUs',freq_4(:,2),'b');
plot(g_GLUs',freq_4(:,3),'g');
 title('6) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron');
     saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = [0.1 0.2 0.3 0.4 0.5];
g_GLUs = [0.3]; % mS/cm^2
m_5 = zeros(5, length(g_GLUs)); % means
freq_5=[];

for j=1:length(g_GLUs)
    freq_spk_5=[];
    for i =i: length(I)
     I_exts_ = [I(i) 0 0 0 0]; % uA/cm^2
     for k =1: length(I_exts_)
    I_exts = @(t) I_exts_ .*(t>=100 && t<=101);
    [t,v] = network(I_exts, [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;...
        0 0 0 0 0], [0 g_GLUs(j) 0 0 0;0 0 g_GLUs(j) 0 0;0 0 0 g_GLUs(j) 0;...
        0 0 0 0 g_GLUs(j);g_GLUs(j) 0 0 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    plot(t(1:end), V(:,1:end));
    title(['8) Coupled Neuron Circuits Feedback -> I_ext -' ...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures',['figure' num2str(fig) '.png']));
fig=fig+1;
     m_5 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
    freq_spk_5=horzcat(freq_spk_5,1000/(m_5));
     end
    end
    freq_5=vertcat(freq_5,freq_spk_5);
end
figure;
hold on;
plot(g_GLUs',freq_5(:,1),'r');
plot(g_GLUs',freq_5(:,2),'b');
plot(g_GLUs',freq_5(:,3),'g');
plot(g_GLUs',freq_5(:,4),'k');
plot(g_GLUs',freq_5(:,5),'c');
 title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
    hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
% I = [10];
% % g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% % We will test the following g_GABA values:
% g_GABAs = [0.1 0.2 0.3 0.4 0.5];
% g_GLUs = [0.3]; % mS/cm^2
% m_5 = zeros(5, length(g_GLUs)); % means
% freq_5=[];
% 
% for j=1:length(g_GLUs)
%     freq_spk_5=[];
%     for i =i: length(I)
%      I_exts_ = [I(i) 0 0 0 0]; % uA/cm^2
%      for k =1: length(I_exts_)
%     I_exts = @(t) I_exts_ .*(t>=100 && t<=101);
%     [t,v] = network(I_exts, [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;...
%         0 0 0 0 0], ...
%         [0 g_GLUs(j) 0 0 0;...
%          0 0 g_GLUs(j) 0 0;...
%          0 0 0 g_GLUs(j) 0;...
%          0 0 0 0 g_GLUs(j);...
%          g_GLUs(j) 0 0 0 0]);
%     v = v';
%     % MATLAB
%     V = v(1:6:end,:);
%     r1=v(5:6:end,:);
%     r2=v(6:6:end,:);
%     figure;
%     plot(t(1:end), V(:,1:end));
%     title(['8) Coupled Neuron Circuits Feedback 5 Neuron-> I_ext -' ...
%         num2str(I(i)) ' GLU-' num2str(g_GLUs(j)) ' GABA-' num2str(g_GABAs(j))]);
%     xlabel('Time (ms)');
%     ylabel('Volatage (mV)');
%     legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
% saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
% fig=fig+1;
%      m_5 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
%     freq_spk_5=horzcat(freq_spk_5,1000/(m_5));
%      end
%     end
%     freq_5=vertcat(freq_5,freq_spk_5);
% end
% figure;
% hold on;
% plot(g_GLUs',freq_5(:,1),'r');
% plot(g_GLUs',freq_5(:,2),'b');
% plot(g_GLUs',freq_5(:,3),'g');
% plot(g_GLUs',freq_5(:,4),'k');
% plot(g_GLUs',freq_5(:,5),'c');
%  title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
%     xlabel('GLU/GABA');
%     ylabel('spike freq (Hz)');
%     legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
%     saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
%     fig=fig+1;
%     hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GABAs = [0.1 0.2 0.3 0.4 0.5];
g_GLUs = [0.3]; % mS/cm^2
m_5 = zeros(5, length(g_GLUs)); % means
freq_5=[];

for j=1:length(g_GLUs)
    freq_spk_5=[];
    for i =i: length(I)
     I_exts_ = [I(i) 0 0 0 0]; % uA/cm^2
     for k =1: length(I_exts_)
         
    I_exts = @(t) I_exts_ .*((t>=100 && t<=101)||(t>=395 && t<=410));
    [t,v] = network(I_exts, [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;...
        0 0 0 0 0], ...
        [0 g_GLUs(j) 0 0 0;...
         0 0 g_GLUs(j) 0 0;...
         0 0 0 g_GLUs(j) 0;...
         0 0 0 0 g_GLUs(j);...
         g_GLUs(j) 0 0 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    plot(t(1:end), V(:,1:end));
    title(['8) Coupled Neuron Circuits Feedback 5 Neuron-> I_ext -' ...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j)) ' GABA-' num2str(g_GABAs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
    m_5 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
    freq_spk_5=horzcat(freq_spk_5,1000/(m_5));
    end
    end
    freq_5=vertcat(freq_5,freq_spk_5);
end
figure;
hold on;
plot(g_GLUs',freq_5(:,1),'r');
plot(g_GLUs',freq_5(:,2),'b');
plot(g_GLUs',freq_5(:,3),'g');
plot(g_GLUs',freq_5(:,4),'k');
plot(g_GLUs',freq_5(:,5),'c');
 title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
    hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
I = [10];
beta_r_e=0.1;
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GLUs = [0.5]; % mS/cm^2
m_7 = zeros(3, length(g_GLUs)); % means
freq_7=[];

for j=1:length(g_GLUs)
    freq_spk_7=[];
    for i =i: length(I)
     I_exts_ = [I(i) 0 0]; % uA/cm^2
     for k =1: length(I_exts_)
    I_exts = @(t) I_exts_ .*(t>=100 && t<=101);
    dr_e_4 = @(V,r) alpha_r_e.*T_e(V).*(1.0-r)-beta_r_e.*r;
     d_single_4 = @(t, x, I_ext, g_GABA, g_Glu) ...
     [dV(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),I_ext(t),g_GABA,g_Glu); ...
     dm(x(1,:),x(2,:)); ...
     dh(x(1,:),x(3,:)); ...
     dn(x(1,:),x(4,:)); ...
     dr_i(x(1,:),x(5,:));...
     dr_e_4(x(1,:),x(6,:))];
    d_4 = @(t, x, I_exts, g_GABA, g_Glu) ...
    reshape(d_single_4(t, reshape(x, 6, length(x)/6), ...
    I_exts, g_GABA, g_Glu), length(x), 1);
    % Shortcut for simulating a whole network
    T = linspace(0, 500, 5000); % each division is 0.1 ms
    
    network_4 = @(I_exts, g_GABA, g_Glu) ...
    ode45(d_4, T, zeros(1, size(g_GABA, 1).*6), [], I_exts, g_GABA, g_Glu);
    [t,v] = network_4(I_exts, [0 0 0;0 0 0;0 0 0], ...
        [0 g_GLUs(j) 0;...
         0 0 g_GLUs(j);...
         g_GLUs(j) 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    plot(t(1:end), V(:,1:end));
    title(['8) Coupled Neuron Circuits Feedback 3 Neuron -> I_ext -' ...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j)) ' GABA-' num2str(g_GABAs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron');
% saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
% fig=fig+1;
     m_7 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
    freq_spk_7=horzcat(freq_spk_7,1000/(m_7));
     end
    end
    freq_7=vertcat(freq_7,freq_spk_7);
end
figure;
hold on;
plot(g_GLUs',freq_7(:,1),'r');
plot(g_GLUs',freq_7(:,2),'b');
plot(g_GLUs',freq_7(:,3),'g');

 title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron');
%     saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
%     fig=fig+1;
    hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
beta_r_e=0.36;
g_GABAs = [0.1 0.2 0.3 0.4 0.5];
% g_GLUs = [0.1 0.2 0.3 0.4 0.5]; % mS/cm^2
g_GLUs=[0.2];
m_6 = zeros(4, length(g_GLUs)); % means
freq_6=[];

for j=1:length(g_GLUs)
    freq_spk_6=[];
    for i =i: length(I)
     I_exts_ = [I(i) 0 0 0]; % uA/cm^2
     for k =1: length(I_exts_)
    I_exts = @(t) I_exts_ .*((t>=100 && t<=101)||(t>=395 && t<=410));
     dr_e_4 = @(V,r) alpha_r_e.*T_e(V).*(1.0-r)-beta_r_e.*r;
     d_single_4 = @(t, x, I_ext, g_GABA, g_Glu) ...
     [dV(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),I_ext(t),g_GABA,g_Glu); ...
     dm(x(1,:),x(2,:)); ...
     dh(x(1,:),x(3,:)); ...
     dn(x(1,:),x(4,:)); ...
     dr_i(x(1,:),x(5,:));...
     dr_e_4(x(1,:),x(6,:))];
    d_4 = @(t, x, I_exts, g_GABA, g_Glu) ...
    reshape(d_single_4(t, reshape(x, 6, length(x)/6), ...
    I_exts, g_GABA, g_Glu), length(x), 1);
    % Shortcut for simulating a whole network
    T = linspace(0, 500, 5000); % each division is 0.1 ms
    
    network_4 = @(I_exts, g_GABA, g_Glu) ...
    ode45(d_4, T, zeros(1, size(g_GABA, 1).*6), [], I_exts, g_GABA, g_Glu);
    
    [t,v] = network_4(I_exts, [0 0 0 0 ;0 0 0 0 ;0 0 0 0 ;0 0 0 0],...
        [0 g_GLUs(j) 0 0 ;...
         0 0 g_GLUs(j) 0 ;...
         0 0 0 g_GLUs(j);...
        g_GLUs(j) 0 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    plot(t(1:end), V(:,1:end));
    title(['8) Coupled Neuron Circuits Feedback 4 Neuron -> I_ext -' ...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron');
saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
     m_6 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
    freq_spk_6=horzcat(freq_spk_6,1000/(m_6));
     end
    end
    freq_6=vertcat(freq_6,freq_spk_6);
end
figure;
hold on;
plot(g_GLUs',freq_6(:,1),'r');
plot(g_GLUs',freq_6(:,2),'b');
plot(g_GLUs',freq_6(:,3),'g');
plot(g_GLUs',freq_6(:,4),'k');
 title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
    hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
I = [10];
beta_r_e=0.36;
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
g_GLUs = [0.2]; % mS/cm^2
m_7 = zeros(3, length(g_GLUs)); % means
freq_7=[];

for j=1:length(g_GLUs)
    freq_spk_7=[];
    for i =i: length(I)
     I_exts_ = [I(i) 0 0]; % 5/cm^2
     for k =1: length(I_exts_)
    I_exts = @(t) I_exts_ .*(t>=100 && t<=101);
    [t,v] = network(I_exts, [0 0 0;0 0 0;0 0 0], ...
        [0 g_GLUs(j) 0;...
         0 0 g_GLUs(j);...
         g_GLUs(j) 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    plot(t(1:end), V(:,1:end));
    title(['8) Coupled Neuron Circuits Feedback 3 Neuron -> I_ext -' ...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron');
saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
     m_7 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
    freq_spk_7=horzcat(freq_spk_7,1000/(m_7));
     end
    end
    freq_7=vertcat(freq_7,freq_spk_7);
end
figure;
hold on;
plot(g_GLUs',freq_7(:,1),'r');
plot(g_GLUs',freq_7(:,2),'b');
plot(g_GLUs',freq_7(:,3),'g');

 title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
    hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: 
I = [10];
% g_GABA has reciprical connections [0 x; x 0], g_Glu is all 0
% We will test the following g_GABA values:
beta_r_e=0.36;
g_GABAs = [0.1 0.2 0.3 0.4 0.5];
% g_GLUs = [0.1 0.2 0.3 0.4 0.5]; % mS/cm^2
g_GLUs=[0.2];
m_6 = zeros(4, length(g_GLUs)); % means
freq_6=[];

for j=1:length(g_GLUs)
    freq_spk_6=[];
    for i =i: length(I)
     I_exts_ = [I(i) 0 0 0]; % uA/cm^2
     for k =1: length(I_exts_)
    I_exts = @(t) I_exts_ .*((t>=100 && t<=101));
     dr_e_4 = @(V,r) alpha_r_e.*T_e(V).*(1.0-r)-beta_r_e.*r;
     d_single_4 = @(t, x, I_ext, g_GABA, g_Glu) ...
     [dV(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),I_ext(t),g_GABA,g_Glu); ...
     dm(x(1,:),x(2,:)); ...
     dh(x(1,:),x(3,:)); ...
     dn(x(1,:),x(4,:)); ...
     dr_i(x(1,:),x(5,:));...
     dr_e_4(x(1,:),x(6,:))];
    d_4 = @(t, x, I_exts, g_GABA, g_Glu) ...
    reshape(d_single_4(t, reshape(x, 6, length(x)/6), ...
    I_exts, g_GABA, g_Glu), length(x), 1);
    % Shortcut for simulating a whole network
    T = linspace(0, 500, 5000); % each division is 0.1 ms
    
    network_4 = @(I_exts, g_GABA, g_Glu) ...
    ode45(d_4, T, zeros(1, size(g_GABA, 1).*6), [], I_exts, g_GABA, g_Glu);
    
    [t,v] = network_4(I_exts, [0 0 0 0 ;0 0 0 0 ;0 0 0 0 ;0 0 0 0],...
        [0 g_GLUs(j) 0 0 ;...
         0 0 g_GLUs(j) 0 ;...
         0 0 0 g_GLUs(j);...
        g_GLUs(j) 0 0 0]);
    v = v';
    % MATLAB
    V = v(1:6:end,:);
    r1=v(5:6:end,:);
    r2=v(6:6:end,:);
    figure;
    plot(t(1:end), V(:,1:end));
    title(['8) Coupled Neuron Circuits Feedback 4 Neuron -> I_ext -' ...
        num2str(I(i)) ' GLU-' num2str(g_GLUs(j))]);
    xlabel('Time (ms)');
    ylabel('Volatage (mV)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron','E - Neuron');
saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
fig=fig+1;
     m_6 = isi(t(2500:end), V(k,2500:end)); % skip the first 250ms
    freq_spk_6=horzcat(freq_spk_6,1000/(m_6));
     end
    end
    freq_6=vertcat(freq_6,freq_spk_6);
end
figure;
hold on;
plot(g_GLUs',freq_6(:,1),'r');
plot(g_GLUs',freq_6(:,2),'b');
plot(g_GLUs',freq_6(:,3),'g');
plot(g_GLUs',freq_6(:,4),'k');
 title('8) Spike Frequency vs GLU/GABA Plot 10 uA/cm^2');
    xlabel('GLU/GABA');
    ylabel('spike freq (Hz)');
    legend('A - Neuron','B - Neuron','C - Neuron','D - Neuron');
    saveas(gcf,fullfile('C:\Users\ashak\Documents\MATLAB\Figures1',['figure' num2str(fig) '.png']));
    fig=fig+1;
    hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 8: