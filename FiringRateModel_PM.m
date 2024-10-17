% FiringRateModel.m: firing rate model for parametric working memory 
% based on negative derivative feedback.
% Created by Sukbin Lim
%
% This program reproduces the experiment in Figs. 2 & 4 in Lim and Goldman
% (2013), Nature Neuroscience

function FiringRateModel_PM

    flag = 1;   % 1 for transient input and 0 for step-like input

    %%% Parameters %%%
    % Intrinsic time constant Ti(ms) for i,j = E (excitatory) or I (inhibitory)
    TE = 20;            TI = 10;            

    % Time constant of exponential filtering of external input
    TEO = 100;       

    % Synaptic time constant Tij (ms) from population j to population i
    TEE = 100;          TIE = 25;           TEI = 10;           TII = 10;           

    % Synaptic connectivity matrix Wij with Gaussian-shaped profiles
    JEE = 150;          JIE = 150;          JEI = 300;          JII = 300; 

%     % Perturbation
%     p = 0.1;
%     JEE = 150*(1+p);    JIE = 150*(1+p);    JEI = 300*(1+p);    JII = 300*(1+p); 

    % Strength of external input
    if (flag ==1)
        JEO = 1500;       % JEO = 1500, 3000 or 4500 for transient input
    else
        JEO = 100;        % JEO = 100, 200 or 300 for step-like input
    end
    JIO = 0*JEO;

    %%% Simulation %%%
    options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-5 1e-5 1e-5 1e-5 1e-5]);
    [T,Y] = ode45(@dtest,0:1:5500,[0 0 0 0 0 0 0],options,flag,JEE,JIE,JEI,JII,TE,TI,TEE,TIE,TEI,TII,TEO,JEO,JIO);

    RE  = Y(:,1); RI  = Y(:,2); SEE = Y(:,3); SIE = Y(:,4); SEI = Y(:,5); SII = Y(:,6); SO = Y(:,7);

    %%% Figure %%%
    LineWidth = 1;
    FontSize = 12;

    % Time course of firing rate of excitatory population 
    figure(1);
    plot(T,RE,'k','LineWidth',LineWidth)
    xlim([0 T(end)]);ylim([0 100]);
    set(gca,'FontSize',FontSize);
    title('Time course of excitatory population')
    xlabel('Time (s)');ylabel('Firing Rate (Hz)');
    set(gca,'Xtick',500:1000:6000)
    set(gca,'XTickLabel',0:1:6)
    set(gca,'Ytick',0:50:100)

    % Time course of inputs to excitatory population
    figure(2);
    plot(T,JEE*SEE,'b','LineWidth',LineWidth);hold on
    plot(T,JEI*SEI,'r','LineWidth',LineWidth);
    plot(T,JEO*SO,'k','LineWidth',LineWidth);
    plot(T,JEE*SEE-JEI*SEI+JEO*SO,'k--','LineWidth',LineWidth);hold off
    xlim([0 T(end)]);
    set(gca,'FontSize',FontSize);
    legend('Excitation','Inhibition','Ext. Input','Total Input');
    xlabel('Time(s)'); ylabel('Input');
    set(gca,'Xtick',500:1000:6000)
    set(gca,'XTickLabel',0:1:6)

function dy = dtest(t,y,flag,JEE,JIE,JEI,JII,TE,TI,TEE,TIE,TEI,TII,TEO,JEO,JIO)
    dy  = zeros(7,1);

    RE  = y(1);
    RI  = y(2);
    SEE = y(3);
    SIE = y(4);
    SEI = y(5);
    SII = y(6);
    SO = y(7);

%     % Nonlinear input-output transfer
%     N = 2;      th = 10;        sig = 30;       maxf = 100;
%     f = @(x) maxf*(x-th)^N/(sig^N+(x-th)^N)*(x>th);
    % Linear input-output transfer
    f = @(x) x;

    if (flag ==1)
        dSO    = 1/TEO*(-SO + (t>500)*(t<600));    % Transient input
    else
        dSO    = 1/TEO*(-SO + (t>500));            % Step-like input
    end

    dRE   = 1/TE * (-RE + f(JEO*SO + JEE*SEE - JEI*SEI));
    dRI   = 1/TI * (-RI + f(JIO*SO + JIE*SIE - JII*SII));
    dSEE  = 1/TEE* (-SEE + RE);
    dSIE  = 1/TIE* (-SIE + RE);
    dSEI  = 1/TEI* (-SEI + RI);
    dSII  = 1/TII* (-SII + RI);

    dy(1) = dRE;
    dy(2) = dRI;
    dy(3) = dSEE;
    dy(4) = dSIE;
    dy(5) = dSEI;
    dy(6) = dSII;
    dy(7) = dSO;
