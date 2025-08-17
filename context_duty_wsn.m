% =========================================================
% Context-Aware Duty Cycling for WSN Lifetime (Ready-to-Run)
% =========================================================
% This script simulates a simple WSN with duty cycling based on
% residual energy, entropy of sensed data, and traffic load.
% It outputs lifetime metrics (FND, HND, LND) and plots them.
%
% Author: ChatGPT Refactor
% =========================================================

clc; clear; close all;

%% PARAMETERS
N = 100;                   % number of nodes
E0 = 0.5;                 % initial energy (J)
field = [100 100];        % area size (m)
BS = [50 150];            % base station location
rounds = 2500;            % maximum rounds

% Radio model params
E_elec = 50e-9;
E_fs   = 10e-12;
E_mp   = 0.0013e-12;
k      = 4000;            % bits per message
d0     = sqrt(E_fs/E_mp);

% Duty control params
Dmin = 0.05; Dmax = 0.6; gamma = 0.5;
betaE = 0.5; betaH = 0.3; betaT = 0.2;
eta   = 0.4; delta = 0.05;

W = 16; B = 10;           % entropy window size, bins
Tmax = 8; alphaGen = 0.2;

%% INITIALIZE NODES
nodes = struct();
for i=1:N
    nodes(i).x = rand()*field(1);
    nodes(i).y = rand()*field(2);
    nodes(i).E = E0;
    nodes(i).D = 0.1;
    nodes(i).Dnext = 0.1;
    nodes(i).xwin = nan(W,1);
    nodes(i).fwdCnt = 0;
    nodes(i).genCnt = 0;
end

%% METRICS
aliveCount = zeros(1,rounds);
deadCount  = zeros(1,rounds);
firstDead = 0; halfDead = 0; lastDead = 0;

%% MAIN LOOP
for r=1:rounds
    alive = sum([nodes.E] > 0);
    dead  = N - alive;
    aliveCount(r) = alive;
    deadCount(r)  = dead;

    % Lifetime markers
    if firstDead==0 && dead>=1, firstDead=r; end
    if halfDead==0 && dead>=N/2, halfDead=r; end
    if lastDead==0 && dead==N, lastDead=r; break; end

    % ---- Sensing & Context Update ----
    for i=1:N
        if nodes(i).E <= 0, continue; end

        % Residual energy normalization
        En = nodes(i).E / E0;

        % Sensing model (synthetic: Gaussian + occasional burst)
        sensed = randn() * 0.1;
        if rand() < 0.01, sensed = sensed + 5*rand(); end % burst event
        nodes(i).xwin = [nodes(i).xwin(2:end); sensed];

        % Entropy
        H = entropyDiscrete(nodes(i).xwin, B);
        Hn = H / log(B);

        % Traffic load (simplified: 1 gen pkt + fwd fraction)
        nodes(i).genCnt = 1;  % always generate one pkt
        nodes(i).fwdCnt = poissrnd(2); % synthetic forwarding
        Tl = nodes(i).fwdCnt + alphaGen*nodes(i).genCnt;
        Tn = min(1, Tl/Tmax);

        % Controller
        S = betaE*(1-En) + betaH*Hn + betaT*Tn;
        Dstar = min(max(Dmin + gamma*S, Dmin), Dmax);

        % Hysteresis + smoothing
        if abs(Dstar - nodes(i).D) > delta
            nodes(i).Dnext = (1-eta)*nodes(i).D + eta*Dstar;
        else
            nodes(i).Dnext = nodes(i).D;
        end

        % ---- Energy consumption ----
        % Node sends to BS directly (simplified model)
        d = sqrt((nodes(i).x-BS(1))^2 + (nodes(i).y-BS(2))^2);
        if d < d0
            Etx = E_elec*k + E_fs*k*d^2;
        else
            Etx = E_elec*k + E_mp*k*d^4;
        end
        Erx = E_elec*k*nodes(i).fwdCnt;

        consume = (Etx + Erx) * nodes(i).D; % scaled by duty factor
        nodes(i).E = nodes(i).E - consume;
        if nodes(i).E < 0, nodes(i).E = 0; end

    end

    % Update duty factors
    for i=1:N
        nodes(i).D = nodes(i).Dnext;
    end
end

%% RESULTS
fprintf('First Node Dead (FND): %d rounds\n', firstDead);
fprintf('Half Nodes Dead (HND): %d rounds\n', halfDead);
fprintf('Last Node Dead (LND): %d rounds\n', lastDead);

%% PLOTS
figure;
plot(1:rounds, aliveCount, 'LineWidth',1.5);
xlabel('Rounds'); ylabel('Alive Nodes');
title('Alive Nodes over Time');
grid on;

figure;
plot(1:rounds, deadCount, 'r','LineWidth',1.5);
xlabel('Rounds'); ylabel('Dead Nodes');
title('Dead Nodes over Time');
grid on;

%% SUPPORT FUNCTION
function H = entropyDiscrete(xwin, B)
    x = xwin(~isnan(xwin));
    if numel(x)<2, H=0; return; end
    xmin=min(x); xmax=max(x);
    if xmax==xmin, H=0; return; end
    edges=linspace(xmin,xmax,B+1);
    cnts=histcounts(x,edges);
    p=cnts/sum(cnts);
    p=p(p>0);
    H=-sum(p.*log(p));
end
