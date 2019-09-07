
% This file provides data descriptions and analyses demos relevant to the
% manuscript: ''Send-receive communication asymmetry in brain networks:
% Inferring directionallity of neural signalling from undirected
% connectomes''.

% Caio Seguin, Adeel Razi, Andrew Zalesky
% Feb 2019

%% Load data

load('sen_rev_data.mat');

%% Data description: asy struct

% The struct asy contains data on related to the measure of send-receive communication
% asymmetry. The struct is organized as folows:
% 
% asy.com - Send-receive communication asymmetry.
% 
%   asy.com.rg - Regional communication asymmetry.
% 
%       asy.com.rg.{n256,n360,n512} - Regional cortical parcellations for
%       N = 256,360,512.
%  
%           asy.com.rg.n256.{d10,d15,d20} - Regional communication asymmetry
%           computed for 10%, 15%, 20% structural connection density.
% 
%               asy.com.rg.n256.d10.{nav,dif,si} - Regional communication asymmetry
%               for navigation efficiency, diffusion effeciency and search
%               information. This are the variables projected onto the
%               cortical surface in Fig 3a,d,g.
% 
%   asy.com.ss - Subsystems communication asymmetry.
% 
%       asy.com.ss.n360.{d10,d15,d20} - Subsystems communication asymmetry for
%       10%, 15%, 20% structural connection density.
% 
%           asy.com.ss.n360.d10.{nav,dif,si} - Subsystems communication asymmetry
%           for navigation efficiency, diffusion effeciency and search
%           information.
% 
%           asy.com.ss.n360.d10.nav.{yeo7,yeo17,hcp22} - Subsystems
%           communication asymmetry for yeo 7, yeo 17, hcp 22 cortical
%           subsystems
% 
%           asy.com.ss.n360.d10.nav.yeo7.{vec,mat} - System-wise and
%           pairwise subsystems communication asymmetry. These are the
%           variables shown in Figs 4a,c (mat) and 4b,d (vec).
% 
% asy.ec - Send-receive effective connectivity asymmetry
% 
%   asy.ec.ss.n360.{yeo7,yeo17,hcp22} - Effective connectivity asymmetry
%   for the yeo 7, yeo 17 and hcp 22 cortical subsystems.
% 
%       asy.ec.ss.n360.yeo7.{ses1,ses2} - Effective connectivity asymmetry
%       for resting-state fMRI sessions 1 and 2
% 
%           asy.ec.ss.n360.yeo7.ses1.mat - Pairwise effective connecitivty
%           asymmetry. These are the variables shown in Figs 5c,d (Eff. conn. asymmetry matrix)

%% Data description: sre struct

% The struct sre constains data on regional Send/Receive efficiencies. The
% the struct field naming convetion follows the one for the asy struct
% described above.

%% Data description: g1 struct

% The g1 struct contains the cortical gradient of functional connectivity
% from Margulies et al 2016, PNAS (https://doi.org/10.1073/pnas.1608282113)
% downsampled from vertice to regional resolution for parcellations with
% N=256,360,512. The original cortical gradient is publicly available at
% https://www.neuroconnlab.org/data/index.html

% g1.n360.val stores the average gradient value of vertices comprised into
% each of 360 cortical regions. Correlations between g1.n360.val and
% communication asymmetries are reported in the Send-receive communication
% asymmetries of the human connectome section.

% g1.n360.group.thr contains a classification of cortical regions into
% unimodal (g1.n360.val <= -2), neutral (-2 < g1.n360.val <= 2) or multimodal
% (g1.n360.val > 2). This classification is used in Fig 3c,f,i.

%% Demo 1: Reproduce analyses used in Fig 3d,e,f (Send-receive communication asymmetries of the human connectome)

% Fig 3d
display('*** Analysis: Fig 3d ***');
display(['Values projected onto the cortical surface are stored in asy.com.rg.n360.d15.si ', ... 
    'and available as an interactive connectome workbench scene in <insert BALSA URL>']);
fprintf('\n');

% Fig 3e
display('*** Analysis: Fig 3e ***');
figure;
scatter(sre.rg.n360.d15.si.rev, sre.rg.n360.d15.si.sen);
refline(1,0);
xlabel('Receiving search information');
ylabel('Sending search information');
title('Fig 3E');
set(gca, 'FontSize', 14);
display('Scatter plot generated.');
fprintf('\n');

% Fig 3f
display('*** Analysis: Fig 3f ***');
% Continuous analysis
[r,p] = corr(asy.com.rg.n360.d15.si, g1.n360.val);
% Categorical analysis
aux = asy.com.rg.n360.d15.si;
bux = g1.n360.group.thr;
[p_u_n, ~, s_u_n] = ranksum(aux(bux == 1), aux(bux == 2), 'tail', 'both');
[p_n_h, ~, s_n_h] = ranksum(aux(bux == 2), aux(bux == 3), 'tail', 'both');
[p_u_h, ~, s_u_h] = ranksum(aux(bux == 1), aux(bux == 3), 'tail', 'both');
fprintf('Relationship between search info asymmetry and G1\n');
fprintf('Correlation: r=%.2f, p=%.0e\n', r, p);
fprintf('Wilcoxon test, Unimodal > Neutral: Z=%.2f, p=%.0e\n', s_u_n.zval , p_u_n);
fprintf('Wilcoxon test, Neutral > Multimodal: Z=%.2f, p=%.0e\n', s_n_h.zval , p_n_h);
fprintf('Wilcoxon test, Unimodal > Multimodal: Z=%.2f, p=%.0e\n', s_u_h.zval , p_u_h);
fprintf('\n');

%% Demo 2: Reproduce analyses in Fig 4c,d (Send-receive communication asymmetries of cortical subsystems)

display('*** Analysis: Fig 4c,d ***');

% Bonferoni corrected significance threshold for pairwise subsystem asymmetries
% alpha=0.05, 22*21/2 multiple comparisons, two-sided test, 200 subjects
pw_sig_thr = -tinv(0.025/(11*21), 199);

% Bonferoni corrected significance threshold for individual subsystem asymmetries
% alpha=0.05, 22 multiple comparisons, two-sided test, 200 subjects
ind_sig_thr = -tinv(0.025/22, 199);

X = asy.com.ss.n360.d15.nav.hcp22.mat;
X(X < pw_sig_thr) = 0;

[asy_sort, idx_sort] = sort(asy.com.ss.n360.d15.nav.hcp22.vec, 'descend');

figure;
subplot(1,2,1);
imagesc(X); colorbar;
title('Fig 4C: Send-receive navigation asymmetry matrix');

subplot(1,2,2);
barh(flipud(asy_sort));
title('Fig 4D: Subsystem send-receive navigation asymmetry');
xlabel('Receivers < Neutral > Senders');

aux = (1:22)';
set(gca, 'YTick', aux);
set(gca, 'YTickLabel', flipud(aux(idx_sort)));
set(gca, 'FontSize', 14);

l_sen = line([ind_sig_thr, ind_sig_thr], get(gca,'YLim'));
l_rev = line(-[ind_sig_thr, ind_sig_thr], get(gca,'YLim'));
l_sen.LineWidth = 2;
l_sen.LineStyle = '--';
l_sen.Color = [1 0.5 0];
l_rev.LineWidth = 2;
l_rev.LineStyle = '--';
l_rev.Color = [0 0.5 1];

display('Plots generated.');
fprintf('\n');

%% Demo 3: Reproduce analyses in Fig 5c (Send-receive communication asymmetry and effective connectivity)

display('*** Analysis: Fig 5c ***');

% Correlation between communication and effective connectivity asymmetry
% for navigation asymmetry, 15% connection density, Yeo 17, fMRI session 1 (Fig
% 5c).

% Upper triangle index
upt = find(triu(ones(17),1));
[r,p] = corr(asy.com.ss.n360.d15.nav.yeo17.mat(upt), asy.ec.ss.n360.yeo17.ses1.mat(upt));

figure;
scatter(asy.com.ss.n360.d15.nav.yeo17.mat(upt), asy.ec.ss.n360.yeo17.ses1.mat(upt));
lsline;
xlabel('Send-receive navigation asymmetry');
ylabel('Effective connectivity asymmetry');
title(sprintf('Fig 5c: r=%.2f, p=%.0e', r, p));
set(gca, 'FontSize', 14);

display('Scatter plot generated');
fprintf('\n');

%% Demo 4: Reproduce analyses in Fig 6 (Send-receive communication asymmetry of directed non-human connectomes)

display('*** Analysis: Fig 6 (mouse results) ***');

display('Load data from Rubinov et al 2015 PNAS. https://doi.org/10.1073/pnas.1420315112');
load('mouse_data');

% Directed (asymmetric) connectivity matrix
W_dir = W;
% Euclidean distance matrix between region centroids
D = squareform(pdist(coords));

display('Step-by-step computation of network communication measures and send-receive asymmetries.');

%%% Step 1: Create an undirected (symmetric) macaque connectome
W_und = (W_dir + W_dir')./2;

%%% Step 2: Weight-to-length transformation

% Directed
L_dir = -log10(W_dir./(max(W_dir(:)) + min(W_dir(W_dir>0))));
L_dir(L_dir == Inf) = 0;

% Undirected
L_und = -log10(W_und./(max(W_und(:)) + min(W_und(W_und>0))));
L_und(L_und == Inf) = 0;

%%% Step 3: Compute weighted communication asymmetries for the directed (original) connectome

% Weighted navigation efficiency
[~, ~, PL_nav] = navigate(L_dir, D);
E_nav = 1./PL_nav;
dir.E_nav_asy = (E_nav - E_nav')./(E_nav + E_nav');

% Weighted search Info
SI = search_information(W_dir, L_dir);
E_si = -SI;
dir.E_si_asy = (E_si - E_si')./(E_si + E_si');

% Weighted diffusion efficiency
[~, E_dif] = diffusion_efficiency(W_dir);
dir.E_dif_asy = (E_dif - E_dif')./(E_dif + E_dif');

%%% Step 4: Compute weighted communication asymmetries for the undirected (original) connectome

% Weighted navigation efficiency
[~, ~, PL_nav] = navigate(L_und, D);
E_nav = 1./PL_nav;
und.E_nav_asy = (E_nav - E_nav')./(E_nav + E_nav');

% Weighted search Info
SI = search_information(W_und, L_und);
E_si = -SI;
und.E_si_asy = (E_si - E_si')./(E_si + E_si');

% Weighted diffusion efficiency
[~, E_dif] = diffusion_efficiency(W_und);
und.E_dif_asy = (E_dif - E_dif')./(E_dif + E_dif');

%%% Step 5: Determine unconnected node pairs in the symmetrized connectome
n = length(W_und);
aux = ones(n);
bux = zeros(n);
bux(triu(aux,1) == 1 & W_und == 0) = 1;
unc_idx = find(bux);

%%% Step 5: Correlation between directed and undirected asymmetries for unconneted node pairs
[r,p] = corr(und.E_dif_asy(unc_idx), dir.E_dif_asy(unc_idx));

figure;
scatter(und.E_dif_asy(unc_idx), dir.E_dif_asy(unc_idx));
lsline;
xlabel('Undirected diffusion asymmetry');
ylabel('Directed diffusion asymmetry');
title(sprintf(['Fig 6f: Directed vs undirected weighted diffusion asymmetry in the\n', ... 
    'macaque connectome: rho=%.2f, p=%.0e'],r,p));
set(gca, 'FontSize', 14);

display('Scatter plot generated.');
fprintf('\n');




