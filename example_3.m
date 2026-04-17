%% EXAMPLE 3: Apply TOD in the absence of velocity data (e.g. pressure or schlieren data)
%  In this example, we demonstrate how to analyze datasets where velocity
%  fields are unavailable (see equation 4.1 in [1]). Modal energy flow
%  analysis requires velocities and their gradients. In their absence, TOD
%  can still detect triads based on statistical correlation. In this
%  setting, physical interpretation of the results is analogous to
%  bispectral mode decomposition (BMD; [2]) and the codes are used
%  similarly. For demonstration, we apply TOD to vorticity fields computed
%  from the same dataset [3] as in example 1.
%
%   References:
%     [1] Yeung, B., Chu, T., and Schmidt, O. T., Triadic orthogonal
%         decomposition reveals nonlinearity in fluid flows, J. Fluid Mech. 1031, A34, 2026
%     [2] Schmidt, O. T., Bispectral mode decomposition of nonlinear flows, 
%         Nonlinear Dyn. 102, 2479-2501, 2020
%     [3] Chu, T. and Schmidt, O. T., RBF-FD discretization of the Navier-Stokes
%         equations on scattered but staggered nodes, J. Comp. Phys. 474, 111756, 2023
%
%  B. Yeung (byeung@ucsd.edu), T. Chu, and O. T. Schmidt
%  Original script: BY
%  Last revision:   16-Apr-2026 (BY)

% Clean up the worspace.
clc, clear
addpath('aux')

% Specify estimation parameters
nDFT = 150;
nOvlp = ceil(nDFT / 2);
window = ones(nDFT, 1);
weight = [];

%% Handing general observables
load('cylinder_Re100.mat')

dt = data{1};
x = data{2};
y = data{3};
u = data{4}; % size(u) = [nt ny nx]
v = data{5};

nx = size(x, 2);
ny = size(x, 1);
nt = size(u, 1);
dx = abs(x(1, 2) - x(1, 1));
dy = abs(y(2, 1) - y(1, 1));

% Build differential operators (only used to compute vorticity, not needed otherwise)
Dx_1D = Dmats_SBP(nx, dx, 4);
Dy_1D = Dmats_SBP(ny, dy, 4);
Dx = kron(Dx_1D, speye(ny));
Dy = sparse(kron(speye(nx), Dy_1D));

% Calculate vorticity
vort = Dx * double(v(:, :).') - Dy * double(u(:, :).');
X = single(reshape(vort.', [nt 1 ny nx])); % size(X) = [nt nvar=1 ny nx]
X = repmat(X, 1, 2, 1, 1); % inflate nvar dimension to prevent it from being squeezed out
clear Dx_1D Dy_1D Dx Dy

% Run TOD
[B, P, f] = tod(X, window, weight, nOvlp, dt);

%% Plot bispectra
% Normalize frequencies by vortex-shedding frequency
f0 = 0.16728;
[fl, fn] = ndgrid(f);
fl = fl / f0;
fn = fn / f0;

modei = 1; % modal rank to plot

figure

% Mode bispectrum
Bi = B(:, :, modei);
eps_display = 5e-2;
B_idx = Bi > eps_display; % threshold for ease of plotting
axis(5 * [-1 1 0 1]), axlims = axis;
clb = plotBispectrum(fl(B_idx), fn(B_idx), Bi(B_idx), axlims);
clim(200 * [0 1])
colormap(gca, flipud(gray))
title('mode bispectrum (singular values)')
ylabel(clb, ['$\sigma_' num2str(modei) '$'], 'interpreter', 'latex')

set(gcf, 'color', 'w')
set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'latex')

%% Plot modes
% Specify which donor and recipient frequency indices to plot
li = 86;
ni = 106;

% Retrieve modes
P_this = P{li, ni}(:, :, :, 1, modei); % left and right singular vectors

% Zoom in to region of interest
axlims = [-2 10 -2.5 2.5];

figure
plotGeneralObservable(x, y, P_this, axlims)
