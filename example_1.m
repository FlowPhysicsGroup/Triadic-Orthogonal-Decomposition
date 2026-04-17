%% EXAMPLE 1: Load data, build differential operators, and calculate and plot TOD.
%  This example corresponds to the cylinder flow analysis by Yeung, Chu, &
%  Schmidt [1] using TOD and reproduces figures 4 and 8. The Re=100
%  cylinder wake dataset was calculated using the RBF DNS solver developed
%  by Chu & Schmidt [2].
%
%   References:
%     [1] Yeung, B., Chu, T., and Schmidt, O. T., Triadic orthogonal
%         decomposition reveals nonlinearity in fluid flows, J. Fluid Mech. 1031, A34, 2026
%     [2] Chu, T. and Schmidt, O. T., RBF-FD discretization of the
%         Navier-Stokes equations on scattered but staggered nodes, J. Comp. Phys. 474, 111756, 2023
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

%% Load data
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

X = permute(cat(4, u, v), [1 4 2 3]); % size(X) = [nt nvar ny nx]

% Build differential operators
Dx_1D = Dmats_SBP(nx, dx, 4);
Dy_1D = Dmats_SBP(ny, dy, 4);
Dx = kron(Dx_1D, speye(ny));
Dy = sparse(kron(speye(nx), Dy_1D));

clear opts
% Provide form of quadratic nonlinearity
% in this case: Q([u1,v1],[u2,v2]) = -[u1*du2/dx + v1*du2/dy
%                                      u1*dv2/dx + v1*dv2/dy]
% As an alternative to calculating velocity gradients in TOD, gradients can
% be precomputed and exported from a numerical solver, then appended to the
% data matrix X, e.g. X(:, 3, :, :) = du/dx, etc. In this case opts.LHS
% also need to be specified because the LHS of the NSE includes only
% velocities, not their gradients
opts.Q = @(q1, q2) -[squeeze(q1(1, :, :)) .* (Dx * double(squeeze(q2(1, :, :)))) + squeeze(q1(2, :, :)) .* (Dy * double(squeeze(q2(1, :, :))));...
                     squeeze(q1(1, :, :)) .* (Dx * double(squeeze(q2(2, :, :)))) + squeeze(q1(2, :, :)) .* (Dy * double(squeeze(q2(2, :, :))))];
opts.nmode = 2; % store two leading modes

% Run TOD
[B, P, f, T, A, Xi] = tod(X, window, weight, nOvlp, dt, opts);

%% Plot bispectra
% Normalize frequencies by vortex-shedding frequency
f0 = 0.16728;
[fl, fn] = ndgrid(f);
fl = fl / f0;
fn = fn / f0;

modei = 1; % modal rank to plot

figure
tiledlayout(2, 1);

% Mode bispectrum
nexttile
Bi = B(:, :, modei);
eps_display = 5e-2;
B_idx = Bi > eps_display; % threshold for ease of plotting
axis(5 * [-1 1 0 1]), axlims = axis;
clb = plotBispectrum(fl(B_idx), fn(B_idx), Bi(B_idx), axlims);
clim(20 * [0 1])
colormap(gca, flipud(gray))
title('mode bispectrum (singular values)')
ylabel(clb, ['$\sigma_' num2str(modei) '$'], 'interpreter', 'latex')

% Modal energy budget
nexttile
realT = real(T(:, :, modei));
eps_display = 1e-3;
T_idx = abs(realT) > eps_display;
clb = plotBispectrum(fl(T_idx), fn(T_idx), realT(T_idx), axlims);
clim([-1 1])
colormap(gca, bluered)
title('modal energy budget')
hold on
plot([0 0], [0 axlims(4)], 'g') % production
plot([0 axlims(4)], [0 axlims(4)], 'm') % linear advection
hold off
ylabel(clb, ['$\hat{\mathcal{T}}_{l\to n,' num2str(modei) '}^{\mathrm{avg},\mathcal{R}}$'], 'interpreter', 'latex')

set(gcf, 'color', 'w')
set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'latex')

%% Plot modes
% Specify which donor and recipient frequency indices and variable index to plot
li = 86;
ni = 106;
fl = f(li) / f0;
fn = f(ni) / f0;
vari = 1;

% Retrieve a bunch of modes
B_this = B(li, ni, modei); % singular value
T_this = T(li, ni, modei); % modal transfer
P_this = P{li, ni}(:, :, :, vari, modei); % convective and recipient modes
Tx_this = B_this * real(squeeze(sum(conj(P{li, ni}(2, :, :, :, modei)) .* P{li, ni}(1, :, :, :, modei), 4))); % modal transfer field
Xi_this = Xi{li, ni}(:, :, :, vari, modei); % donor and catalyst modes

% Zoom in to region of interest
axlims = [-2 10 -2.5 2.5];

figure
plotModes2D(x, y, [fl, fn], T_this, P_this, Tx_this, Xi_this, axlims)
