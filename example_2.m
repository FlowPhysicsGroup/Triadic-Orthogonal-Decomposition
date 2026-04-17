%% EXAMPLE 2: Apply TOD to 3D data.
%  This example demonstrates the application of TOD [1] to 3D data. The
%  dataset [2] is the same as in example 1. For demonstration only, we make
%  the 2D velocity fields artifically 3D by extruding them in z and
%  appending w=0. Nothing fundamentally has changed in the data, so we
%  should recover the results from example 1 (excepting scaling
%  differences). Could take a while to run.
%
%  This script can be adapted to your own 3D flow data by arranging the
%  velocity components u, v, and w into a data matrix X of size 
%  [nt nvar nz ny nx], such that X(:, 1, :, :, :) = u, X(:, 2, :, :, :) = v,
%  and X(:, 3, :, :, :) = w. Usage is similar for data in cylindrical or
%  spherical coordinates.
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

%% Handling 3D data
load('cylinder_Re100.mat')

dt = data{1};
x = data{2};
y = data{3};
u = data{4};
v = data{5};

nx = size(x, 2);
ny = size(x, 1);
nt = size(u, 1);
dx = abs(x(1, 2) - x(1, 1));
dy = abs(y(2, 1) - y(1, 1));

% Make 3D data
nz = 10;
dz = 0.1;
x = x(1, :);
y = y(:, 1);
z = 0:dz:dz*(nz-1);
u = permute(repmat(u, 1, 1, 1, nz), [1 4 2 3]); % size(u) = [nt nz ny nx]
v = permute(repmat(v, 1, 1, 1, nz), [1 4 2 3]);
w = zeros(size(u), 'like', u);

X = permute(cat(5, u, v, w), [1 5 2 3 4]); % size(X) = [nt nvar nz ny nx]

% Build differential operators
Dx_1D = Dmats_SBP(nx, dx, 4);
Dy_1D = Dmats_SBP(ny, dy, 4);
Dz_1D = Dmats_SBP(nz, dz, 4);
Dx = kron(kron(Dx_1D, speye(ny)), speye(nz));
Dy = kron(kron(speye(nx), Dy_1D), speye(nz));
Dz = kron(kron(speye(nx), speye(ny)), Dz_1D);

clear opts
% Provide form of quadratic nonlinearity
% in this case: Q([u1,v1,w1],[u2,v2,w2]) = -[u1*du2/dx + v1*du2/dy + w1*du2/dz
%                                            u1*dv2/dx + v1*dv2/dy + w1*dv2/dz
%                                            u1*dw2/dx + v1*dw2/dy + w1*dw2/dz]
% As an alternative to calculating velocity gradients in TOD, gradients can
% be precomputed and exported from a numerical solver, then appended to the
% data matrix X, e.g. X(:, 4, :, :, :) = du/dx, etc. In this case opts.LHS
% also need to be specified because the LHS of the NSE includes only
% velocities, not their gradients
opts.Q = @(q1, q2) -[squeeze(q1(1, :, :)) .* (Dx * double(squeeze(q2(1, :, :)))) + squeeze(q1(2, :, :)) .* (Dy * double(squeeze(q2(1, :, :)))) + squeeze(q1(3, :, :)) .* (Dz * double(squeeze(q2(1, :, :))));...
                     squeeze(q1(1, :, :)) .* (Dx * double(squeeze(q2(2, :, :)))) + squeeze(q1(2, :, :)) .* (Dy * double(squeeze(q2(2, :, :)))) + squeeze(q1(3, :, :)) .* (Dz * double(squeeze(q2(2, :, :))));...
                     squeeze(q1(1, :, :)) .* (Dx * double(squeeze(q2(3, :, :)))) + squeeze(q1(2, :, :)) .* (Dy * double(squeeze(q2(3, :, :)))) + squeeze(q1(3, :, :)) .* (Dz * double(squeeze(q2(3, :, :))))];
opts.nmode = 1; % store only the leading mode
opts.nfreq = 35; % store 35 frequencies

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
axis(3.5 * [-1 1 0 1]), axlims = axis;
clb = plotBispectrum(fl(B_idx), fn(B_idx), Bi(B_idx), axlims);
clim(200 * [0 1])
colormap(gca, flipud(gray))
title('mode bispectrum (singular values)')
ylabel(clb, ['$\sigma_' num2str(modei) '$'], 'interpreter', 'latex')

% Modal energy budget
nexttile
realT = real(T(:, :, modei));
eps_display = 1e-3;
T_idx = abs(realT) > eps_display;
clb = plotBispectrum(fl(T_idx), fn(T_idx), realT(T_idx), axlims);
clim(10 * [-1 1])
colormap(gca, bluered)
title('modal energy budget')
hold on
plot([0 0], [0 axlims(4)], 'g') % production
plot([0 axlims(4)], [0 axlims(4)], 'm') % linear advection
hold off
ylabel(clb, ['$\hat{\mathcal{T}}_{l\to n,' num2str(modei) '}^{\mathrm{avg},\mathcal{R}}$'], 'interpreter', 'latex')

set(gcf, 'color', 'w')
set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'latex')

%% Plot 3D modes
% Specify which donor and recipient frequency indices and variable index to plot
li = 86;
ni = 106;
fl = f(li) / f0;
fn = f(ni) / f0;
vari = 1;

% Retrieve a bunch of modes
B_this = B(li, ni, modei); % singular value
T_this = T(li, ni, modei); % modal transfer
P_this = P{li, ni}(:, :, :, :, vari, modei); % convective and recipient modes
Tx_this = B_this * real(squeeze(sum(conj(P{li, ni}(2, :, :, :, :, modei)) .* P{li, ni}(1, :, :, :, :, modei), 5))); % modal transfer field
Xi_this = Xi{li, ni}(:, :, :, :, vari, modei); % donor and catalyst modes

% Zoom in to region of interest
axlims = [-2 10 -2.5 2.5 0 0.9];

figure
plotModes3D(x, y, z, [fl, fn], T_this, P_this, Tx_this, Xi_this, axlims)
