function plotGeneralObservable(x, y, P, varargin)
% PLOTGENERALOBSERVABLE plots TOD modes from non-velocity data

if nargin > 3
    axlims = varargin{1};
else
    axlims = [min(x(1, :)) max(x(1, :)) min(y(:, 1)) max(y(:, 1))];
end

tiledlayout(2, 1)

% right singular vector
nexttile
cmode = real(squeeze(P(2, :, :)));
cmode = cmode / max(abs(cmode(:)));
pcolor(x, y, cmode), shading interp
axis equal tight, axis(axlims)
drawCylinder
text(axlims(1)-2, 0, '$\mbox{\boldmath$\hat{\phi}$}_{n}$')

% left singular vector
nexttile
cmode = real(squeeze(P(1, :, :)));
cmode = cmode / max(abs(cmode(:)));
pcolor(x, y, cmode), shading interp
axis equal tight, axis(axlims)
drawCylinder
text(axlims(1)-2, 0, '$\mbox{\boldmath$\hat{\psi}$}_{(n-l)\circ l}$')

set(gcf, 'color', 'w')
set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'latex')

end

function drawCylinder
rectangle('pos', [-0.5 -0.5 1 1], 'cur', [1 1], 'facecol', 0.9 * [1 1 1], 'edgecol', 'no')
end