function plotModes2D(x, y, fln, T, P, Tx, Xi, varargin)
% PLOTMODES2D is a helper function that demonstrates how to index into and
% plot 2D TOD modes

if nargin > 7
    axlims = varargin{1};
else
    axlims = [min(x(1, :)) max(x(1, :)) min(y(:, 1)) max(y(:, 1))];
end
xmini = findnearest(axlims(1), x(1, :));
xmaxi = findnearest(axlims(2), x(1, :));
contourlevels = linspace(0.1, 0.9, 6);
colorGray = 0.5 * [1 1 1];
wavyLineLength = 10;

fl = fln(1);
fn = fln(2);

tiledlayout(4, 4, 'tilespacing', 'tight', 'padding', 'tight')

% catalyst mode
nexttile(1, [2 1])
cmode = real(squeeze(Xi(2, :, :)));
cmode = cmode / max(abs(cmode(:)));
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), contourlevels, 'k-'), hold on
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), -contourlevels, 'k-.'), hold off
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder
lineStart = [axlims(2) 0];
lineEnd = lineStart + wavyLineLength * [1 -0.5];
lineMid = (lineStart + lineEnd) / 2;
drawWavyLine(lineStart, lineEnd, fn - fl, 0, false, colorGray);
text(lineMid(1), lineMid(2) + 2, ['$n-l=' num2str(fn - fl) '$'], 'rot', -45)
text(-3, 0, '$\mbox{\boldmath$\hat{\xi}$}_{n-l}$')

% donor mode
nexttile(9, [2 1])
cmode = real(squeeze(Xi(1, :, :)));
cmode = cmode / max(abs(cmode(:)));
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), contourlevels, 'k-'), hold on
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), -contourlevels, 'k-.'), hold off
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder
if real(T) > 0
    lineStart = [axlims(2) 0];
    lineEnd = lineStart + wavyLineLength * [1 0.5];
    wavyLineColor = 'r';
else
    lineEnd = [axlims(2) 0];
    lineStart = lineEnd + wavyLineLength * [1 0.5];
    wavyLineColor = 'b';
end
lineMid = (lineStart + lineEnd) / 2;
drawWavyLine(lineStart, lineEnd, fl, 0, true, wavyLineColor);
text(lineMid(1), lineMid(2) + 2, ['$' num2str(real(T), 2) '$'], 'col', wavyLineColor, 'rot', 45)
text(lineMid(1), lineMid(2) - 2, ['$l=' num2str(fl) '$'], 'rot', 45)
text(-3, 0, '$\mbox{\boldmath$\hat{\xi}$}_{l}$')

% recipient mode
nexttile(8, [2 1])
cmode = real(squeeze(P(2, :, :)));
cmode = cmode / max(abs(cmode(:)));
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), contourlevels, 'k-'), hold on
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), -contourlevels, 'k-.'), hold off
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder
text(-3, 0, '$\mbox{\boldmath$\hat{\phi}$}_{n}$')

% modal transfer field
nexttile(6, [2 1])
Tx = Tx / max(abs(Tx(:)));
pcolor(x, y, Tx), shading interp, axis equal tight
colormap(bluered)
clim(max(abs(clim)) * [-1 1] / 2)

% convective mode
cmode = real(squeeze(P(1, :, :)));
cmode = cmode / max(abs(cmode(:)));
hold on, contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), contourlevels, 'k-')
contour(x(:, xmini:xmaxi), y(:, xmini:xmaxi), cmode(:, xmini:xmaxi), -contourlevels, 'k-.'), hold off
axis equal tight, axis(axlims), axis off
drawCylinder
if real(T) > 0
    lineStart = [axlims(2) 0];
    lineEnd = lineStart + wavyLineLength * [1 0];
    wavyLineColor = 'r';
else
    lineEnd = [axlims(2) 0];
    lineStart = lineEnd + wavyLineLength * [1 0];
    wavyLineColor = 'b';
end
lineMid = (lineStart + lineEnd) / 2;
drawWavyLine(lineStart, lineEnd, fn, 0, true, wavyLineColor);
text(lineMid(1), lineMid(2) + 2, ['$' num2str(real(T), 2) '$'], 'col', wavyLineColor)
text(lineMid(1), lineMid(2) - 2, ['$n=' num2str(fn) '$'])
text(-3, 0, '$\mbox{\boldmath$\hat{\psi}$}_{l\rightarrow n}$')

set(gcf, 'color', 'w')
set(findall(gcf, '-property', 'interpreter'), 'interpreter', 'latex')

end

function drawCylinder
rectangle('pos', [-0.5 -0.5 1 1], 'cur', [1 1], 'facecol', 0.9 * [1 1 1], 'edgecol', 'no')
end