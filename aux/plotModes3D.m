function plotModes3D(x, y, z, fln, T, P, Tx, Xi, varargin)
% PLOTMODES3D is a helper function that demonstrates how to index into and
% plot 3D TOD modes as isosurfaces

if nargin > 8
    axlims = varargin{1};
else
    axlims = [min(x) max(x) min(y) max(y) min(z) max(z)];
end
xmini = findnearest(axlims(1), x);
xmaxi = findnearest(axlims(2), x);
isoval = 0.1 * [-1 1];
colorGray = 0.5 * [1 1 1];
facecolPos = colorGray;
facecolNeg = 0.875 * [1 1 1];
wavyLineLength = 10;

fl = fln(1);
fn = fln(2);

tiledlayout(4, 3, 'tilespacing', 'tight', 'padding', 'tight')

% catalyst mode
nexttile(1, [2 1])
cmode = permute(real(squeeze(Xi(2, :, :, :))), [1 3 2]);
cmode = cmode / max(abs(cmode(:)));
for i = 1:length(isoval)
    if isoval(i) > 0
        enclose = 'above';
        facecol = facecolPos;
    else
        enclose = 'below';
        facecol = facecolNeg;
    end
    patch(isosurface(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i)), 'facecol', facecol, 'edgecol', 'no'),
    patch(isocaps(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i), enclose), 'facecol', facecol, 'edgecol', 'no')
end
view(3), camlight left, camlight(45, -15), camlight(45, -15), lighting gouraud
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder(z(end) - z(1))
lineStart = [axlims(2) 0];
lineEnd = lineStart + wavyLineLength * [0 -1.1];
lineMid = (lineStart + lineEnd) / 2;
drawWavyLine(lineStart, lineEnd, fn - fl, 0, false, colorGray);
text(lineMid(1), lineMid(2) + 2, ['$n-l=' num2str(fn - fl) '$'], 'rot', -45)
text(-3, 0, '$\mbox{\boldmath$\hat{\xi}$}_{n-l}$')

% donor mode
nexttile(7, [2 1])
cmode = permute(real(squeeze(Xi(1, :, :, :))), [1 3 2]);
cmode = cmode / max(abs(cmode(:)));
for i = 1:length(isoval)
    if isoval(i) > 0
        enclose = 'above';
        facecol = facecolPos;
    else
        enclose = 'below';
        facecol = facecolNeg;
    end
    patch(isosurface(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i)), 'facecol', facecol, 'edgecol', 'no'),
    patch(isocaps(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i), enclose), 'facecol', facecol, 'edgecol', 'no')
end
view(3), camlight left, camlight(45, -15), camlight(45, -15), lighting gouraud
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder(z(end) - z(1))
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
nexttile(6, [2 1])
cmode = permute(real(squeeze(P(2, :, :, :))), [1 3 2]);
cmode = cmode / max(abs(cmode(:)));
for i = 1:length(isoval)
    if isoval(i) > 0
        enclose = 'above';
        facecol = facecolPos;
    else
        enclose = 'below';
        facecol = facecolNeg;
    end
    patch(isosurface(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i)), 'facecol', facecol, 'edgecol', 'no'),
    patch(isocaps(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i), enclose), 'facecol', facecol, 'edgecol', 'no')
end
view(3), camlight left, camlight(45, -15), camlight(45, -15), lighting gouraud
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder(z(end) - z(1))
text(-3, 0, '$\mbox{\boldmath$\hat{\phi}$}_{n}$')

% convective mode colored by modal transfer field
nexttile(5, [2 1])
Tx = permute(Tx / max(abs(Tx(:))), [1 3 2]);
cmode = permute(real(squeeze(P(1, :, :, :))), [1 3 2]);
cmode = cmode / max(abs(cmode(:)));
for i = 1:length(isoval)
    if isoval(i) > 0
        enclose = 'above';
    else
        enclose = 'below';
    end
    patch(isosurface(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i), Tx(:, xmini:xmaxi, :)), 'facecol', 'interp', 'edgecol', 'no'),
    fvc = isocaps(x(xmini:xmaxi), z, y, cmode(:, xmini:xmaxi, :), isoval(i), enclose);
    TxIsocaps = interp3(x(xmini:xmaxi), z, y, Tx(:, xmini:xmaxi, :), fvc.vertices(:, 1), fvc.vertices(:, 2), fvc.vertices(:, 3));
    patch(fvc, 'FaceVertexCData', TxIsocaps, 'facecol', 'interp', 'edgecol', 'no')
end
colormap(bluered)
clim(max(abs(clim)) * [-1 1] / 2)
view(3), camlight left, camlight(45, -15), camlight(45, -15), lighting gouraud
set(gca, 'clipping', 'off')
axis equal tight, axis(axlims), axis off
drawCylinder(z(end) - z(1))
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

function drawCylinder(varargin)
cylinderHeight = 1;
if nargin > 0
    cylinderHeight = varargin{1};
end
radius = 0.5;
circumference = 64;
[X, Y, Z] = cylinder(radius, circumference);
Z = Z * cylinderHeight;

% draw hollow cylinder
cylinderCol = 0.35 * [1 1 1];
hold on, surf(X, Z, Y, 'facecol', cylinderCol, 'edgecol', 'no')

% draw cylinder caps
for i = 1:2
    fill3(X(i, :), Z(i, :), Y(i, :), cylinderCol, 'edgecol', 'no')
end
hold off

end