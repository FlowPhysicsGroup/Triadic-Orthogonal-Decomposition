function clb = plotBispectrum(fl, fn, B, varargin)
% PLOTBISPECTRUM is a helper function that demonstrates how to plot a TOD
% mode bispectrum (i.e., bispectrum of singular values) or modal energy
% budget

scatter(fl, fn, 100, B, 'filled')
axis equal tight
if nargin > 3
    axlims = varargin{1};
else
    axlims = [min(fl) max(fl) min(fn) max(fn)];
end
axis(axlims)
xlabel('$f_l/f_0$'), ylabel('$f_n/f_0$')
hold on, plot([0 -axlims(4)], [0 axlims(4)], 'k'), hold off
text(-axlims(4) * 0.75, axlims(4) * 0.75, '$f_k/f_0$', 'rot', -45, 'ver', 'bot')
set(gca, 'box', 'on', 'layer', 'top')
clb = colorbar('location', 'eastoutside'); % return colorbar handle

end