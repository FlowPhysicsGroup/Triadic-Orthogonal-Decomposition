function [cmap] = bluered(N)

rgb = [linspace(0,1,32).' linspace(0,1,32).' linspace(1,1,32).'
       linspace(1,1,32).' linspace(1,0,32).' linspace(1,0,32).'];

if nargin==0
    cmap   = rgb;
else
    x_old       = linspace(0,1,length(rgb));
    x_new       = linspace(0,1,N);
    cmap        = zeros(N,3);
    cmap(:,1)   = interp1(x_old,rgb(:,1),x_new); 
    cmap(:,2)   = interp1(x_old,rgb(:,2),x_new);
    cmap(:,3)   = interp1(x_old,rgb(:,3),x_new);
end

