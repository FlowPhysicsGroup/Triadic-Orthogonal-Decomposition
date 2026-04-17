function obj = drawWavyLine(p1, p2, frequency, phase, varargin)
% DRAWWAVYLINE draws a sinusoidal line between two points
%   drawWavyLine([x1, y1], [x2, y2])
    drawArrow = true;
    col = 'k';
    linestyle = '-';
    linewidth = 0.5;
    if nargin > 4
        drawArrow = varargin{1};
        if nargin > 5
            col = varargin{2};
            if nargin > 6
                linestyle = varargin{3};
                if nargin > 7
                    linewidth = varargin{4};
                end
            end
        end
    end

    % Number of points along the line
    n = 100;
    
    % Line vector
    t = linspace(0, 1, n);
    x = p1(1) + (p2(1) - p1(1)) * t;
    y = p1(2) + (p2(2) - p1(2)) * t;

    % Direction vector and perpendicular unit vector
    v = p2 - p1;
    v = v / norm(v);
    nvec = [-v(2), v(1)];

    % Add sinusoidal displacement
    amplitude = 1;   % relative amplitude
    y = y + amplitude * sin(2 * pi * (frequency * t + phase)) * nvec(2);
    x = x + amplitude * sin(2 * pi * (frequency * t + phase)) * nvec(1);

    % Plot
    hold on
    obj = plot(x, y, 'color', col, 'linestyle', linestyle, 'linewidth', linewidth);
    hold off
    set(obj, 'clipping', 'off')

    % Arrowhead at midpoint
    mid = round(n/2);
    dir = [x(mid+1)-x(mid-1), y(mid+1)-y(mid-1)];
    dir = dir / norm(dir);
    perp = [-dir(2), dir(1)];

    L = 0.3;  % arrow length (scaled by linewidth)
    W = L * 2;  % arrow width

    % Triangle coordinates
    tip = [x(mid), y(mid)] + L * dir;
    base = [x(mid), y(mid)] - L * dir;
    p1tri = base + 0.5 * W * perp;
    p2tri = base - 0.5 * W * perp;
    
    if drawArrow
        hold on
        fill([tip(1) p1tri(1) p2tri(1)], [tip(2) p1tri(2) p2tri(2)], ...
            col, 'EdgeColor', 'none', 'Clipping', 'off');
        hold off
    end
end