function earthPlot(body)
% earthPlot Generate a plot of the Earth in the current axes
%
% Inputs:
%   body - value indicating body to plot
%          1 = Earth (Default)
%          2 = Moon


if ~exist('body', 'var')
    body = 1; % Set default body to Earth
end

% Options
npanels = 180;
alpha   = 1; % globe transparency level
GMST0 = 4.89496121282306; % Epoch corresponding to J2000

if body == 1
    
    % 2:1 unprojected Earth texture image (must be in current path)
    image_file = 'EarthSatelliteImage.jpg';
    
    % Treat the Earth as spherical
    erad = 6371.0087714; % equatorial radius [km]
    prad = erad;         % polar radius [km]

elseif body == 2
    
    % 2:1 unprojected Moon texture image (must be in current path)
    image_file = 'MoonSatelliteImage.jpg';
    
    % Treat the moon as spherical
    erad = 1738.1; % equatorial radius [km]
    prad = erad;   % polar radius [km]
    
end
    
% Turn off the normal axes
set(gca, 'NextPlot','add', 'Visible','on');

axis equal;
axis auto;

% Set initial view
view(3);

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

% Load body image for texture map
cdata = imread(image_file);

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end