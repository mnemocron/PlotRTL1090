%% PlotRTL1090
%  3D visualization of air traffic through RTL-SDR (dump1090) and MATLAB
%  Copyright (C) 2014  Jorge Garcia Tiscar
%  Modified for Switzerland DEM 2020 github.com/mnemocron
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 3 of the License, or
%  (at your option) any later version (see LICENSE).

%% Initialize
clear all; clc;
if exist('coords.mat','file')
    load coords
else
    lat = []; % latitude
    lon = []; % longitude 
    alt = []; % altitude
    spd = []; % speed
    flg = {}; % flight
    tim = []; % time
end

%% Acquisition loop
while true % adjust duration here or stop with Ctrl+C

    % Get data from server
    data   = urlread('http://192.168.1.81/tar1090/data/aircraft.json');
    planes = fromjson(data); % https://github.com/christianpanton/matlab-json
    
    % Parse data
    N   = length(planes.aircraft);
    val = 0;
    for i = 1:length(planes.aircraft)
        plane = planes.aircraft{i};
        % changed sanity check to "isfield"
        if isfield(plane,'lat')
            if isfield(plane,'alt_geom')
                if isfield(plane,'mach')
                    if isfield(plane,'flight')
                        val = val + 1;
                        lat = [lat plane.lat];
                        lon = [lon plane.lon];
                        alt = [alt plane.alt_geom];
                        spd = [spd plane.mach];
                        tim = [tim now];
                        flg{length(lat)} = plane.flight;
                    end
                end
            end
        end
    end
    disp([num2str(N) ' planes detected with ' num2str(val) ' valid coords'])
    
    % Save and continue
    save('coords.mat','lat','lon','alt','spd','flg','tim')
    pause(1)
end

%---------------------------------------------------- After stopping the loop execute the following:
%%
% import Digital Elevation Model (DEM)
% takes ~1 minute
% 200m free DEM dataset from:
% https://shop.swisstopo.admin.ch/en/products/height_models/dhm25200
data = importdata('DHM200.xyz');

% convert x y z rows to z(x,y) matrix
X = data(:,1); Y = data(:,2); Z = data(:,3);
Xs = unique(X);
Ys = unique(Y);
Xi = arrayfun( @(x) find(Xs==x), X );
Yi = arrayfun( @(y) find(Ys==y), Y );
Li = Yi + (Xi-1) * numel(Ys);
XYZ = nan(numel(Ys), numel(Xs));  % empty NaN matrix
XYZ( Li ) = Z;  % fill with Z valzes

% altitude feet to meter
altm = alt./3.281;

%%
% Settings
% calculation requires ~1 minute
centerLoc = [47.1, 7.6]; % lat / long
lineColor = [0.7,0.7,0.7];
vertiExag = 5; % vertical exaggeration
markSize  = 21;
maxRadius = 200e3; % bullseye max radius (m)

% Prepare UTM scenario
mstruct      = defaultm('utm');
mstruct.zone = utmzone(centerLoc(1),centerLoc(2));
mstruct      = defaultm(mstruct);

% SHPs from Natural Earth: http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
SHPdir    = '.\SHPs\';
% Change 'es' for the ISO_A2 code of your country
% For all countries delete selector or input non-existent field (ISO_A2 => foo)
countries = shaperead([SHPdir 'ne_10m_admin_0_countries.shp'],...
            'Selector',{@(x) strcmpi(x,'CHE'),'ISO_A3'},'UseGeoCoords', true);
% Change 'ES.VC' for the provinces/states of your preference or use a RegExp
% for all provinces: @(x) strcmpi(x,'ES.VC') => @(x) ~isempty(regexpi(x,'^ES.*$'))
provinces = shaperead([SHPdir 'ne_10m_admin_1_states_provinces.shp'],...
            'Selector',{@(x) ~isempty(regexpi(x,'^CH.*$')),'code_hasc'},'UseGeoCoords', true);
[x,y]     = mfwdtran(mstruct,[countries.Lat],[countries.Lon]); % exclude provinces
x = [x, x(1)]; % close the polygon shape
y = [y, y(1)];

% TODO: these offsets require some tweaking - still not 100% aligned
Xd = Xs-2.2e5;
Yd = Ys+5e6 -0.3e4;
XYZd = XYZ + 0;
zoffs = 2000;

% Find only the DEM points that are inside the national borders
yq = repelem(Yd, length(Xd));    % repeat x vals for all y vals
xq = repmat(Xd', 1, length(Yd)); % repeat y vals for all x vals
in = inpolygon(xq,yq,x,y);  % pick points inside of 

dstep = xq(2)-xq(1);
XYZb = nan(numel(Ys), numel(Xs));
for k=1:length(in)
    if(in(k))
        xi = (xq(k) - xq(1))/dstep +1;
        yi = (yq(k) - yq(1))/dstep +1;
        XYZb(yi,xi) = XYZd(yi,xi);
    end
end
XYZd = XYZb;
clear XYZb;
% mesh(Xd, Yd, XYZb)

% generate a colormap for relief data
% start at blue (water) --> green --> grey (mountain) --> white (snow)
rmin = min(XYZd(:));
rmax = max(XYZd(:));
vec = [rmax;      rmax-1500; (rmax-rmin)*2/6; (rmax-rmin)/6; rmin+150;  rmin];
vec = ((double(vec)./double(max(vec)).*100 ));
hex = ['#dedede'; '#dedede'; '#a3a3a3';       '#419600';     '#419600'; '#0064b5'];
hex = (hex);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 256;
relievColors = interp1(vec,raw,linspace(100,0,N),'pchip');
relievColors(relievColors<0) = 0;
relievColors = flip(uint8(relievColors.*255));

% Test Colormap
% rgbplot(relievColors)
% hold on
% colormap(relievColors)
% colorbar('Ticks',[])

%% Render the recorded data
% Prepare figure
load   coords %#ok<*UNRCH>
close  all
opengl software % Hardware OpenGL rendering is unreliable for exporting images
figure('Renderer','opengl',...
       'DefaultTextFontName', 'Miryad Pro', ...
       'DefaultTextFontSize', 10,...
       'DefaultTextHorizontalAlignment', 'right')
hold on

mesh(Xd, Yd, XYZd);
light();
material dull;
colormap(relievColors);

% Plot land contours
% this time also include the provinces
[x,y]     = mfwdtran(mstruct,[countries.Lat provinces.Lat],[countries.Lon provinces.Lon]);
[xc,yc]   = mfwdtran(mstruct,centerLoc(1),centerLoc(2));
z = x.*0 + zoffs; % z offset of where to plot the borders

plot3(x,y,z,'-k')

% Plot bullseye
t = linspace(0,2*pi);
for i = 0:50e3:maxRadius
    plot3(xc+i.*cos(t),yc+i.*sin(t), t.*0+zoffs,'-','color',lineColor)
    if i>0
        text(xc+i-15e3, yc-15e3, zoffs, [num2str(i/1e3) 'km'],'color',lineColor,'HorizontalAlignment','center');
    end
end

% Plot lines
plot3([xc xc],[yc-i-10e3 yc+i+10e3], [zoffs zoffs],'-','color',lineColor)
plot3([xc-i-10e3 xc+i+10e3],[yc yc], [zoffs zoffs], '-','color',lineColor)

% Plot tracks
filter    = (altm < 150000); % You can filter per altitude (meter), speed, etc.
                          % To filter per airline use a cell function:
                          % For Ryanair: filter = cellfun(@(x) ~isempty(regexpi(x,'^RYR.*$')),flg);
color     = altm; % Color by altitude or by speed, squawk... 
[x,y]     = mfwdtran(mstruct,lat(filter),lon(filter));

% Make your 'colormap' based on z or whatever you want...
c_dat    = double(altm');
c_map = flip(spring(length(c_dat)));     % pick your map
trackColors = interp1(linspace(min(c_dat), max(c_dat), size(c_dat, 1)), c_map, c_dat(:));
scatter3(x,y,vertiExag.*altm(filter).*0.3048,markSize,trackColors,'Marker','.')

% Some more figure settings
axis equal
extra     = 40e3; % some extra space around the bullseye
axis([xc-i-extra xc+i+extra yc-i-extra yc+i+extra])
set(gcf, 'Color', 'white');
axis off

% Prepare 3D view
% This gives a 1920 x 1080 mp4/gif whith the selected cropping
pos       = [0, 0, 1920+2, floor((1080+2)/0.8)];
set(gcf, 'Position', pos);
axis vis3d
view(-7,26) % Some fiddling required!
camzoom(1.9) % Same here!
set(gca, 'LooseInset', [0,0,0,0]);

%%
% Prepare animation
% rendering takes ~10 minutes
tic
filename  = 'plot1090_3D';
duration  = 12; % Adjust at will!
frameRate = 25; 

% Prepare MP4
try
    myVideo = VideoWriter(filename,'MPEG-4');
catch error % Video didn't close!
    close(myVideo);
    myVideo = VideoWriter(filename,'MPEG-4');
end
myVideo.FrameRate = frameRate;  
myVideo.Quality   = 100; % High quality! 
open(myVideo)

% Rotation loop
frames = frameRate*duration; 
for i = 1:frames

  % Move camera
  camorbit(360/frames,0,'data',[0 0 1])
  % new Matlab Version does not support cropping anymore idk.
  % frame = getframe(gcf,pos-[0 0 0 round(0.3*pos(4))]); % Cropping as needed (top 20% here)
  % [left, bottom, width, height]
  frame = getframe(gcf);
  % Write GIF
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im,256);
  if i == 1;
      imwrite(imind,cm,[filename '.gif'],'gif','Loopcount',inf,'DelayTime',1/frameRate);
  else
      imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',1/frameRate);
  end
  
  % Write MP4
  writeVideo(myVideo,frame);
end
close(myVideo);
toc

% This is the end...
