%% PlotRTL1090
%  3D visualization of air traffic through RTL-SDR (dump1090) and MATLAB
%  Copyright (C) 2014  Jorge Garcia Tiscar
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 3 of the License, or
%  (at your option) any later version (see LICENSE).

%% Initialize
clear all
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
    data = urlread('http://trckr.ch:8888/persist/data/aircraft.json');
%     data   = urlread('http://192.168.1.81/persist/data/aircraft.json');
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
data = importdata('DHM200.xyz');

% convert x y z rows to z(x,y) matrix
X = data(:,1); Y = data(:,2); Z = data(:,3);
Xs = unique(X);
Ys = unique(Y);
Xi = arrayfun( @(x) find(Xs==x), X );
Yi = arrayfun( @(y) find(Ys==y), Y );
Li = Yi + (Xi-1) * numel(Ys);
XYZ = nan(numel(Ys), numel(Xs));
XYZ( Li ) = Z;

% TODO: these offsets require some tweaking
Xd = Xs-2.3e5;
Yd = Ys+5e6;
XYZd = XYZ + 0;

%% Render the recorded data
% Prepare figure
load   coords %#ok<*UNRCH>
close  all
opengl software % Hardware OpenGL rendering is unreliable for exporting images
figure('Renderer','opengl',...
       'DefaultTextFontName', 'Miryad Pro', ...
       'DefaultTextFontSize', 10,...
       'DefaultTextHorizontalAlignment', 'right')

% Settings
centerLoc = [47.1, 7.6]; % lat / long
lineColor = [0.7,0.7,0.7];
vertiExag = 5; % vertical exaggeration
markSize  = 21;
maxRadius = 200e3; % bullseye max radius (m)
hold on

mesh(Xd, Yd, XYZ);
light

% Prepare UTM scenario
mstruct      = defaultm('utm');
mstruct.zone = utmzone(centerLoc(1),centerLoc(2));
mstruct      = defaultm(mstruct);

% Plot land contours
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
[x,y]     = mfwdtran(mstruct,[countries.Lat provinces.Lat],[countries.Lon provinces.Lon]);
[xc,yc]   = mfwdtran(mstruct,centerLoc(1),centerLoc(2));
plot(x,y,'-k')

% Plot bullseye
t = linspace(0,2*pi);
for i = 0:50e3:maxRadius
    plot(xc+i.*cos(t),yc+i.*sin(t),'-','color',lineColor)
    if i>0
        text(xc+i-15e3, yc-15e3, [num2str(i/1e3) 'km'],'color',lineColor,'HorizontalAlignment','center');
    end
end

% Plot lines
plot([xc xc],[yc-i-10e3 yc+i+10e3],'-','color',lineColor)
plot([xc-i-10e3 xc+i+10e3],[yc yc],'-','color',lineColor)

% Plot traces
filter    = (alt < 150000); % You can filter per altitude (feet), speed, etc.
                          % To filter per airline use a cell function:
                          % For Ryanair: filter = cellfun(@(x) ~isempty(regexpi(x,'^RYR.*$')),flg);
color     = alt; % Color by altitude or by speed, squawk... 
[x,y]     = mfwdtran(mstruct,lat(filter),lon(filter));
scatter3(x,y,vertiExag.*alt(filter).*0.3048,markSize,color(filter),'Marker','.')
colormap cool            % change color map to e.g. "jet" or "spring"

% Some more figure settings
axis equal
extra     = 40e3; % some extra space around the bullseye
axis([xc-i-extra xc+i+extra yc-i-extra yc+i+extra])
set(gcf, 'Color', 'white');
axis off

% Prepare 3D view
pos       = [0, 0, 1920+2, floor((1080+2)/0.8)]; % This gives a 854 x 480 mp4/gif whith the selected cropping
set(gcf, 'Position', pos);
axis vis3d
view(-7,10) % Some fiddling required!
camzoom(1.4) % Same here!
set(gca, 'LooseInset', [0,0,0,0]);

% Prepare animation
filename  = 'plot1090_3D';
duration  = 7; % Adjust at will!
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
  frame = getframe(gcf,pos-[0 0 0 round(0.3*pos(4))]); % Cropping as needed (top 20% here)
  % [left, bottom, width, height]
  
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

% This is the end...
