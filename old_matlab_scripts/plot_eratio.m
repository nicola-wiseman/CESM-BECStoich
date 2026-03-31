%% Plot CESM e-ratio (POCexp/NPP) for VarAll and FixAll models
clearvars; close all

addpath ../data/
addpath ../utils/
file1 = ('/home/mv23682/Documents/MATLAB/data/varAll.pop.h.0281_0300.nc');
file2 = ('/home/mv23682/Documents/MATLAB/data/fixAll.pop.h.0281_0300.nc');
files = {file1,file2};

% General Variables
lat = ncread(file1,'TLAT');
lon = ncread(file1,'TLONG');
tarea = ncread(file1,'TAREA').*0.0001; % convert to meters
depth = ncread(file1,'z_t').*1e-2; % convert to meters
depth_150m = ncread(file1,'z_t_150m').*1e-2; % convert to meters
vol = reshape(tarea,[size(tarea, 1), size(tarea, 2), 1]) .* reshape(depth,[1, 1, size(depth, 1)]);

NPP = NaN([100 116 15 length(files)]);
NPP_100m = NaN([100 116 length(files)]);
POC_100m = NaN([100 116 length(files)]);
zonal_mean = NaN([180 length(files)]);

for i = 1:length(files)
    photoC_sp = ncread(files{i},'photoC_sp');
    photoC_diat = ncread(files{i},'photoC_diat');
    photoC_diaz = ncread(files{i},'photoC_diaz');
    % Unpack data
    % Carbon
    NPP(:,:,:,i) = photoC_sp + photoC_diat + photoC_diaz; % units mmol/m3/s
    NPP_100m(:,:,i) = sum(NPP(:,:,1:10,i).*1000,3); % sum over top 100m
    NPP_100m(:,:,i) = NPP_100m(:,:,i).*(365.25 * 86400.0 * 12.011 * 0.00001); % convert to molC/m2/yr
    POC_FLUX_IN = ncread(files{i},'POC_FLUX_IN'); % units mmol/m3/s
    POC_100m(:,:,i) = POC_FLUX_IN(:,:,depth==105).*(365.25 * 86400.0 * 12.011 * 0.00001); % convert to molC/m2/yr
end
e_ratio = POC_100m./NPP_100m; % calculate e-ratio (EP/NPP), where EP is POC flux @100m,...
                             % and NPP is integrated NPP in the top 100m
%% Plot E-ratio
% % Interpolate to 180 x 360 degree grid
cmap = cbrewer('seq','YlGnBu',9);
cmap = flipud(cmap);

for i = 1:length(files)
    eratio_180x360 = interpTo90x180(lon,lat,e_ratio(:,:,i),'nearest');
    load woa_grid.mat
    eratio_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
    % Create Figure
    figure(i);
    set(gcf,'Position',[500 100 1000 500],'Color','white')
    v = [0.00,0.10,0.1125,0.125,0.1375,0.15,0.1625,0.175,0.20,0.30];
    contourfnu(eratio_180x360.x,eratio_180x360.y,eratio_180x360.m,v,cmap,'eastoutside','true','pcolor')
    % add coastline
    hold on
    [~,c] = contour(grid.XT,grid.YT,M3d(:,:,1),[1 1],'k');
    c.LineWidth = 2;
    % Format figure text and axis labels
    fig = gca;
    fig.FontSize = 12;
    fig.FontWeight = 'bold';
    fig.TickDir = 'out';
    fig.TickLength = [0.01 0.01];
    fig.Layer = 'top';
    fig.YLim = [-90 90];
    fig.YTick = [-89.99, -45, 0, 45, 89.99];
    fig.YTickLabel = {'90\circS','45\circS','0\circ','45\circN','90\circN'};
    fig.XTick = [0, 90, 180, 270, 360];
    fig.XTickLabel = {};
    fig.XGrid = 'on';
    fig.YGrid = 'on';
    fig.GridLineStyle = '--';
    fig.GridAlpha = 0.5;
    fig.Color = [0.9 0.9 0.9];
    cb = colorbar;
    cb.FontSize = 14;
    cb.Ticks = [1,2,3,4,5,6,7,8,9,10];
    cb.TickLabels = {'0.00','0.10','0.1125','0.125','0.1375','0.15','0.1625','0.175','0.20','0.30'};
    cb.Label.String = 'e-ratio';
    cb.Label.FontSize = 14;

    figure(3)
    Area = grid.Areat;
    Area(isnan(eratio_180x360.m)) = NaN;
    zonal_mean(:,i) = sum(eratio_180x360.m.*Area,2,'omitnan')./sum(Area,2,'omitnan');
    plot(eratio_180x360.y,zonal_mean(:,i),"LineWidth",3)
    hold on
    fig = gca;
    fig.FontSize = 12;
    fig.FontWeight = 'bold';
    fig.XLim = [-90 90];
    fig.YLim = [0.05 0.25];
    fig.XLabel.String = 'Latitude';
    fig.YLabel.String = 'Zonal average e-ratio';
    variance(i) = var(eratio_180x360.m(:),'omitnan');
end
legend('VarAll','FixAll','Location','northwest')
camroll(90)
set(gca,'XDir','reverse');

%% Plot e-ratio
% % Interpolate to 90 x 180 degree grid
cmap = cbrewer('div','RdBu',17);
cmap = [cmap(1:8,1:3);cmap(10:17,1:3)];
% cmap = cbrewer('div','RdBu',15);
% cmap = [cmap(1:7,1:3);cmap(9:15,1:3)];
cmap = flipud(cmap);

%eratio_change = ((e_ratio(:,:,1) - e_ratio(:,:,2))./e_ratio(:,:,2)).*100;
eratio_change = e_ratio(:,:,1) - e_ratio(:,:,2);
eratio_change(eratio_change<=-0.07) = -0.069;
%eratio_change(eratio_change>=50) = 49.9;
eratio_change_180x360 = interpTo90x180(lon,lat,eratio_change,'nearest');
load woa_grid.mat
eratio_change_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

% Create Figure
eratio = figure(4);
set(gcf,'Position',[500 100 1000 500],'Color','white')
%v = [-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08];
%v = [-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07];
v = [-0.08,-0.065,-0.05,-0.035,-0.02,-0.015,-0.01,-0.005,0,0.005,0.01,0.015,0.02,0.035,0.05,0.065,0.08];
%v = [-50,-37.5,-25,-18.75,-12.5,-9.375,-6.25,-3.125,0,3.125,6.25,9.375,12.5,18.75,25,37.5,50];
contourfnu(eratio_change_180x360.x,eratio_change_180x360.y,eratio_change_180x360.m,v,cmap,'eastoutside','true','pcolor')

% add coastline
hold on
[~,c] = contour(grid.XT,grid.YT,M3d(:,:,1),[1 1],'k');
c.LineWidth = 2;

% Format figure text and axis labels
fig = gca;
fig.FontSize = 12;
fig.FontWeight = 'bold';
fig.TickDir = 'out';
fig.TickLength = [0.01 0.01];
fig.Layer = 'top';
fig.YLim = [-90 90];
fig.YTick = [-89.99, -45, 0, 45, 89.99];
fig.YTickLabel = {'90\circS','45\circS','0\circ','45\circN','90\circN'};
fig.XTick = [0, 90, 180, 270, 360];
fig.XTickLabel = {};
fig.XGrid = 'on';
fig.YGrid = 'on';
fig.GridLineStyle = '--';
fig.GridAlpha = 0.5;
fig.Color = [0.9 0.9 0.9];

cb = colorbar;
cb.FontSize = 14;
%cb.Ticks = [2,4,6,8,10,12,14];
%v = [-0.08,-0.07,-0.05,-0.03,-0.02,-0.015,-0.01,-0.005,0,0.005,0.01,0.015,0.02,0.03,0.05,0.07,0.08];
cb.Ticks = [1,3,5,7,9,11,13,15,17];
cb.TickLabels = {'-0.08','-0.05','-0.02','-0.01','0.0','0.01','0.02','0.05','0.08'};
%cb.TickLabels = {'-0.08','-0.06','-0.04','-0.02','0.0','0.02','0.04','0.06','0.08'};
%cb.TickLabels = {'-50.0','-25.0','-12.5','-6.25','0.00','6.25','12.5','25.0','50.0'};
cb.Label.String = 'C) Change in e-ratio, varAll - fixAll';
cb.Label.FontSize = 14;

saveas(eratio,'figures\eratio_change_fixAll_varAll_abs.png')