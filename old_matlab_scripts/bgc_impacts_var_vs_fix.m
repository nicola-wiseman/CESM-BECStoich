%% Biogeochemical Impacts of Variable vs Fixed Stoichiometry
clearvars; close all

addpath ../data/
addpath ../utils/
file1 = ('/home/mv23682/Documents/CNP_Manuscript/model_output/varAll.pop.h.0281_0300.nc');
file2 = ('/home/mv23682/Documents/CNP_Manuscript/model_output/fixAll.pop.h.0281_0300.nc');
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
Nfix = NaN([100 116 length(files)]);
WCdenit = NaN([100 116 length(files)]);
POC_100m = NaN([100 116 length(files)]);
Fe_scavenge_100m = NaN([100 116 length(files)]);
Fe_uptake_100m = NaN([100 116 length(files)]);
PON_100m = NaN([100 116 length(files)]);
POP_100m = NaN([100 116 length(files)]);
P_iron_100m = NaN([100 116 length(files)]);
SiO2_100m = NaN([100 116 length(files)]);
vol_low_O2 = NaN([length(files),1]);
for i = 1:length(files)
    % Unpack data
    % Carbon
    photoC_sp = ncread(files{i},'photoC_sp');
    photoC_diat = ncread(files{i},'photoC_diat');
    photoC_diaz = ncread(files{i},'photoC_diaz');
    photoFe_sp = ncread(files{i},'photoFe_sp');
    photoFe_diat = ncread(files{i},'photoFe_diat');
    photoFe_diaz = ncread(files{i},'photoFe_diaz');
    NPP(:,:,:,i) = photoC_sp + photoC_diat + photoC_diaz; % units mmol/m3/s
    NPP_100m(:,:,i) = sum(NPP(:,:,1:10,i),3); % sum over top 100m
    NPP_100m(:,:,i) = NPP_100m(:,:,i).*(12.01*1e-3*10.0*365.25*24*60*60); % convert to gC/m2/yr
    Fe_uptake = photoFe_sp + photoFe_diat + photoFe_diaz;
    Fe_uptake_100m(:,:,i) = sum(Fe_uptake(:,:,1:10),3);
    Fe_scavenge = ncread(files{i},'Fe_scavenge');
    Fe_scavenge_100m(:,:,i) = sum(Fe_scavenge(:,:,1:10),3);
    POC_FLUX_IN = ncread(files{i},'POC_FLUX_IN');
    POC_100m(:,:,i) = POC_FLUX_IN(:,:,depth==105);
    PON_FLUX_IN = ncread(files{i},'PON_FLUX_IN');
    PON_100m(:,:,i) = PON_FLUX_IN(:,:,depth==105);
    POP_FLUX_IN = ncread(files{i},'POP_FLUX_IN');
    POP_100m(:,:,i) = POP_FLUX_IN(:,:,depth==105);
    P_iron_FLUX_IN = ncread(files{i},'P_iron_FLUX_IN');
    P_iron_100m(:,:,i) = P_iron_FLUX_IN(:,:,depth==105);
    SiO2_FLUX_IN = ncread(files{i},'SiO2_FLUX_IN');
    SiO2_100m(:,:,i) = SiO2_FLUX_IN(:,:,depth==105);
    %Nitrogen
    diaz_Nfix = ncread(files{i},'diaz_Nfix');
    diaz_Nfix(diaz_Nfix<=0) = NaN;
    Nfix(:,:,i) = sum(diaz_Nfix(:,:,:),3,'omitnan'); % sum over all depths
    DENITRIF = ncread(files{i},'DENITRIF');
    WCdenit(:,:,i) = sum(DENITRIF(:,:,:),3,'omitnan'); % sum over all depths
    O2 = ncread(files{i},'O2');
    
    vol_low_O2(i) = sum(vol(O2<30),'all');


    clear photoC_sp photoC_diat photoC_diaz diaz_Nfix DENITRIF
end

NPP_change = ((NPP_100m(:,:,1) - NPP_100m(:,:,2))./NPP_100m(:,:,2)).*100;
POC_export_change = ((POC_100m(:,:,1) - POC_100m(:,:,2))./POC_100m(:,:,2)).*100;
Nfix_change = ((Nfix(:,:,1) - Nfix(:,:,2))./Nfix(:,:,2)).*100;
WCdenit_change = ((WCdenit(:,:,1) - WCdenit(:,:,2))./WCdenit(:,:,2)).*100;
WCdenit_change(WCdenit(:,:,1)==0) = 0;

POP_export_change = ((POP_100m(:,:,1) - POP_100m(:,:,2))./POP_100m(:,:,2)).*100;
PON_export_change = ((PON_100m(:,:,1) - PON_100m(:,:,2))./PON_100m(:,:,2)).*100;
Si_export_change = ((SiO2_100m(:,:,1) - SiO2_100m(:,:,2))./SiO2_100m(:,:,2)).*100;
P_iron_export_change = ((P_iron_100m(:,:,1) - P_iron_100m(:,:,2))./P_iron_100m(:,:,2)).*100;
frac_bio = Fe_uptake_100m./(Fe_uptake_100m+Fe_scavenge_100m);
P_iron_100m_bio = P_iron_100m.*frac_bio;
P_iron_bio_export_change = ((P_iron_100m_bio(:,:,1) - P_iron_100m_bio(:,:,2))./P_iron_100m_bio(:,:,2)).*100;

%% Plot NPP
% % Interpolate to 90 x 180 degree grid
cmap = cbrewer('div','RdBu',21);
cmap = [cmap(1:10,1:3);cmap(12:21,1:3)];
cmap = flipud(cmap);

NPP_change(NPP_change>=1000) = 999.9;
npp_180x360 = interpTo90x180(lon,lat,NPP_change,'nearest');
load woa_grid.mat
npp_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

% Create Figure
NPP = figure(1);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
%v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,25.0,37.5,50,75.0,100,150,200,350,500];
contourfnu(npp_180x360.x,npp_180x360.y,npp_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
%cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','25.0','50.0','100.','200.','500.'};
cb.Label.String = 'A) Change in NPP, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot Carbon Export
% % Interpolate to 90 x 180 degree grid

POC_export_change(POC_export_change>=1000) = 999.9;
export_180x360 = interpTo90x180(lon,lat,POC_export_change,'nearest');
export_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

% Create Figure
C_export = figure(2);
set(gcf,'Position',[500 100 1000 500],'Color','white')
%v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,62.5,75,87.5,100];
contourfnu(export_180x360.x,export_180x360.y,export_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
%cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','75.0','100.'};
cb.Label.String = 'B) Change in POC Export, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot Nfix
% % Interpolate to 90 x 180 degree grid

Nfix_change(Nfix_change>=300) = 299.9;
nfix_180x360 = interpTo90x180(lon,lat,Nfix_change,'nearest');
nfix_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

% Create Figure
Nfixation = figure(3);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,200,300];
%v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,62.5,75,87.5,100];
contourfnu(nfix_180x360.x,nfix_180x360.y,nfix_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','300.'};
cb.Label.String = 'Change in N Fixation, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot WC Denit
% % Interpolate to 90 x 180 degree grid

WCdenit_change(WCdenit_change == 0) = NaN;
wcdenit_180x360 = interpTo90x180(lon,lat,WCdenit_change,'nearest');
wcdenit_180x360.m(wcdenit_180x360.m==0) = NaN;
wcdenit_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
wcdenit_180x360.m(wcdenit_180x360.m>=1000) = 999.9;

% Create Figure
WCdenitrif = figure(4);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
contourfnu(wcdenit_180x360.x,wcdenit_180x360.y,wcdenit_180x360.m,v,cmap,'eastoutside','true','pcolor');

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
%cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','75.0','100.'};
cb.Label.String = 'Change in WC Denit, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot Nitrogen Export
% % Interpolate to 90 x 180 degree grid

n_export_180x360 = interpTo90x180(lon,lat,PON_export_change,'nearest');
n_export_180x360.m(n_export_180x360.m == 0) = NaN; % exclude inland seas
n_export_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
n_export_180x360.m(n_export_180x360.m >= 1000.0) = 999.9; %exclude land

% Create Figure
N_export = figure(5);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
contourfnu(n_export_180x360.x,n_export_180x360.y,n_export_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
cb.Label.String = 'Change in PON Export, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot Phosphorus Export
% % Interpolate to 90 x 180 degree grid

p_export_180x360 = interpTo90x180(lon,lat,POP_export_change,'nearest');
p_export_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
p_export_180x360.m(p_export_180x360.m >= 1000.0) = 999.9; %exclude land

% Create Figure
P_export = figure(6);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
contourfnu(p_export_180x360.x,p_export_180x360.y,p_export_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
cb.Label.String = 'Change in POP Export, varAll - fixAll (%)';
cb.Label.FontSize = 14;


%% Plot Iron Export
% % Interpolate to 90 x 180 degree grid

fe_export_180x360 = interpTo90x180(lon,lat,P_iron_export_change,'nearest');
fe_export_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
fe_export_180x360.m(fe_export_180x360.m >= 1000.0) = 999.9;

% Create Figure
Fe_export = figure(7);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
contourfnu(fe_export_180x360.x,fe_export_180x360.y,fe_export_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
cb.Label.String = 'Change in Iron Export, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot Bio Iron Export
% % Interpolate to 90 x 180 degree grid

fe_bio_export_180x360 = interpTo90x180(lon,lat,P_iron_bio_export_change,'nearest');
fe_bio_export_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
fe_bio_export_180x360.m(fe_bio_export_180x360.m >= 1000.0) = 999.9;

% Create Figure
Fe_bio_export = figure(8);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
contourfnu(fe_bio_export_180x360.x,fe_bio_export_180x360.y,fe_bio_export_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
cb.Label.String = 'Change in Bio Iron Export, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%% Plot Si Export
% % Interpolate to 90 x 180 degree grid

si_export_180x360 = interpTo90x180(lon,lat,Si_export_change,'nearest');
si_export_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
si_export_180x360.m(si_export_180x360.m >= 1000.0) = 999.9;
%si_export_180x360.m(si_export_180x360.m <= -100.0) = -100.0;

% Create Figure
Si_export = figure(9);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,75.0,100,550,1000];
contourfnu(si_export_180x360.x,si_export_180x360.y,si_export_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.Ticks = [1,3,5,7,9,11,13,15,17,19,21];
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','100.','1000.'};
cb.Label.String = 'Change in SiO_2 Export, varAll - fixAll (%)';
cb.Label.FontSize = 14;

%%
saveas(NPP,'figures\NPP_change_fixAll_varAll.png')
saveas(C_export,'figures\C_export_change_fixAll_varAll.png')
saveas(Nfixation,'figures\Nfix_change_fixAll_varAll.png')
saveas(WCdenitrif,'figures\WCdenit_change_fixAll_varAll.png')
saveas(N_export,'figures\N_export_change_fixAll_varAll.png')
saveas(P_export,'figures\P_export_change_fixAll_varAll.png')
saveas(Fe_export,'figures\Fe_export_change_fixAll_varAll.png')
saveas(Fe_bio_export,'figures\Fe_bio_export_change_fixAll_varAll.png')
saveas(Si_export,'figures\Si_export_change_fixAll_varAll.png')