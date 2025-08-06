%% Biodiversity Impacts of Variable vs Fixed Stoichiometry
clearvars; close all

addpath ../data/
addpath ../utils/
file1 = ('/home/mv23682/Documents/CNP_Manuscript/model_output/varAll.pop.h.0281_0300.nc');
file2 = ('/home/mv23682/Documents/CNP_Manuscript/model_output/fixAll.pop.h.0281_0300.nc');
%file3 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixN.pop.h.0281_0300.nc');
%file4 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixP.pop.h.0281_0300.nc');
%file5 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixFe.pop.h.0281_0300.nc');
%file6 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixSi.pop.h.0281_0300.nc');
%file7 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/varN.pop.h.0281_0300.nc');
%file8 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/varP.pop.h.0281_0300.nc');
%file9 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/varFe.pop.h.0281_0300.nc');
%file10 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/varSi.pop.h.0281_0300.nc');
%file11 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixAll_red.pop.h.0281_0300.nc');
%file12 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixAll_mean.pop.h.0281_0300.nc');
%file13 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixAll_CMIP.pop.h.0281_0300.nc');
%file14 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/varAll_dynCO2.pop.h.0281_0300.nc');
%file15 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/fixAll_dynCO2.pop.h.0281_0300.nc');
%file16 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/aVarAll.pop.h.0281_0300.nc');
%file17 = ('C:/Users/niciw/OneDrive/Documents/MATLAB/data/aFixAll.pop.h.0281_0300.nc');


%files = {file1};
%files = {file1,file2,file3,file4,file5,file6,file7,file8,file9,file10,file11,file12,file13,file14,file15,file16,file17};
files = {file1,file2};
% General Variables
lat = ncread(file1,'TLAT');
lon = ncread(file1,'TLONG');
tarea = ncread(file1,'TAREA').*0.0001; % convert to meters
depth = ncread(file1,'z_t').*1e-2; % convert to meters
depth_150m = ncread(file1,'z_t_150m').*1e-2; % convert to meters

spC = NaN([100 116 15 length(files)]);
diatC = NaN([100 116 15 length(files)]);
diazC = NaN([100 116 15 length(files)]);
diaz_Nfix = NaN([100 116 15 length(files)]);
spC_zint = NaN([100 116 length(files)]);
diatC_zint = NaN([100 116 length(files)]);
diazC_zint = NaN([100 116 length(files)]);
NPP_zint = NaN([100 116 length(files)]);

spCtot = NaN([1 length(files)]);
diatCtot = NaN([1 length(files)]);
diazCtot = NaN([1 length(files)]);
POCfluxtot = NaN([1 length(files)]);
POCfluxtotPg = NaN([1 length(files)]);
PONfluxtot = NaN([1 length(files)]);
POPfluxtot = NaN([1 length(files)]);
Sifluxtot = NaN([1 length(files)]);
Fefluxtot = NaN([1 length(files)]);
NPPtot = NaN([1 length(files)]);
Nfixtot = NaN([1 length(files)]);
WCdenittot = NaN([1 length(files)]);
for n = 1:length(files)
    % Unpack data
    % Carbon
    sp_C = ncread(files{n},'spC'); % units mmol/m^3
    sp_C(sp_C<0) = 0;
    spC(:,:,:,n) = sp_C;
    clear sp_C
    diat_C = ncread(files{n},'diatC'); % units mmol/m^3
    diat_C(diat_C<0) = 0;
    diatC(:,:,:,n) = diat_C;
    clear diat_C;
    diaz_C = ncread(files{n},'diazC'); % units mmol/m^3
    diaz_C(diaz_C<=0) = NaN;
    diazC(:,:,:,n) = diaz_C;
    clear diaz_C
    POC_FLUX_IN = ncread(files{n},'POC_FLUX_IN'); % units nmol/cm^2/s
    PON_FLUX_IN = ncread(files{n},'PON_FLUX_IN'); % units nmol/cm2/s
    POP_FLUX_IN = ncread(files{n},'POP_FLUX_IN'); % units nmol/cm^2/s
    SiO2_FLUX_IN = ncread(files{n},'SiO2_FLUX_IN'); % units
    P_iron_FLUX_IN = ncread(files{n},'P_iron_FLUX_IN'); 
    photoC_sp = ncread(files{n},'photoC_sp');
    photoC_diat = ncread(files{n},'photoC_diat');
    photoC_diaz = ncread(files{n},'photoC_diaz');
    NPP = photoC_sp + photoC_diat + photoC_diaz; % units mmol/m3/s
    for i = 1:100
        for j = 1:116
            spC_zint(i,j,n) = sum(squeeze(spC(i,j,:,n)).*depth_150m,'omitnan');
            diatC_zint(i,j,n) = sum(squeeze(diatC(i,j,:,n)).*depth_150m,'omitnan');
            diazC_zint(i,j,n) = sum(squeeze(diazC(i,j,:,n)).*depth_150m,'omitnan');
            NPP_zint(i,j,n) = sum(squeeze(NPP(i,j,1:10)),'omitnan');
        end
    end
    spCtot(n) = sum(spC_zint(:,:,n).*tarea,'all','omitnan');
    diatCtot(n) = sum(diatC_zint(:,:,n).*tarea,'all','omitnan'); 
    diazCtot(n) = sum(diazC_zint(:,:,n).*tarea,'all','omitnan');
    NPPtot(n) = sum(NPP_zint(:,:,n).*tarea,'all','omitnan').*(12.01*1e-3*10.0*365.25*24*60*60*1e-15);
    POCfluxtotPg(n) = sum(POC_FLUX_IN(:,:,depth == 105).*tarea,'all','omitnan').*(12.01*1e-3*10.0*365.25*24*60*60*1e-18); % global sum, conver to TmolC/yr
    POCfluxtot(n) = sum(POC_FLUX_IN(:,:,depth == 105).*tarea,'all','omitnan').*(1e4*1e-21*365.25*24*60*60); % global sum, conver to TmolC/yr
    PONfluxtot(n) = sum(PON_FLUX_IN(:,:,depth == 105).*tarea,'all','omitnan').*(1e4*1e-21*365.25*24*60*60); % global sum, conver to TmolN/yr;
    POPfluxtot(n) = sum(POP_FLUX_IN(:,:,depth == 105).*tarea,'all','omitnan').*(1e4*1e-21*365.25*24*60*60); % global sum, conver to TmolP/yr;
    Sifluxtot(n) = sum(SiO2_FLUX_IN(:,:,depth == 105).*tarea,'all','omitnan').*(1e4*1e-21*365.25*24*60*60); % global sum, conver to TmolFe/yr
    Fefluxtot(n) = sum(P_iron_FLUX_IN(:,:,depth == 105).*tarea,'all','omitnan').*(1e4*1e-21*365.25*24*60*60*1e3); % global sum, conver to TmolSi/yr;
    diaz_Nfix(:,:,:,n) = ncread(files{n},'diaz_Nfix'); % mmol/m^3/s, 15 layers
    Nfix = sum(diaz_Nfix(:,:,:,n),3,'omitnan'); % sum over all depths
    Nfixtot(n) = sum(Nfix.*tarea,'all','omitnan').*(14.0067*1e4*1e-18*365.25*24*60*60); % global sum, TgN/yr;
    DENITRIF = ncread(files{n},'DENITRIF');
    WCdenit = sum(DENITRIF(:,:,:),3,'omitnan'); % sum over all depths
    WCdenittot(n) = sum(WCdenit.*tarea,'all','omitnan').*(14.0067*1e4*1e-18*365.25*24*60*60);
end
bioCtot = spCtot + diatCtot + diazCtot;
fracsp = spCtot./bioCtot;
fracdiat = diatCtot./bioCtot;
fracdiaz = diazCtot./bioCtot;


sp_bio_change = ((spC_zint(:,:,1) - spC_zint(:,:,2))./spC_zint(:,:,2)).*100;
%sp_bio_change = ((spC_zint(:,:,6) - spC_zint(:,:,1))./spC_zint(:,:,1)).*100;
sp_bio_change(sp_bio_change>=1000) = 999.9;
sp_bio_change(sp_bio_change<=-100) = NaN;
diat_bio_change = ((diatC_zint(:,:,1) - diatC_zint(:,:,2))./diatC_zint(:,:,2)).*100;
%diat_bio_change = ((diatC_zint(:,:,6) - diatC_zint(:,:,1))./diatC_zint(:,:,1)).*100;
diat_bio_change(diat_bio_change>=1000) = 999.9;
diat_bio_change(diat_bio_change<=-100) = NaN;
diaz_bio_change = ((diazC_zint(:,:,1) - diazC_zint(:,:,2))./diazC_zint(:,:,2)).*100;
%diaz_bio_change = ((diazC_zint(:,:,6) - diazC_zint(:,:,1))./diazC_zint(:,:,1)).*100;
diaz_bio_change(diaz_bio_change>=100) = 99.9;
diaz_bio_change(diaz_bio_change<=-100) = NaN;
%% Plot sp_bio_change
% % Interpolate to 90 x 180 degree grid

sp_bio_180x360 = interpTo90x180(lon,lat,sp_bio_change,'nearest');
sp_bio_180x360.m(sp_bio_180x360.m == 0) = NaN; % exclude inland seas
load woa_grid.mat
sp_bio_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
sp_bio_180x360.m(sp_bio_180x360.m>=1000) = 999.9;

%% Create Figure
sp_bio = figure(1);
set(gcf,'Position',[500 100 1000 500],'Color','white')
cmap = cbrewer('div','RdBu',21);
cmap = [cmap(1:10,1:3);cmap(12:21,1:3)];
cmap = flipud(cmap);
%v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,25.0,37.5,50,75.0,100,150,200,350,500];
v = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,25.0,37.5,50,75.0,100,300,500,750,1000];


contourfnu(sp_bio_180x360.x,sp_bio_180x360.y,sp_bio_180x360.m,v,cmap,'eastoutside','true','pcolor')

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
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','25.0','50.0','100.','500.','1000.'};
cb.Label.String = 'A) Change in Sp Biomass, varAll - fixAll (%)';
%cb.Label.String = 'Change in Sp Biomass, fixSi - varAll (%)';
cb.Label.FontSize = 14;

%% Plot diat_bio_change
% % Interpolate to 90 x 180 degree grid
diat_bio_180x360 = interpTo90x180(lon,lat,diat_bio_change,'nearest');
diat_bio_180x360.m(diat_bio_180x360.m == 0) = NaN; % exclude inland seas
diat_bio_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
diat_bio_180x360.m(diat_bio_180x360.m>=1000) = 999.9;

%% Create Figure
diat_bio = figure(2);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v2 = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,25.0,37.5,50,75.0,100,300,500,750,1000];
contourfnu(diat_bio_180x360.x,diat_bio_180x360.y,diat_bio_180x360.m,v2,cmap,'eastoutside','true','pcolor')

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
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','25.0','50.0','100.','500.','1000.'};
cb.Label.String = 'B) Change in Diat Biomass, varAll - fixAll (%)';
%cb.Label.String = 'Change in Diat Biomass, fixSi - varAll (%)';
cb.Label.FontSize = 14;

%% Plot diaz_bio_change
% % Interpolate to 90 x 180 degree grid
diaz_bio_180x360 = interpTo90x180(lon,lat,diaz_bio_change,'nearest');
diaz_bio_180x360.m(diaz_bio_180x360.m == 0) = NaN; % exclude inland seas
diaz_bio_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
diaz_bio_180x360.m(diaz_bio_180x360.m>=100) = 99.9;

%% Create Figure
diaz_bio = figure(3);
set(gcf,'Position',[500 100 1000 500],'Color','white')
v3 = [-100,-87.5,-75,-62.5,-50,-37.5,-25,-18.75,-12.5,-6.25,0,6.25,12.5,18.75,25,37.5,50,62.5,75,87.5,100];
contourfnu(diaz_bio_180x360.x,diaz_bio_180x360.y,diaz_bio_180x360.m,v3,cmap,'eastoutside','true','pcolor')

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
cb.TickLabels = {'-100.','-75.0','-50.0','-25.0','-12.5','0.00','12.5','25.0','50.0','75.0','100.'};
cb.Label.String = 'C) Change in Diaz Biomass, varAll - fixAll (%)';
%cb.Label.String = 'Change in Diaz Biomass, fixSi - varAll (%)';
cb.Label.FontSize = 14;
set(gcf,'Position',[500 100 1000 500],'Color','white')

saveas(sp_bio,'figures\sp_biomass_change_fixAll_varAll.png')
saveas(diat_bio,'figures\diat_biomass_change_fixAll_varAll.png')
saveas(diaz_bio,'figures\diaz_biomass_change_fixAll_varAll.png')