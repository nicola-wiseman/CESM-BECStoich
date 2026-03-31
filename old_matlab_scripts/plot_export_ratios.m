%% Plot CESM-BGC Output for Wiseman et al., 2025
% This script generates the main figures and suplemental figures
% from Wiseman et al., 2025
close all; clearvars
addpath ../data/
addpath ../utils/
file_varAll = ('varAll.pop.h.0281_0300.nc');
file_inverse_NP = ('OCIM_N_model_results.nc');
file_inverse_CP = ('POC2POP_100m_weilei2023_gx3v7.nc');

%% CESM Model Data
lat = ncread(file_varAll,'TLAT');
lon = ncread(file_varAll,'TLONG');
depth = ncread(file_varAll,'z_t').*1e-2; % meters

% POM
pop_flux_in = ncread(file_varAll,'POP_FLUX_IN');
poc_flux_in = ncread(file_varAll,'POC_FLUX_IN');
pon_flux_in = ncread(file_varAll,'PON_FLUX_IN');
pfe_flux_in = ncread(file_varAll,'P_iron_FLUX_IN');

% flux at 100m
poc_flux_out_100m = poc_flux_in(:,:,depth==105);
pop_flux_out_100m = pop_flux_in(:,:,depth==105);
pon_flux_out_100m = pon_flux_in(:,:,depth==105);
pfe_flux_out_100m = pfe_flux_in(:,:,depth==105);

clear poc_flux_in pop_flux_in pon_flux_in pfe_flux_in depth

% stoichiometry at 100m
cn_flux_100m = poc_flux_out_100m./pon_flux_out_100m;
cp_flux_100m = poc_flux_out_100m./pop_flux_out_100m;
np_flux_100m = pon_flux_out_100m./pop_flux_out_100m;
fec_flux_100m = 1e6*pfe_flux_out_100m./poc_flux_out_100m;

clear poc_flux_out_100m pon_flux_out_100m pop_flux_out_100m pfe_flux_out_100m

%% Inverse Model Data
inverse_np = ncread(file_inverse_NP,'N2P_exp');
inverse_cp = ncread(file_inverse_CP,'POC2POP_100m');

%% Regrid Data to 180x360m
% Load WOA grid
load woa_grid.mat

% Interpolate to 180 x 360 degree grid
cn_180x360 = interpTo90x180(lon,lat,cn_flux_100m,'nearest');
cn_180x360.m(cn_180x360.m == 0) = NaN; % exclude inland seas
cn_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

cp_180x360 = interpTo90x180(lon,lat,cp_flux_100m,'nearest');
cp_180x360.m(cp_180x360.m == 0) = NaN; % exclude inland seas
cp_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

np_180x360 = interpTo90x180(lon,lat,np_flux_100m,'nearest');
np_180x360.m(np_180x360.m == 0) = NaN; % exclude inland seas
np_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land

fec_180x360 = interpTo90x180(lon,lat,fec_flux_100m,'nearest');
fec_180x360.m(fec_180x360.m == 0) = NaN; % exclude inland seas
fec_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
fec_180x360.m(fec_180x360.m > 300.0) = 300;

inverse_np_180x360 = interpTo90x180(lon,lat,inverse_np,'nearest');
inverse_np_180x360.m(inverse_np_180x360.m == 0) = NaN; % exclude inland seas
inverse_np_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
inverse_np_180x360.m(inverse_np_180x360.m > 50) = 50;

inverse_cp_180x360 = interpTo90x180(lon,lat,inverse_cp,'nearest');
inverse_cp_180x360.m(inverse_cp_180x360.m == 0) = NaN; % exclude inland seas
inverse_cp_180x360.m(M3d(:,:,1) == 0) = NaN; %exclude land
inverse_cp_180x360.m(inverse_cp_180x360.m < 90) = 90;

clear cn_flux_100m cp_flux_100m np_flux_100m fec_flux_100m ...
    inverse_cp inverse_np lat lon MSKS

%% Set up caxis
v_cn = [5.6,5.9,6.2,6.5,6.8,7.1,7.4,7.7,8.0,8.3,8.6,9.0];
v_cp = [90,100,110,120,130,140,160,180,200,225,250,300];
v_np = [8,10,12,14,16,18,20,22,24,26,30,50];
v_fec = [3.0,4.0,6.0,9.0,12.0,16.0,20.0,30.0,50.0,100,300];

c_cn = {'5.6','5.9','6.2','6.5','6.8','7.1','7.4','7.7','8.0','8.3','8.6','9.0'};
c_cp = {'90.','100','110','120','130','140','160','180','200','225','250','300'};
c_np = {'8.','10','12','14','16','18','20','22','24','26','30','50'};
c_fec = {'3.0','4.0','6.0','9.0','12.','16.','20.','30.','50.','100','300'};

%% Plot models

% plot CESM
f_cn = plot_ratio(cn_180x360,v_cn,'eastoutside',c_cn,'C/N',grid,M3d);
f_cp = plot_ratio(cp_180x360,v_cp,'eastoutside',c_cp,'C/P',grid,M3d);
f_np = plot_ratio(np_180x360,v_np,'eastoutside',c_np,'N/P',grid,M3d);
f_fec = plot_ratio(fec_180x360,v_fec,'eastoutside',c_fec,'Fe/C',grid,M3d);

% plot inverse model
f_inverse_cp = plot_ratio(inverse_cp_180x360,v_cp,'eastoutside',c_cp,'C/P',grid,M3d);
f_inverse_np = plot_ratio(inverse_np_180x360,v_np,'eastoutside',c_np,'N/P',grid,M3d);

% export files
saveas(f_cn,'figures\sinking_CN_varAll.png')
saveas(f_cp,'figures\sinking_CP_varAll.png')
saveas(f_np,'figures\sinking_NP_varAll.png')
saveas(f_fec,'figures\sinking_FeC_varAll.png')
saveas(f_inverse_cp,'figures\sinking_CP_inverse.png')
saveas(f_inverse_np,'figures\sinking_NP_inverse.png')
%% Plotting function
function f = plot_ratio(ratio,v,cbar_loc,cbar_labels,cbar_unit,grid,M3d)
    % create figure
    f = figure();
    set(gcf,'Position',[500 100 1000 500],'Color','white')
    contourfnu(ratio.x,ratio.y,ratio.m,v,parula,'none','true','pcolor')
    hold on
    % add coastline
    [~,c] = contour(grid.XT,grid.YT,M3d(:,:,1),[1 1],'k');
    c.LineWidth = 2;
    % format figure text and axis labels
    fig = gca;
    fig.FontSize = 18;
    fig.FontWeight = 'bold';
    fig.FontName = 'Arial';
    fig.TickDir = 'out';
    fig.TickLength = [0.01 0.01];
    fig.Layer = 'top';
    fig.YTick = [-89.5, -45, 0, 45, 89.5];
    fig.YTickLabel = {'90\circS','45\circS','0\circ','45\circN','90\circN'};
    fig.XTick = [0, 90, 180, 270, 360];
    fig.XTickLabel = {};
    fig.XGrid = 'on';
    fig.YGrid = 'on';
    fig.GridLineStyle = '--';
    fig.GridAlpha = 0.5;
    fig.Color = [0.9 0.9 0.9];

    cb = colorbar;
    cb.Location = cbar_loc;
    cb.FontSize = 22;
    cb.Ticks = 1:length(cbar_labels);
    cb.TickLabels = cbar_labels;
    cb.Label.String = cbar_unit;
    cb.Label.FontSize = 24;

end
