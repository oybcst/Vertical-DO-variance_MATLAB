function [ierr] = O2_prod(parm)
%% compute colum NEM
parm

ncpath = '/expanse/lustre/projects/ncs124/liuz1/ROMS_result/mobile_bio_Fennel/mobile_201905_org/';
gname = '/expanse/lustre/projects/ncs124/liuz1/ROMS_result/mobile_bio_Fennel/mobile_whole_grid11202021.nc';

nclist = dir([ncpath,'mobile_dia_*']);


lon_rho = ncread(gname,'lon_rho');
lat_rho = ncread(gname,'lat_rho');

xv = [-88.3 -88.3 -87.7 -87.7 -88.3];
yv = [30.1 30.75 30.75 30.1 30.1];

in = inpolygon(lon_rho,lat_rho,xv,yv);

lon_rho_bay = lon_rho(in);
lat_rho_bay = lat_rho(in);

[x1,~]=find((lon_rho.*double(in))==min(lon_rho_bay));
[x2,~]=find((lon_rho.*double(in))==max(lon_rho_bay));
[~,y1]=find((lat_rho.*double(in))==min(lat_rho_bay));
[~,y2]=find((lat_rho.*double(in))==max(lat_rho_bay));

x1 = unique(x1);
x2 = unique(x2);
y1 = unique(y1);
y2 = unique(y2);

nci = parm;
fname = [ncpath,nclist(nci).name];

    ocean_time = ncread(fname,'ocean_time');
    NO3_uptake = ncread(fname,'NO3_uptake',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]);
    P_Production = ncread(fname,'P_Production',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]);

    rOxNO3 = 8.625;
    rOxNH4 = 6.625;

    NH4_uptake = P_Production - NO3_uptake;
    O2_so = NO3_uptake*rOxNO3 + NH4_uptake*rOxNH4;

    matname = ['O2_so_',nclist(nci).name(end-4:end-3),'.mat'];
    save(matname,'ocean_time','O2_so');
end
