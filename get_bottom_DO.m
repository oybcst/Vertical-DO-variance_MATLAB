clear

ncpath = '/expanse/lustre/projects/ncs124/liuz1/ROMS_result/mobile_bio_Fennel/mobile_201905_org/';
gname = '/expanse/lustre/projects/ncs124/liuz1/ROMS_result/mobile_bio_Fennel/mobile_whole_grid11202021.nc';

nclist = dir([ncpath,'mobile_g01_his_*']);
lnc = length(nclist);

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

x1=unique(x1);
x2=unique(x2);
y1=unique(y1);
y2=unique(y2);

lon = lon_rho(x1:x2,y1:y2);
lat = lat_rho(x1:x2,y1:y2);

[M,N]=size(lon);

time = [];
DO_btm = [];
totlen = 0;

for nci = 1:lnc
    ocean_time =  ncread([ncpath,nclist(nci).name],'ocean_time');
    do_btm0 = ncread([ncpath,nclist(nci).name],'oxyg',[x1 y1 1 1],[M N 1 Inf]);
    lentime = length(ocean_time);
    time((totlen+1):(totlen+lentime),1) = ocean_time(:);
    DO_btm(1:M,1:N,(totlen+1):(totlen+lentime)) = do_btm0(1:M,1:N,1:lentime);
    totlen = length(time);
end

totlen

DO_btm_lp = NaN*DO_btm;

for i = 1:M
    for j = 1:N
         x(1:totlen,1) = DO_btm(i,j,1:totlen);
         if ~isnan(sum(x))
            [x_lp,~,~,~,~] = lanczosfilter(x,3600,1/48/3600,[],'low');
            DO_btm_lp(i,j,1:totlen) = x_lp(1:totlen,1);
         end
     end
end
save('Bay_wide_DO_btm_201905.mat','time','lon','lat','DO_btm','DO_btm_lp');
