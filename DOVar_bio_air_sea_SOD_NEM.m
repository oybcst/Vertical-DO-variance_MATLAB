function [ierr] = DOVar_bio_air_sea_SOD_NEM(parm)
parm

ncpath = '/expanse/lustre/projects/ncs124/liuz1/ROMS_result/mobile_bio_Fennel/mobile_201905_org/';
gname = '/expanse/lustre/projects/ncs124/liuz1/ROMS_result/mobile_bio_Fennel/mobile_whole_grid11202021.nc';
nclist = dir([ncpath,'mobile_g01_his_*']);

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


nci = parm;
fname = [ncpath,nclist(nci).name];

hc = double(ncread(fname,'hc'));
s_rho = double(ncread(fname,'s_rho'));
Cs_r = double(ncread(fname,'Cs_r'));

h = ncread(gname,'h',[x1 y1],[x2-x1+1 y2-y1+1]);
zeta = double(ncread(fname,'zeta',[x1 y1 1],[x2-x1+1 y2-y1+1 Inf]));
DO = double(ncread(fname,'oxyg',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]));

[M,N,tNum] = size(zeta);
sNum = 16;

DO_var_name = ['DOVar_',nclist(nci).name(end-4:end-3),'.mat'];
O2_so_name = ['O2_so_',nclist(nci).name(end-4:end-3),'.mat'];
O2_btm_name = ['O2_btm',nclist(nci).name(end-4:end-3),'.mat'];
O2_airsea_name = ['O2_airsea_',nclist(nci).name(end-4:end-3),'.mat'];

load(DO_var_name,'DOMean');
%time1 = ocean_time;

load(O2_btm_name,'O2_btm');
load(O2_airsea_name,'O2_airsea');
%time2 = ocean_time;
%lt2 = length(time2);

load(O2_so_name,'O2_so');
%time3 = ocean_time;
%lt3 = length(time3);

%%
air_sea_exchange = zeros(M,N,tNum)*NaN;
air_sea_exchange_bar = zeros(M,N,tNum)*NaN;

SOD = zeros(M,N,tNum)*NaN;
SOD_bar = zeros(M,N,tNum)*NaN;

column_NEM = zeros(M,N,tNum)*NaN;
column_NEM_bar = zeros(M,N,tNum)*NaN;

prod = O2_so;
prodMean = zeros(M,N,tNum)*NaN;

surf = O2_airsea;
btm = O2_btm;

%%
for ti = 1:tNum
    for loni = 1:M
        for lati = 1:N
            if ~isnan(DOMean(loni,lati,ti))

                z(1:sNum) = sigma2z(hc,s_rho,Cs_r,h(loni,lati),zeta(loni,lati,ti),sNum);
                DOz(1:sNum) = DO(loni,lati,1:sNum,ti);
                DOPrime = DOz-DOMean(loni,lati,ti);

                prodz(1:sNum) = prod(loni,lati,1:sNum,ti);
                prodMean(loni,lati,ti) = simps(z,prodz)/abs(z(end)-z(1));
                prodPrime = prodz-prodMean(loni,lati,ti);

                column_NEM(loni,lati,ti) = simps(z,2*DOPrime.*prodPrime);
                column_NEM_bar(loni,lati,ti) = column_NEM(loni,lati,ti)/abs(z(end)-z(1));

                air_sea_exchange(loni,lati,ti) = 2*DOPrime(end)*surf(loni,lati,ti);
                air_sea_exchange_bar(loni,lati,ti) = air_sea_exchange(loni,lati,ti)/abs(z(end)-z(1));

                SOD(loni,lati,ti) = -2*DOPrime(1)*btm(loni,lati,ti);
                SOD_bar(loni,lati,ti) = SOD(loni,lati,ti)/abs(z(end)-z(1));
            end
        end
    end
end

%ocean_time = time1;
matname = ['Bio_terms_',nclist(nci).name(end-4:end-3),'.mat'];
save(matname,'air_sea_exchange','air_sea_exchange_bar',...
    'SOD','SOD_bar','column_NEM','column_NEM_bar');

end

