function [ierr] = DOVar(parm)
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


lon_DOv = lon_rho(x1:x2,y1:y2);
lat_DOv = lat_rho(x1:x2,y1:y2);
pm = ncread(gname,'pm',[x1 y1],[x2-x1+1 y2-y1+1]);
pn = ncread(gname,'pn',[x1 y1],[x2-x1+1 y2-y1+1]);
h = ncread(gname,'h',[x1 y1],[x2-x1+1 y2-y1+1]);

nci = parm;

fname = [ncpath,nclist(nci).name];

hc = double(ncread(fname,'hc'));
s_rho = double(ncread(fname,'s_rho'));
s_w = double(ncread(fname,'s_w'));
Cs_r = double(ncread(fname,'Cs_r'));
Cs_w = double(ncread(fname,'Cs_w'));

ocean_time = ncread(fname,'ocean_time');
zeta =double(ncread(fname,'zeta',[x1 y1 1],[x2-x1+1 y2-y1+1 Inf]));
DO = double(ncread(fname,'oxyg',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]));
urho = double(ncread(fname,'u_eastward',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]));
vrho = double(ncread(fname,'v_northward',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]));
Akt = double(ncread(fname,'AKt',[x1 y1 1 1],[x2-x1+1 y2-y1+1 Inf Inf]));
%Akt_bak = double(ncread(fname,'Akt_bak',[2],[1]));

wetdry_mask_rho = double(ncread(fname,'wetdry_mask_rho',[x1 y1 1],[x2-x1+1 y2-y1+1 Inf]));
wetdry_mask_rho(wetdry_mask_rho<1)=NaN;
[M,N,sNum,tNum] = size(DO);

%%
depth = zeros(M,N,tNum)*NaN;
DOMean  = zeros(M,N,tNum)*NaN;
DOVar  = zeros(M,N,tNum)*NaN;
DOVar_bar  = zeros(M,N,tNum)*NaN;
usp2  = zeros(M,N,tNum)*NaN;
vsp2  = zeros(M,N,tNum)*NaN;
upsp = zeros(M,N,tNum)*NaN;
vpsp = zeros(M,N,tNum)*NaN;
dissip = zeros(M,N,tNum)*NaN;


%%
for ti = 1:tNum

    for loni = 1:M
        for lati = 1:N

            if ~isnan(wetdry_mask_rho(loni, lati,ti))
                DOz(1:sNum) = DO(loni,lati,1:sNum,ti);
                z(1:sNum) = sigma2z(hc,s_rho,Cs_r,h(loni,lati),zeta(loni,lati,ti),sNum);
                z_w(1:sNum+1) = sigma2z(hc,s_w,Cs_w,h(loni,lati),zeta(loni,lati,ti),sNum+1);
                depth(loni,lati,ti)= abs(z(end)-z(1));
                
                DOMean(loni,lati,ti) = simps(z,DOz)/abs(z(end)-z(1));
                DOPrime = DOz-DOMean(loni,lati,ti);
                DOPrimeSq = DOPrime.^2;
                DOVar(loni,lati,ti) = simps(z,DOPrimeSq);
                DOVar_bar(loni,lati,ti) = DOVar(loni,lati,ti)/abs(z(end)-z(1));

                uz(1:sNum) = urho(loni,lati,1:sNum,ti);
                vz(1:sNum) = vrho(loni,lati,1:sNum,ti);
                uPrime = uz-simps(z,uz)/abs(z(end)-z(1));
                vPrime = vz-simps(z,vz)/abs(z(end)-z(1));

                usp2(loni,lati,ti) = simps(z,uz.*DOPrimeSq);
                vsp2(loni,lati,ti) = simps(z,vz.*DOPrimeSq);

                upsp(loni,lati,ti) = simps(z,-2*uPrime.*DOPrime);
                vpsp(loni,lati,ti) = simps(z,-2*vPrime.*DOPrime);

                Akt0(1:sNum+1) = Akt(loni,lati,1:sNum+1,ti);

                dDOdz_z_w(1:sNum+1) = 0;
                dDOdz_z_w(2:sNum)= (DOz(2:sNum)-DOz(1:sNum-1))./(z(2:sNum)-z(1:sNum-1));
               % dDOdz_z_w(1) = O2_btm(loni,lati,ti)/24/3600/Akt0(1);
               % dDOdz_z_w(sNum+1) = O2_airsea(loni,lati,ti)/24/3600/Akt0(end);

                dDOdz_z_w(1) = dDOdz_z_w(2);
                dDOdz_z_w(sNum+1) = dDOdz_z_w(sNum);
                dissip(loni,lati,ti) = simps(z_w,2*Akt0.*(dDOdz_z_w).^2);
            end
        end
    end
end


%% calcuate graident d /dx

dusp2dx  = zeros(M,N,tNum)*NaN;
dvsp2dy  = zeros(M,N,tNum)*NaN;
dDOMeandx  = zeros(M,N,tNum)*NaN;
dDOMeandy  = zeros(M,N,tNum)*NaN;

deltax1 = 1./pm(1,1:end);
deltax2 = 1./pm(2,1:end);
deltax3 = 1./pm(3,1:end);

dusp2dx(1,1:N,1:tNum) = ((usp2(2,1:N,1:tNum)-usp2(1,1:N,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4 -((usp2(3,1:N,1:tNum)-usp2(1,1:N,1:tNum)).*(deltax1+deltax2).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

dDOMeandx(1,1:N,1:tNum) = ((DOMean(2,1:N,1:tNum)-DOMean(1,1:N,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4 -((DOMean(3,1:N,1:tNum)-DOMean(1,1:N,1:tNum)).*(deltax1+deltax2).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

for i = 2 : M-1
    deltax1 = 1./pm(i-1,1:end);
    deltax2 = 1./pm(i,1:end);
    deltax3 = 1./pm(i+1,1:end);

    dusp2dx(i,1:N,1:tNum) = ((usp2(i+1,1:N,1:tNum)-usp2(i,1:N,1:tNum)).*(deltax1+deltax2).^2/4 -((usp2(i-1,1:N,1:tNum)-usp2(i,1:N,1:tNum)).*(deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

    dDOMeandx(i,1:N,1:tNum) =((DOMean(i+1,1:N,1:tNum)-DOMean(i,1:N,1:tNum)).*(deltax1+deltax2).^2/4 -((DOMean(i-1,1:N,1:tNum)-DOMean(i,1:N,1:tNum)).*(deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

end

deltax1 = 1./pm(end-2,1:end);
deltax2 = 1./pm(end-1,1:end);
deltax3 = 1./pm(end,1:end);

dusp2dx(end,1:N,1:tNum) = ((usp2(end-2,1:N,1:tNum)-usp2(end,1:N,1:tNum)).*(deltax2+deltax3).^2/4 -((usp2(end-1,1:N,1:tNum)-usp2(end,1:N,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

dDOMeandx(end,1:N,1:tNum) = ((DOMean(end-2,1:N,1:tNum)-DOMean(end,1:N,1:tNum)).*(deltax2+deltax3).^2/4 -((DOMean(end-1,1:N,1:tNum)-DOMean(end,1:N,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

%% calcuate graident d /dy

deltax1 = 1./pn(1:end,1);
deltax2 = 1./pn(1:end,2);
deltax3 = 1./pn(1:end,3);

dvsp2dy(1:M,1,1:tNum) = ((vsp2(1:M,2,1:tNum)-vsp2(1:M,1,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4 -((vsp2(1:M,3,1:tNum)-vsp2(1:M,1,1:tNum)).*(deltax1+deltax2).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

dDOMeandy(1:M,1,1:tNum) = ((DOMean(1:M,2,1:tNum)-DOMean(1:M,1,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4 -((DOMean(1:M,3,1:tNum)-DOMean(1:M,1,1:tNum)).*(deltax1+deltax2).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

for i = 2 : N-1
    deltax1 = 1./pn(1:end,i-1);
    deltax2 = 1./pn(1:end,i);
    deltax3 = 1./pn(1:end,i+1);

    dvsp2dy(1:M,i,1:tNum) = ((vsp2(1:M,i+1,1:tNum)-vsp2(1:M,i,1:tNum)).*(deltax1+deltax2).^2/4 -((vsp2(1:M,i-1,1:tNum)-vsp2(1:M,i,1:tNum)).*(deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

    dDOMeandy(1:M,i,1:tNum) = ((DOMean(1:M,i+1,1:tNum)-DOMean(1:M,i,1:tNum)).*(deltax1+deltax2).^2/4 -((DOMean(1:M,i-1,1:tNum)-DOMean(1:M,i,1:tNum)).*(deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);
end

deltax1 = 1./pn(1:end,end-2);
deltax2 = 1./pn(1:end,end-1);
deltax3 = 1./pn(1:end,end);

dvsp2dy(1:M,end,1:tNum) = ((vsp2(1:M,end-2,1:tNum)-vsp2(1:M,end,1:tNum)).*(deltax2+deltax3).^2/4 -((vsp2(1:M,end-1,1:tNum)-vsp2(1:M,end,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);
dDOMeandy(1:M,end,1:tNum) = ((DOMean(1:M,end-2,1:tNum)-DOMean(1:M,end,1:tNum)).*(deltax2+deltax3).^2/4 -((DOMean(1:M,end-1,1:tNum)-DOMean(1:M,end,1:tNum)).*(deltax1+2*deltax2+deltax3).^2/4))./((deltax1+deltax2).*(deltax1+2*deltax2+deltax3).*(deltax2+deltax3)/8);

DOVar_bar = DOVar./depth;
dissip_bar = dissip./depth;

%% adv term, straing and dissipation terms
adv = dusp2dx + dvsp2dy;
adv_bar = adv./depth;

strain = upsp.*dDOMeandx + vpsp.*dDOMeandy;
strain_bar = strain./depth;

matname = ['DOVar_',nclist(nci).name(end-4:end-3),'.mat'];

save(matname,'ocean_time','lon_DOv','lat_DOv','depth',...
        'DOVar','DOVar_bar',...
        'DOMean','urho','vrho',...
        'adv','adv_bar','usp2','vsp2','dusp2dx','dvsp2dy',...
        'strain','strain_bar','upsp','vpsp','dDOMeandx','dDOMeandy',...
        'dissip','dissip_bar');

end
