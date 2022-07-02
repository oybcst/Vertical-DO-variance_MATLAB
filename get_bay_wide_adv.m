clear

matlist = dir('DOVar_*.mat');
lmat = length(matlist);
load(matlist(1).name,'lon_DOv','lat_DOv');
[M,N]=size(lon_DOv);
lon = lon_DOv;
lat = lat_DOv;

time = [];
BW_adv_bar = [];
totlen = 0;
for mati = 1:lmat
    load(matlist(mati).name,'ocean_time');
    load(matlist(mati).name,'adv_bar');
    lentime = length(ocean_time);
    time((totlen+1):(totlen+lentime),1) = ocean_time(:);
    BW_adv_bar(1:M,1:N,(totlen+1):(totlen+lentime)) = adv_bar(1:M,1:N,1:lentime);
    totlen = length(time);
end
totlen
BW_adv_bar_lp = NaN*BW_adv_bar;
for i = 1:M
    for j = 1:N
         x(1:totlen,1) = BW_adv_bar(i,j,1:totlen);
         if ~isnan(sum(x))
            [x_lp,~,~,~,~] = lanczosfilter(x,3600,1/48/3600,[],'low');
            BW_adv_bar_lp(i,j,1:totlen) = x_lp(1:totlen,1);
         end
     end
end
save('Bay_wide_201905_adv_bar.mat','time','lon','lat','BW_adv_bar','BW_adv_bar_lp');


