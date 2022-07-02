clear

matlist = dir('Bio_terms_*.mat');
lmat = length(matlist);
load(matlist(1).name,'column_NEM_bar');
[M,N,~]=size(column_NEM_bar);
%lon = lon_sv;
%lat = lat_sv;

time = [];
BW_column_NEM_bar = [];
totlen = 0;
for mati = 1:lmat
%    load(matlist(mati).name,'ocean_time');
    load(matlist(mati).name,'column_NEM_bar');
    lentime = 24;
   % time((totlen+1):(totlen+lentime),1) = ocean_time(:);
    BW_column_NEM_bar(1:M,1:N,(totlen+1):(totlen+lentime)) = column_NEM_bar(1:M,1:N,1:lentime);
    totlen = totlen+24;
end
totlen
BW_column_NEM_bar_lp = NaN*BW_column_NEM_bar;
for i = 1:M
    for j = 1:N
         x(1:totlen,1) = BW_column_NEM_bar(i,j,1:totlen);
         if ~isnan(sum(x))
            [x_lp,~,~,~,~] = lanczosfilter(x,3600,1/48/3600,[],'low');
            BW_column_NEM_bar_lp(i,j,1:totlen) = x_lp(1:totlen,1);
         end
     end
end
save('Bay_wide_201905_column_NEM_bar.mat','BW_column_NEM_bar','BW_column_NEM_bar_lp');


