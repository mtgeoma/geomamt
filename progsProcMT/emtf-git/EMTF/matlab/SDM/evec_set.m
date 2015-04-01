%  sets up for eigenvector plotting ... call once to open file/intialize

[Hp,Ep,Hz] = ch_pair_def(nch,chid);

% some calculations for scaling plot sizes ...
lat_max =max(max(stcor(1,:)));
lat_min =min(min(stcor(1,:)));
lon_max =max(max(stcor(2,:)));
lon_min =min(min(stcor(2,:)));
lat_av = (lat_max+lat_min)/2.;
lat_range = (lat_max-lat_min);
lon_range = cos(lat_av*pi/180)*(lon_max-lon_min);
marg = .35*max([lat_range,lon_range]);
if(marg == 0) 
%  user forgot to define site locations in sp file ...  yet again !!!
%   define dummy site locations ...
   nn = size(stcor);nsta = nn(2);
   stcor = [ 1:nsta; 1:nsta ];
   lat_max =max(max(stcor(1,:)));
   lat_min =min(min(stcor(1,:)));
   lon_max =max(max(stcor(2,:)));
   lon_min =min(min(stcor(2,:)));
   lat_av = (lat_max+lat_min)/2.;
   lat_range = (lat_max-lat_min);
   lon_range = cos(lat_av*pi/180)*(lon_max-lon_min);
   marg = .35*max([lat_range,lon_range]);
end

asp = (2*marg+lon_range)/(2*marg+lat_range);
lat_max = lat_max+marg;
lon_max = lon_max+marg;
lat_min = lat_min-marg;
lon_min = lon_min-marg;
ll_lim = [ lat_min,lat_max,lon_min,lon_max];
