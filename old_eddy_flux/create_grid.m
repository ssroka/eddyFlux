function [X,Y] = create_grid(lat,lon)
% INPUT
%     lat: vector of latitudes (degrees)
%     lon: vector of longitude (degrees)
% OUTPUT
%      [X,Y] set of coordinates in cartesian with the origin at the left
%      and


X = zeros(length(lat),length(lon));
Y = zeros(length(lat),length(lon));

for i = 1:length(lat)
    pt1x = [lon(1) lat(i)];
    for j = 2:length(lon)
        pt2x = [lon(j) lat(i)];
        X(i,j) = deg_2_m(pt1x,pt2x);
    end
end

for i = 1:length(lon)
    pt1x = [lon(i) lat(1)];
    for j = 2:length(lat)
        pt2x = [lon(i) lat(j)];
        Y(j,i) = deg_2_m(pt1x,pt2x);
    end
end





end