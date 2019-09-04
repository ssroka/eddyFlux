function [d] = deg_2_m(pt1,pt2)
% pt1 = [lon lat] in degrees
% pt2 = [lon lat] in degrees
% [x_m,y_m] the horizontal and vertical distance between pt1 and pt2
% from https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
R = 6378.137; %km
d_lat = (pt2(2) - pt1(2));
d_lon = (pt2(1) - pt1(1));

a = sind(d_lat/2).^2 + cosd(pt1(2))*cosd(pt2(2))*sind(d_lon/2)*sind(d_lon/2);
c = 2*atan2(sqrt(a),sqrt(1-a));
d = R*c * 1000;

end