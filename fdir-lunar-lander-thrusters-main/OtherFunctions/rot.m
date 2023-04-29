function [R]= rot(theta,type)

% anticlockwise +ve

if type==3
R = [cosd(theta) sind(theta) 0;-sind(theta) cosd(theta) 0;0 0 1];           %C

elseif type==2
R = [cosd(theta) 0 -sind(theta);0 1 0;sind(theta) 0 cosd(theta)];  

elseif type==1
R = [1 0 0;0 cosd(theta) sind(theta);0 -sind(theta) cosd(theta)];           %AC
end

end