function unfolded = unfold_volume_Jasmin(V,angle_resolution,radius_resolution,cent)
% unfolding x-y-height volume to cylindric coordinates 
% and permuting dimensions so we end up with angle-height-radius
V = double(V);
[X,Y,Z] = size(V);

angles = 2*pi*(0:1/angle_resolution:1-1/angle_resolution);
radius_max = ceil(norm([max(X-cent(2), cent(2)) max(Y-cent(1), cent(1))])); % maximal radius from axis to corner
radii =  0:radius_max/(radius_resolution-1):radius_max;

x_matrix = radii(:)*sin(angles)+cent(1); 
y_matrix = radii(:)*cos(angles)+cent(2);

% unfolded = zeros(radius_resolution,angle_resolution);
unfolded = zeros(radius_resolution,length(angles));
for z=1:Z
    unfolded(:,:,z) = interp2(V(:,:,z),x_matrix,y_matrix,'linear',0);
end
unfolded = permute(unfolded,[2,3,1]); % radius pointing upwards