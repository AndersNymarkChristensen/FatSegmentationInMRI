function XYZ = fold_back(s,dimV,radius_resolution,angle_resolution,cent)


radius_max = ceil(norm([max(dimV(1)-cent(2), cent(2)) max(dimV(2)-cent(1), cent(1))])); % max radius used when unfolding
k = size(s,3); % number of detected surfaces
% angle_resolution = size(s,1); % number of angles used when unfolding
z_resolution = size(s,2); % number of MR slices
% angles = 2*pi*(0:1/angle_resolution:1-1/angle_resolution);
angles = [linspace(0,angle_resolution*1.5/12,angle_resolution*2/12),...
    linspace((angle_resolution)*1.5/12,(angle_resolution)*4.5/12,angle_resolution*2/12+1),... 
    linspace((angle_resolution)*4.5/12,(angle_resolution)*7.5/12,angle_resolution*4/12+1),... 
    linspace((angle_resolution)*7.5/12,(angle_resolution)*10.5/12,angle_resolution*2/12+1),... 
   linspace((angle_resolution)*10.5/12,(angle_resolution),angle_resolution*2/12+1)];
angles = 2*pi*(angles./angle_resolution);
angles = unique(angles);

angles_mat = repmat(angles(:),[1,z_resolution,k]);

r = (s-1)*radius_max/(radius_resolution-1); % from voxel indices to radius
XYZ(:,:,1,:) = r.*sin(angles_mat)+cent(1); 
XYZ(:,:,2,:) = r.*cos(angles_mat)+cent(2);
XYZ(:,:,3,:) = repmat(1:size(s,2),[size(s,1),1,k]);


%% original
% radius_max = ceil(0.5*norm(dimV(1:2))); % max radius used when unfolding
% k = size(s,3); % number of detected surfaces
% angle_resolution = size(s,1); % number of angles used when unfolding
% z_resolution = size(s,2); % number of MR slices
% angles = 2*pi*(0:1/angle_resolution:1-1/angle_resolution);
% angles_mat = repmat(angles(:),[1,z_resolution,k]);
% 
% r = (s-1)*radius_max/(radius_resolution-1); % from voxel indices to radius
% XYZ(:,:,1,:) = r.*sin(angles_mat)+(dimV(2)+1)/2; 
% XYZ(:,:,2,:) = r.*cos(angles_mat)+(dimV(1)+1)/2;
% XYZ(:,:,3,:) = repmat(1:size(s,2),[size(s,1),1,k]);

