function unfolded = unfold_volume(V,angle_resolution,radius_resolution,cent)
% unfolding x-y-height volume to cylindric coordinates 
% and permuting dimensions so we end up with angle-height-radius
V = double(V);
[X,Y,Z] = size(V);

% % Use center of mass found by Otsu threshold
% cent = zeros(Z,2);
% for z = 1:Z
%     level = graythresh(V(:,:,z))*max(max(V(:,:,z)));
% [y, x] = find( V(:,:,z) > level );
% cent(z,:) = [mean(x) mean(y)];
% end
% cent = mean(cent,1);

% TEST
angles = [linspace(0,angle_resolution*1.5/12,angle_resolution*2/12),...
    linspace((angle_resolution)*1.5/12,(angle_resolution)*4.5/12,angle_resolution*2/12+1),... 
    linspace((angle_resolution)*4.5/12,(angle_resolution)*7.5/12,angle_resolution*4/12+1),... 
    linspace((angle_resolution)*7.5/12,(angle_resolution)*10.5/12,angle_resolution*2/12+1),... 
   linspace((angle_resolution)*10.5/12,(angle_resolution),angle_resolution*2/12+1)];
angles = 2*pi*(angles./angle_resolution);
angles = unique(angles);

% original
% angles = 2*pi*(0:1/angle_resolution:1-1/angle_resolution);
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


%% Plot

% figure
% hold on
% imagesc(V(:,:,5)');colormap('gray');axis off;axis image; set(gca,'position',[0 0 1 1],'units','normalized');
% plot(y_matrix(:,1:5:end),x_matrix(:,1:5:end),'y')
% hold off

% figure;imagesc(squeeze(unfolded(:,10,:)));colormap('gray');axis off;axis image; set(gca,'position',[0 0 1 1],'units','normalized');


% figure;plot(x_matrix(:,5),y_matrix(:,5),'y')


%% Original
% angles = 2*pi*(0:1/angle_resolution:1-1/angle_resolution);
% radius_max = ceil(0.5*norm([X Y])); % maximal radius from axis to corner
% radii =  0:radius_max/(radius_resolution-1):radius_max;
% x_matrix = radii(:)*sin(angles)+(Y+1)/2; 
% y_matrix = radii(:)*cos(angles)+(X+1)/2;
% unfolded = zeros(radius_resolution,angle_resolution);
% for z=1:Z
%     unfolded(:,:,z) = interp2(V(:,:,z),x_matrix,y_matrix,'linear',0);
% end
% unfolded = permute(unfolded,[2,3,1]); % radius pointing upwards

% NOTE: interp2 assumes that x_matrix and y_matrix are in meshgrid format 
%   (as opposed to ngrid format) so the data in unfolded is flipped. Still,
%   as long as it is handled correctly when folding back, all is ok.
%   However, consider using griddedInterpoland instead.