function unfolded = unfold_volumeANC(V, S, angle_resolution, radius_resolution)
% unfolding x-y-height volume to cylindric coordinates 
% and permuting dimensions so we end up with angle-height-radius
V = double(V);
[X,Y,Z] = size(V);



unfolded = zeros(radius_resolution,angle_resolution,Z);
for z=1:Z
    x = S(:,z,1);
    y = S(:,z,2);
    
    pathXY = [x, y];
    stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
    stepLengths = [0; stepLengths]; % add the starting point
    cumulativeLen = cumsum(stepLengths);
    finalStepLocs = linspace(0,cumulativeLen(end), angle_resolution);
    finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);
    center = mean(finalPathXY);
    center = [center(1)*ones(length(finalPathXY),1), center(2)*ones(length(finalPathXY),1)];
    
    outer = (finalPathXY - center).* 1.1 + center;
    inner = (finalPathXY - center).* 0.1 + center;
    if z == 10
        figure
        hold on
        imagesc(V(:,:,z))
        plot(outer(:,1),outer(:,2))
        plot(inner(:,1),inner(:,2))
    end
    
    
    for angle = 1:length(outer)
        unfolded(angle,z,:) = improfile(V(:,:,z),...
            [inner(angle,1),outer(angle,1)], [inner(angle,2),outer(angle,2)], radius_resolution, 'bilinear');
    end
end

figure
hold on
imagesc(V(:,:,end))
plot(outer(:,1),outer(:,2))
plot(inner(:,1),inner(:,2))

% unfolded = permute(unfolded,[2,3,1]); % radius pointing upwards

% NOTE: interp2 assumes that x_matrix and y_matrix are in meshgrid format 
%   (as opposed to ngrid format) so the data in unfolded is flipped. Still,
%   as long as it is handled correctly when folding back, all is ok.
%   However, consider using griddedInterpoland instead.



%%
% radii =  0:radius_max/(radius_resolution-1):radius_max;
% 
% x_angle = (ef.a * cos(angles) * cos(-ef.phi) - ef.b * sin(angles) * sin(-ef.phi))./max((ef.a * cos(angles) * cos(-ef.phi) - ef.b * sin(angles) * sin(-ef.phi)));
% y_angle = (ef.b * sin(angles) * cos(-ef.phi) + ef.a * cos(angles) * sin(-ef.phi))./max((ef.b * sin(angles) * cos(-ef.phi) + ef.a * cos(angles) * sin(-ef.phi)));
% 
% x_matrix = radii(:)*x_angle+(Y+1)/2; 
% y_matrix = radii(:)*y_angle+(X+1)/2;