function E = flatten_volume(V,s,q)
% for 10 voxels below surface use q = -9.5:-0.5
% for 10 voxels above surface use q = 0.5:9.5

[X,Y,Z] = size(V); 
Q = numel(q);
E = zeros(X,Y,Q);

for x=1:X
    for y=1:Y
        h = s(x,y)+q; % s is between regions, i.e. 
        h = min(max(1,h),Z); % not to fall out of V
        E(x,y,:) = V(x,y,h);
    end
end
