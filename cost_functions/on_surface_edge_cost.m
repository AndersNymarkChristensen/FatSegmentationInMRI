function cost = on_surface_edge_cost(terrain,mask,sigma)
%cost Generates on-surface cost from volume date by edge detection
%   COST = detect_edge(TERRAIN,MASK,SIGMA)
%   Usage, e.g.:
%   MASK = [1 -1] for detecting transition to darker
%   MASK = [-1 1] for detecting transition to lighter
%   MASK = [1 1 1 0 0 0] for detecting edge with 3 dark pixels below
%   MASK = [-1 -1 -1 0 0 0] for detecting edge 3 light pixels below
%   MASK = [0 0 0 1 1 1] for detecting edge 3 dark pixels above
%   MASK = [0 0 0 -1 -1 -1] for 3 light pixels above
%   MASK = 1; for detecting darkest
%   MASK = -1; for detecting lightest

if nargin<3 || isempty(sigma)
    sigma = 1;
end

if sigma>0
    n = max(6*sigma-1,5);
    z_mask = permute(fspecial('gaussian',[n 1],sigma),[3 2 1]);
    % a bit of smoothing
    terrain = imfilter(terrain,z_mask,'conv','replicate');
end

% now filtering with the mask, maximal response is at an edge
mask = reshape(mask,[1 1 numel(mask)]);
cost = imfilter(terrain,mask,'conv','replicate');
cost = (cost-min(cost(:)))/(max(cost(:))-min(cost(:)));
