function cost = on_surface_cumulative_cost(terrain,change_direction,change_sign)
%surface_cumulative Generates on-surface cost from volume date using cumulative sum
%   COST = surface_cumulative(TERRAIN,CHANGE_DIRECTION,CHANGE_SIGN)
%   Usage
%   CHANGE_DIRECTION = true for summing from above, else summing from below
%   CHANGE_SIGN = true for summing over negative values, penalizing dark voxels 

if change_direction
    terrain = flipdim(terrain,3);
end

if change_sign
    terrain = max(terrain(:))-terrain;
end

cost = cumsum(terrain,3);
if change_direction
    cost = flipdim(cost,3);
end

cost = (cost-min(cost(:)))/(max(cost(:))-min(cost(:)));
