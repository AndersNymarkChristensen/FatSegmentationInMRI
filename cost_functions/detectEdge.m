function edge = detectEdge(volume, direction)
% FInd edges/gradients in image
% INPUT:
%       volume: volume data in the orientation used på Graph-Cut
%       direction: defines when an edge is positive or negative
% OUTPUT:
%       Edge: Same size as volume
% Copyright Anders Nymark Christensen, anym@dtu.dk, DTU Compute 20170607

[~,Z,~] = size(volume);
edge = zeros(size(volume));
for z = 1:Z
    [FX, FY] = gradient(squeeze(volume(:,z,:)));
    edge(:,z,:) = direction*sign(FX).*(abs(FX)+abs(FY));
end