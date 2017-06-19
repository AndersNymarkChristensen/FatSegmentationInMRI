function [V,R,cost_s,cost_r] = peaks_test_volume(mu,sigma,dim)
% returns test volume with three terren-like layers (two surfaces)

if nargin<3
    dim = [20,24,30];
end
% surface 1 is peaks
[x,y] = ndgrid(linspace(-2,2,dim(1)),linspace(-2,2,dim(2)));
s1 = peaks(x,y);
s1_normalized = (s1-min(s1(:)))/(max(s1(:))-min(s1(:)));
S1 = round(dim(3)*(2*s1_normalized + 1)/5); % from 1/5 to 3/5 volume height 
% surface 2 is peaks rotated for 90 deg
[y,x] = ndgrid(linspace(-2,2,dim(2)),linspace(-2,2,dim(1)));
s2 = peaks(y,x)';
s2_normalized = (s2-min(s2(:)))/(max(s2(:))-min(s2(:)));
S2 = round(dim(3)*(2*s2_normalized + 2)/5); % from 2/5 to 4/5 volume height
% volume of heighs 
h_vol = repmat(permute(1:dim(3),[1 3 2]),[dim(1),dim(2),1]);
S1_vol = repmat(S1,[1 1 dim(3)]);
S2_vol = repmat(S2,[1 1 dim(3)]);

% on-surface ground truth cost (i.e. zeros where the surface is)
cost_s(:,:,:,1) = ~(h_vol==S1_vol);
cost_s(:,:,:,2) = ~(h_vol==S2_vol);

% in-region ground truth cost (i.e. zeros where the region is)
cost_r(:,:,:,1) = ~(h_vol<=S1_vol);
cost_r(:,:,:,2) = ~(S1_vol<h_vol & h_vol<=S2_vol);
cost_r(:,:,:,3) = ~(S2_vol<h_vol);

% test volume with regions having distributions N(mu,sigma)
% and ground-truth volume with regions 1 2 3
V = zeros(dim);
R = zeros(dim);
for i=1:3   
    n = sum(~reshape(cost_r(:,:,:,i),[prod(dim),1])); 
    V(~cost_r(:,:,:,i)) = sigma(i)*randn(n,1)+mu(i);
    R(~cost_r(:,:,:,i)) = i;
end

