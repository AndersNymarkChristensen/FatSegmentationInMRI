function cost_R = in_region_membership_cost(V,div,delta)
% region cost based on a fuzzy membership function
% div is a ordered division between regions
% delta is a fuzziness parameter
% vand@dtu.dk

K = numel(div);
cost_R = zeros([size(V),K+1]);

% assumes input is a volume, so cost is 4D
cost_R(:,:,:,1) = 1 - fuzzy_membership_function(V,[min(V(:)),div(1)-delta(1)],delta([1,1]));
for r = 2:K
    cost_R(:,:,:,r) = 1 - fuzzy_membership_function(V,...
        [div(r-1)+delta(r-1),div(r)-delta(r)],delta(r-1:r));
end
cost_R(:,:,:,K+1) = 1 - fuzzy_membership_function(V,[div(K)+delta(K),max(V(:))],delta([K,K]));

