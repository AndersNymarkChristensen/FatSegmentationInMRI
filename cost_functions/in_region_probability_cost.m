function cost_R = in_region_probability_cost(V,mu,sigma)
% in region cost based on normalized probability
% vand@dtu.dk

R = numel(mu);
probability = zeros([size(V),R]);

% assumes input is a volume, so probability is 4D
for r = 1:R
    probability(:,:,:,r) = 1/(sigma(r)*(2*pi).^0.5)*exp(-(V-mu(r)).^2/(2*sigma(r).^2));    
end
% normalizing probability
probability = probability./repmat(sum(probability,4),[1 1 1 R]);
cost_R = 1-probability;

