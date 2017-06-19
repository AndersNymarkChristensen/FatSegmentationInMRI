function out = fuzzy_membership_function(in,interval,delta)
% out is 1 in interval, and smoothly goes to 0 outside the interval
% delta is the width of decrease to 0.5

sigma = delta/(2*log(2))^0.5;
out = ones(size(in));
out(in<interval(1)) = exp(-(in(in<interval(1))-interval(1)).^2/(2*sigma(2)^2));
out(in>interval(2)) = exp(-(in(in>interval(2))-interval(2)).^2/(2*sigma(2)^2));
