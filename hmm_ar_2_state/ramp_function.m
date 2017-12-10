function data=ramp_function(N, y, sd)
% generate ramp data
% Input:
%   N - ramp length
%   y - ramp means
%   sd - standard deviation
%
%   Returns ramp data

n = floor(N/3);
data = zeros(1, n*3);

k = 1;
for i=1:3
    for j = 1:n
        data(1,k) = random('norm', y(i), sd, 1, 1);
        k = k + 1;
    end
end
