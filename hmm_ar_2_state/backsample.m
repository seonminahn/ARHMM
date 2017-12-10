% function states = backsample(seq, num_samples, a, means, cov1, cov2)
% sample from forward matrix to get state probailities
% input:
%   seq - the data
%   num_samples - how many samplez to take
%   a - transition matrix
%   means - data means - nstates X nvars
%   cov1, cov2 - covariance matrix fo states 1 and 2
%
% returns;
%   states - a num_samples x length of data matrix
%           each row contains the state sampled at each position
%           dividing by num_samples will give the state probability

function states = backsample(seq, num_samples, sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means)
% clear all
% load('4DdataWtime_norm.mat');   % fourvariablesCEL
% seq = fourvariablesCEL';
% seq = fliplr(seq);
% seq(1,:) = 2*abs(seq(1,:) - seq(1,1));
% num_samples = 1000;
% load('iter27.mat')

timeDiff = [0 diff(seq(1,:))];
maxTdiff = max(timeDiff);

nvars = size(seq,1)-1;
L = size(seq, 2);
num_states = size(a,1);

% we could have pass this stuff in instead of calculating
% % % em = zeros(num_states, L+1);
% % % em(:,1) = 1;
% % % for count = 2:L+1
% % %      em(1, count) = emit(seq(:, count-1)', means(1, :), cov1);
% % %      em(2, count) = emit(seq(:, count-1)', means(2, :), cov2);
% % % end

tTheta1 = zeros(nvars, nvars, maxTdiff);
tTheta2 = zeros(nvars, nvars, maxTdiff);
tSigma1 = zeros(nvars, nvars, maxTdiff);
tSigma2 = zeros(nvars, nvars, maxTdiff);

tTheta1(:,:,1) = eye(nvars);
tTheta2(:,:,1) = eye(nvars);
tSigma1(:,:,1) = sigmaE1;
tSigma2(:,:,1) = sigmaE2;

for i = 2 : maxTdiff
    tTheta1(:,:,i) = tTheta1(:,:,i-1) + thetaE1^(i-1);
    tTheta2(:,:,i) = tTheta2(:,:,i-1) + thetaE2^(i-1);
    tSigma1(:,:,i) = tTheta1(:,:,i)*sigmaE1*tTheta1(:,:,i)';
    tSigma2(:,:,i) = tTheta2(:,:,i)*sigmaE2*tTheta2(:,:,i)';
end
em = zeros(num_states, L);
em(1,:) = emit(seq, thetaE1, tSigma1, means(:,1));
em(2,:) = emit(seq, thetaE2, tSigma2, means(:,2));

f = n_forward(sa, em, a, timeDiff);

% sample - this is SLOW
states = zeros(num_samples, L);
for n = 1:num_samples
    states(n,L) = randsample(1:num_states, 1, true, f(:,L));
    
    for i = L-1:-1:2
        aCurrent = a^timeDiff(i+1);  
        states(n,i) = randsample(1:num_states, 1, true, f(:,i).*aCurrent(:,states(n, i+1)));
    end
end

% plot(sum(states==1)/num_samples)