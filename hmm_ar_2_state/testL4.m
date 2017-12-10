function p = testL4(x, thetaE1, sigmaE1, thetaE2, sigmaE2, a, pr, M)

nStates = size(a,1);
timeDiff = [0 diff(x(1,:))];
nvars = size(x,1)-1;
L=size(x,2);
maxTdiff = max(timeDiff);

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

em = zeros(nStates, L);
em(1,:) = emit(x, thetaE1, tSigma1, M(:,1));
em(2,:) = emit(x, thetaE2, tSigma2, M(:,2));

% [f, p1, scale]=n_forward(em, a, timeDiff);
% b=n_backward(em, a, scale, timeDiff);
% pr = f.*b;
% pr(:,1) = [];
p = sum(sum(pr.*log(em(:,2:end))));



