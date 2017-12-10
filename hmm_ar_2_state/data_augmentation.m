function x =  data_augmentation(orgX, thetaE1, thetaE2, sigmaE1, sigmaE2, means, pr)

yt = orgX';
thetaE(:,:,1) = thetaE1;
thetaE(:,:,2) = thetaE2;
sigmaE(:,:,1) = sigmaE1;
sigmaE(:,:,2) = sigmaE2;

% clear all
% load('../0403_hmm_hmmar/hmmar2_0228_randStart_parameter')
% load('samples/hmmar2_nan_sample1.mat')
[n, m] = size(yt);

mp_state = ones(n-1,1);
mp_state(find(pr(1,:) < 0.5)) = 2;
mp_state = [0; mp_state];

m = m-1;
x1 = zeros(m,n);
xt = zeros(m,n);
xn = zeros(m,n); 
P1 = zeros(m,m,n);
Pt = zeros(m,m,n);
Pn = zeros(m,m,n);
K = zeros(m,m,n);
J = zeros(m,m,n);

% [~, indy1, ~] = intersect(sampleDfull(:,1), sampleD(:,1));
% [~, indy2, ~] = setxor(sampleDfull(:,1), sampleD(:,1));
% yt = NaN(n,m+1);
% yt(indy1,:) = sampleD;
% yt(:,1) = 0:0.5:(n-1)*0.5;

% A1 = diag([1 1 1 1]); A2 = diag([1 1 0 0]); A3 = diag([0 0 1 1]); A4 = diag([0 0 0 0]);
% B1 = diag([0 0 0 0]); B2 = diag([0 0 1 1]); B3 = diag([1 1 0 0]); B4 = diag([1 1 1 1]);

xt(:,1) = yt(1,2:end)';
Pt(:,:,1) = eye(m);
for i = 2 : n
    s = mp_state(i);
    temp = means(:,s) + thetaE(:,:,s)*xt(:,i-1);
    x1(:,i) = [temp];
    P1(:,:,i) = thetaE(:,:,s)*Pt(:,:,i-1)*thetaE(:,:,s)' + sigmaE(:,:,s);
    
%     if isnan(yt(i,1)) == true
%         A = A4; B = B4;
%     elseif isnan(yt(i,2)) == true
%         A = A3; B = B3;
%     elseif isnan(yt(i,4)) == true
%         A = A2; B = B2;
%     else
%         A = A1; B = B1;
%     end

    B = isnan(yt(i,2:end));
    A = (B == false);
    A = diag(A); B = diag(B);

    K(:,:,i) = P1(:,:,i)*A*inv(A*P1(:,:,i)*A + B);
    cy = yt(i,2:end);
    cy(find(isnan(cy))) = 0;
    xt(:,i) = x1(:,i) + K(:,:,i)*(cy'-A*x1(:,i));
    Pt(:,:,i) = (eye(m,m) - K(:,:,i)*A)*P1(:,:,i);
end

xn(:,n) = xt(:,n);
Pn(:,:,n) = Pt(:,:,n);
for i = n : -1 : 3
    s = mp_state(i);
    J(:,:,i-1) = Pt(:,:,i-1)*thetaE(:,:,s)*inv(P1(:,:,i));
    xn(:,i-1) = xt(:,i-1) + J(:,:,i-1)*(xn(:,i) - x1(:,i));
    Pn(:,:,i-1) = Pt(:,:,i-1) + J(:,:,i-1)*(Pn(:,:,i) - P1(:,:,i))*J(:,:,i-1)';    
end
xn(:,1) = xt(:,1);

x = [0:n-1; xn];

% figure
% hold on
% plot(0.5:0.5:(n-1)*0.5, xn(1,2:end),'o')
% plot(0.5:0.5:(n-1)*0.5, xn(1,2:end)+sqrt(abs(squeeze(Pn(1,1,2:end))')), ':')
% plot(0.5:0.5:(n-1)*0.5, xn(1,2:end)-sqrt(abs(squeeze(Pn(1,1,2:end))')), ':')

% error1 = (x1 - sampleDfull(:,2:end)');%./sampleDfull(:,2:end)';
% errort = (xt - sampleDfull(:,2:end)');%./sampleDfull(:,2:end)';
% errorn = (xn - sampleDfull(:,2:end)');%./sampleDfull(:,2:end)';
% 
% error1(:,1) = [];
% errort(:,1) = [];
% errorn(:,1) = [];
% 
% er1 = error1*error1'
% ert = errort*errort'
% ern = errorn*errorn'