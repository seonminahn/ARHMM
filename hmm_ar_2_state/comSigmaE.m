% function sigmaE = comSigmaE(theta,x,y,diffT,pr)
% 
% [n, m] = size(x);
% 
% thetaTilda = cell(4,2);
% thetaTilda{1,1} = eye(n);
% thetaTilda{2,1} = eye(n) + theta;
% thetaTilda{3,1} = eye(n) + theta + theta^2;
% thetaTilda{4,1} = eye(n) + theta + theta^2 + theta^3;
% for i = 1 : 4
%    thetaTilda{i,2} = inv(thetaTilda{i,1});
% end
% 
% diffYX = zeros(n,m);
% tempSigma = zeros(n,n);
% for i = 1 : m
%     diffYX(:,i) = y(:,i) - (theta^diffT(i))*x(:,i);
%     temp = thetaTilda{diffT(i),2}*diffYX(:,i);
%     tempSigma = tempSigma + pr(i)*temp*temp';
% end
% sigmaE = tempSigma/sum(pr);


function sigmaE = comSigmaE(ctheta,x,y,diffT,pr)
[n, m] = size(x);

theta = ctheta(:,2:end);
M = ctheta(:,1);

thetaTilda = cell(4,2);
thetaTilda{1,1} = eye(n);
thetaTilda{2,1} = eye(n) + theta;
thetaTilda{3,1} = eye(n) + theta + theta^2;
thetaTilda{4,1} = eye(n) + theta + theta^2 + theta^3;
for i = 1 : 4
   thetaTilda{i,2} = inv(thetaTilda{i,1});
end

diffYX = zeros(n,m);
tempSigma = zeros(n,n);
for i = 1 : m
    diffYX(:,i) = y(:,i) - [thetaTilda{diffT(i),1}*M theta^diffT(i)]*[1; x(:,i)];
    temp = thetaTilda{diffT(i),2}*diffYX(:,i);
    tempSigma = tempSigma + pr(i)*(temp*temp');
end
sigmaE = tempSigma/sum(pr);

