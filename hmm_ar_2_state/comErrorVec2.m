% function mLL = comErrorVec2(theta,x,y,diffT,sigmaE,pr)
% 
% [n, m] = size(x);
% 
% thetaTilda = cell(4,3);
% thetaTilda{1,1} = eye(n);
% thetaTilda{2,1} = eye(n) + theta;
% thetaTilda{3,1} = eye(n) + theta + theta^2;
% thetaTilda{4,1} = eye(n) + theta + theta^2 + theta^3;
% 
% for i = 1 : 4
%    thetaTilda{i,2} = inv(thetaTilda{i,1}*sigmaE*thetaTilda{i,1}');
%    thetaTilda{i,3} = log(det(thetaTilda{i,1}*thetaTilda{i,1}'));
% end
% 
% mLL1 = zeros(m,1);
% mLL2 = zeros(m,1);
% mLLsum = zeros(m,1);
% for i = 1 : m
%     mLL1(i) = thetaTilda{diffT(i),3};
%     
%     temp1 = y(:,i) - (theta^diffT(i))*x(:,i);
%     temp2 = temp1'*thetaTilda{diffT(i),2}*temp1;
%     mLL2(i) = temp2;
%     mLLsum(i) = pr(i)*(mLL1(i)+mLL2(i));
% end
% mLL = sum(mLLsum);

function mLL = comErrorVec2(ctheta,x,y,diffT,sigmaE,pr)

[n, m] = size(x);

theta = ctheta(:,2:end);
M = ctheta(:,1);

thetaTilda = cell(4,3);
thetaTilda{1,1} = eye(n);
thetaTilda{2,1} = eye(n) + theta;
thetaTilda{3,1} = eye(n) + theta + theta^2;
thetaTilda{4,1} = eye(n) + theta + theta^2 + theta^3;

for i = 1 : 1
   thetaTilda{i,2} = inv(thetaTilda{i,1}*sigmaE*thetaTilda{i,1}');
   thetaTilda{i,3} = log(det(thetaTilda{i,1}*thetaTilda{i,1}'));
end

mLL1 = zeros(m,1);
mLL2 = zeros(m,1);
mLLsum = zeros(m,1);
for i = 1 : m
    mLL1(i) = thetaTilda{diffT(i),3};
    
    temp1 = y(:,i) - [thetaTilda{diffT(i),1}*M theta^diffT(i)]*[1; x(:,i)];
    temp2 = temp1'*thetaTilda{diffT(i),2}*temp1;
    mLL2(i) = temp2;
    mLLsum(i) = pr(i)*(mLL1(i)+mLL2(i));
end
mLL = sum(mLLsum);