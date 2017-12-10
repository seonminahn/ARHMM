function error = comErrorVec1(a,x,y,diffT)

% load var4DXYT
% diffT = diffT*2;

[n, m] = size(x);

error = 0;
for i = 1 : m
    temp = y(:,i) - (a^diffT(i))*x(:,i);
    error = error + temp'*temp;
end