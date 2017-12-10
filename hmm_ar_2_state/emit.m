% % % function p=emit(x, m, cov)
% % % % emission probability
% % % % input:
% % % %   x - data value, nvars by 1
% % % %   m - mean
% % % %   cov - covariance matrix
% % % %
% % % % output: emission probability (pdf)
% % %     p=mvnpdf(x, m, cov);
% % %     if p < 1.0e-300  % this is a cheat to keep us out of numerical trouble
% % %         p = 1e-300;
% % %     end
% % % end

% function p = emit(x, theta, tSigma, M)
% % emission probability
% % input: 
% %   x - data values, (1+nvars) by L
% % output: emission probability
% 
% [nvars, L] = size(x);
% nvars = nvars-1;
% m = zeros(1,nvars);
% p = zeros(1,L);
% p(1,1) = 1;
% for i = 2 : L
%     timeDiff = x(1,i) - x(1,i-1);
%     data = x(2:5,i) - (theta^timeDiff)*x(2:5,i-1) - M;
%     p(1,i) = mvnpdf(data', m, tSigma(:,:,timeDiff)); 
%     if p(1,i) < 1.0e-300  % this is a cheat to keep us out of numerical trouble
%         p(1,i) = 1e-300;
%     end
% end

function p = emit(x, theta, tSigma, M)
% emission probability
% input: 
%   x - data values, (1+nvars) by L
% output: emission probability

[nvars, L] = size(x);
nvars = nvars-1;
m = zeros(1,nvars);
p = zeros(1,L);
p(1,1) = 1;

tTheta(:,:,1) = eye(nvars);
timeDiff = [0 diff(x(1,:))];
maxTdiff = max(timeDiff);
for i = 2 : maxTdiff
    tTheta(:,:,i) = tTheta(:,:,i-1) + theta^(i-1);
end

for i = 2 : L
    timeDiff = x(1,i) - x(1,i-1);
    cTheta = [tTheta(:,:,timeDiff)*M theta^timeDiff];
    data = x(2:end,i) - cTheta*[1; x(2:end,i-1)];
    p(1,i) = mvnpdf(data', m, tSigma(:,:,timeDiff)); 
    if p(1,i) < 1.0e-300  % this is a cheat to keep us out of numerical trouble
        p(1,i) = 1e-300;
    end
end