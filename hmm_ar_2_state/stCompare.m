clear all
close all

load('4DdataWtime_norm.mat');   % fourvariablesCEL
data = fourvariablesCEL';
data = fliplr(data);
data(1,:) = 2*abs(data(1,:) - data(1,1));

load('0228_randStart/0228_randStart_parameter.mat')
nSample = 1;

st=backsample(data, nSample, sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means);

diffT = diff(data(1,:));

Ey1 = zeros(4, 601);
Ey2 = zeros(4, 601);

nvars = 4;
maxTdiff = 4;
tTheta1 = zeros(nvars, nvars, maxTdiff);
tTheta2 = zeros(nvars, nvars, maxTdiff);
tTheta1(:,:,1) = eye(nvars);
tTheta2(:,:,1) = eye(nvars);
for i = 2 : maxTdiff
    tTheta1(:,:,i) = tTheta1(:,:,i-1) + thetaE1^(i-1);
    tTheta2(:,:,i) = tTheta2(:,:,i-1) + thetaE2^(i-1);
end

for i = 1 : 600
    Ey1(:,i+1) = tTheta1(:,:,diffT(i))*means(:,1) + (thetaE1^diffT(i))*data(2:5,i);
    Ey2(:,i+1) = tTheta2(:,:,diffT(i))*means(:,2) + (thetaE2^diffT(i))*data(2:5,i);
end

st = fliplr(st);
data = fliplr(data);
Ey1 = fliplr(Ey1);
Ey2 = fliplr(Ey2);


%%
find1D = find(st(1,:)== 1);
find2D = find(st(1,:)== 2);


for i = 1 : 4
    figure(1)
    subplot(2,2,1)
    plot(find1D, Ey1(i,find1D), 'bo-', find1D, data(i+1,find1D), 'k*-')
    legend('Estimation', 'Exact value')
    title('State 1')
    subplot(2,2,3)
    plot(find1D, data(i+1,find1D)-Ey1(i,find1D), 'rs')
    legend('Exact value - Estimation')
    xlabel('depth')
    a = axis;
    title(['std(Exact value - Estimation) = ' num2str(std(data(i+1,find1D)-Ey1(i,find1D))) '      sqrt(sigmaE1(', num2str(i), ',' num2str(i) ')) = ', num2str(sqrt(sigmaE1(i,i)))])
    
    subplot(2,2,2)
    plot(find2D, Ey2(i,find2D), 'bo-', find2D, data(i+1,find2D), 'k*-')
    legend('Estimation', 'Exact value')
    title('State 2')
    subplot(2,2,4)
    plot(find2D, data(i+1,find2D)-Ey2(i,find2D), 'rs')    
    legend('Exact value - Estimation')
    xlabel('depth')
    axis(a)
    title(['std(Exact value - Estimation) = ' num2str(std(data(i+1,find2D)-Ey2(i,find2D))) '      sqrt(sigmaE2(', num2str(i), ',' num2str(i) ')) = ', num2str(sqrt(sigmaE2(i,i)))])
    
    export_fig(['stCompare' num2str(i) '.pdf'])
    
    pause
end


% x1 = [];
% x2 = [];
% y1 = [];
% y2 = [];
% yy1 = [];
% yy2 = [];
% dim  = 2;
% for i = 1 : 1
%     find1D = find(st(i,:)== 1);
%     x1 = [x1 find1D];
%     y1 = [y1 data(dim,find1D)];
%     yy1 = [yy1 Ey1(dim-1,find1D)];
%     
%     find2D = find(st(i,:)== 2);
%     x2 = [x2 find2D];
%     y2 = [y2 data(dim,find2D)];
%     yy2 = [yy2 Ey2(dim-1,find2D)];
% end
% 
% % Color encoding
% figure 
% hold on
% 
% A=[x1' y1'];
% [Auniq,~,IC] = unique(A,'rows');
% cnt = accumarray(IC,1);
% scatter(Auniq(:,1), Auniq(:,2), 50, cnt, 'fill');
% colormap(summer(max(cnt)));
% A=[x1' yy1'];
% [Auniq,~,IC] = unique(A,'rows');
% cnt = accumarray(IC,1);
% scatter(Auniq(:,1), Auniq(:,2), 50, cnt);
% colormap(summer(max(cnt)));
% 
% % Color encoding
% figure 
% hold on
% 
% A=[x2' y2'];
% [Auniq,~,IC] = unique(A,'rows');
% cnt = accumarray(IC,1);
% scatter(Auniq(:,1), Auniq(:,2), 50, cnt, 'fill');
% colormap(summer(max(cnt)));
% A=[x2' yy2'];
% [Auniq,~,IC] = unique(A,'rows');
% cnt = accumarray(IC,1);
% scatter(Auniq(:,1), Auniq(:,2), 50, cnt);
% colormap(summer(max(cnt)));


% % jitter
% figure
% scatter(x,y, 'jitter','on', 'jitterAmount', 1);
%
% % Number next to marker
% figure
% scatter(Auniq(:,1), Auniq(:,2));
% for ii=1:numel(cnt)
%     if cnt(ii)>1
%         text(Auniq(ii,1)+0.2,Auniq(ii,2),num2str(cnt(ii)), ...
%             'HorizontalAlignment','left', ...
%             'VerticalAlignment','middle', ...
%             'FontSize', 6);
%     end
% end
% 
% 
% % Number inside marker
% figure
% scatter(Auniq(:,1), Auniq(:,2), (6+2*(cnt>1)).^2); % make the ones where we'll put a number inside a bit bigger
% for ii=1:numel(cnt)
%     if cnt(ii)>1
%         text(Auniq(ii,1),Auniq(ii,2),num2str(cnt(ii)), ...
%             'HorizontalAlignment','center', ...
%             'VerticalAlignment','middle', ...
%             'FontSize', 6);
%     end
% end



