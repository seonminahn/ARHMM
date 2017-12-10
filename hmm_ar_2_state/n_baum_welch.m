% % % function [loglik, a, means, cov1, cov2]=n_baum_welch(x, m1, m2, cv1, cv2)

function [loglik, sa, a, thetaE1, thetaE2, sigmaE1, sigmaE2, means, x]=n_baum_welch(fnum, x, thetaE1, thetaE2, sigmaE1, sigmaE2)
% clear all
% load('4DdataWtime_norm.mat');   % fourvariablesCEL
% x = fourvariablesCEL';
% x = fliplr(x);
% x(1,:) = 2*abs(x(1,:) - x(1,1));
% load('parameterE.mat'); % sigmaE1, sigmaE2, thetaE1, thetaE2
% fnum = 0;

% debugON = false;
% if debugON
%     displayOption = 'iter';
% else
%     displayOption = 'final';
% end
displayOption = 'final';

% meansUpdate = true;

timeDiff = [0 diff(x(1,:))];
maxTdiff = max(timeDiff);
% thetaE1 = eye(4);
% thetaE2 = eye(4);
% preThetaE1 = thetaE1;
% preThetaE2 = thetaE2;
% preSigmaE1 = sigmaE1;
% preSigmaE2 = sigmaE2;

% Baum-Welch algorithm for HMM with multivariate-normal emissions
%
% input:
% x - data matrix N x L, n observations of length L
% m1, m2 - initial means for state 1 and state 2
% cv1, cv2 - initial covariance matrices for states 1 and 2
%
% return:
% loglik - log likelihood of data
% a - estimated transition probailities
% means - estimated emissions for each state
% cov1, cov2 - estimated covariannce for ststes 1 and 2
%
% requires emit.m, n_forward.m, n_backwards,m

iter_max = 1000;
nStates = 2;
% nvars = size(x, 1);
nvars = size(x,1)-1;

L=size(x,2);

pseudo = 0.01;
a_pseudo = zeros(nStates, nStates);
a_pseudo(1,1) = pseudo;
a_pseudo(1,2) = pseudo;
a_pseudo(2,2) = pseudo;
a_pseudo(2,1) = pseudo;

% % generate some initial guess for the transition matix
a=zeros(nStates,nStates);
arandSize = 0.2;
a(1,1) = 1 - arandSize * random('unif', 0, 1, 1, 1);
a(1,2) = random('unif', 0, arandSize, 1, 1);
a(2,1) = random('unif', 0, arandSize, 1, 1);
a(2,2) = 1 - arandSize * random('unif', 0, 1, 1, 1);
for i=1:nStates
    a(i,:)=a(i,:)/sum(a(i,:));
end
sa(1,1) = rand;
sa(2,1) = 1-sa(1,1);
% if debugON
% %     a = [0.9047 0.0953; 0.0526 0.9474];
%     a = [0.8570 0.1430; 0.1511 0.8489];
%     sa = [0.5825; 0.4175];
% end

dataD1 = zeros(nvars,L-1);
dataD2 = zeros(nvars,L-1);
for l = 2 : L
    dataD1(:,l) = x(2:end,l) - (thetaE1^timeDiff(l))*x(2:end,l-1);
    dataD2(:,l) = x(2:end,l) - (thetaE2^timeDiff(l))*x(2:end,l-1);
end
% if meansUpdate == true
%     means(:,1) = nanmean(dataD1,2);
%     means(:,2) = nanmean(dataD2,2);
% else
%     means = zeros(nvars,2);
% end
means(:,1) = nanmean(dataD1,2);
means(:,2) = nanmean(dataD2,2);

% % % means(1, :) = m1;
% % % means(2, :) = m2;
% % % cov1 = cv1;
% % % cov2 = cv2;


% make sure the covariances are legit
[~, err] = cholcov(sigmaE1, 0);
if err ~= 0
    disp('cov1 is not positive definite')
    loglik = -realmax;
    return;
end

[~, err] = cholcov(sigmaE2, 0);
if err ~= 0
    disp('cov2 is not positive definite')
    loglik = -realmax;
    return;
end

loglik= -realmax;
comdiff=100000;
iter = 0;
% this is the main EM loop
% we use teh difference between log likleihoods as a stopping criterie
% we could use a ratio, but this works ok
orgX = x;
pr = zeros(nStates, L-1);

pr(1,:) = 1;
x01 = data_augmentation(orgX, thetaE1, thetaE2, sigmaE1, sigmaE2, means, pr);
pr(1,:) = 0;
pr(2,:) = 1;
x02 = data_augmentation(orgX, thetaE1, thetaE2, sigmaE1, sigmaE2, means, pr);
x = (x01+x02)/2;
timeDiff = [0 diff(x(1,:))];
maxTdiff = max(timeDiff);

while (abs(comdiff)>0.0001) && (iter < iter_max)
    toc  
    
    % initialize
% % %     A=zeros(nStates,nStates);
    A=a_pseudo;
% % %     M = zeros(nStates, nvars);
% % %     COV1 = zeros(nvars, nvars);
% % %     COV2 = zeros(nvars, nvars);
    
    % get the emission probabilities give the current parameters
% % %     em = zeros(nStates, L+1);
% % %     em(:,1) = 1;
% % %     for count = 2:L+1
% % %         em(1, count) = emit(x(:, count-1)', means(1, :), cov1);
% % %         em(2, count) = emit(x(:, count-1)', means(2, :), cov2);
% % %     end
    
    
%     tTheta1 = zeros(nvars, nvars, maxTdiff);
%     tTheta2 = zeros(nvars, nvars, maxTdiff);
%     tSigma1 = zeros(nvars, nvars, maxTdiff);
%     tSigma2 = zeros(nvars, nvars, maxTdiff);
%     
%     tTheta1(:,:,1) = eye(nvars);
%     tTheta2(:,:,1) = eye(nvars);
%     tSigma1(:,:,1) = sigmaE1;
%     tSigma2(:,:,1) = sigmaE2;
%     
%     for i = 2 : maxTdiff
%         tTheta1(:,:,i) = tTheta1(:,:,i-1) + thetaE1^(i-1);
%         tTheta2(:,:,i) = tTheta2(:,:,i-1) + thetaE2^(i-1);
%         tSigma1(:,:,i) = tTheta1(:,:,i)*sigmaE1*tTheta1(:,:,i)';
%         tSigma2(:,:,i) = tTheta2(:,:,i)*sigmaE2*tTheta2(:,:,i)';
%     end
%     
%     em = zeros(nStates, L);
%     em(1,:) = emit(x, thetaE1, tSigma1, means(:,1));
%     em(2,:) = emit(x, thetaE2, tSigma2, means(:,2));
    

    m = zeros(1,nvars);
    em = zeros(nStates, L);
    em(:,1) = 1;
    for i = 2 : L
        data1 = x(2:end,i) - [means(:,1) thetaE1]*[1; x(2:end,i-1)]; 
        data2 = x(2:end,i) - [means(:,2) thetaE2]*[1; x(2:end,i-1)]; 
        em(1,i) = mvnpdf(data1', m, sigmaE1);
        em(2,i) = mvnpdf(data2', m, sigmaE2);
        if em(1,i) < 1e-300
            em(1,i) = 1e-300;
        end
        if em(2,i) < 1e-300
            em(2,i) = 1e-300;
        end        
    end
    
% %     [f, p1, scale]=n_forward(em, a, timeDiff);
% %     b=n_backward(em, a, scale, timeDiff);
% %      % probability of each state at each position
% %     pr = f.*b;
% %     pr(:,1) = [];
% %     
% %     % get expected values
% %     logf = log(f);
% %     logb = log(b);
% %     logGE = log(em);
% %     logGTR = log(a);
% %     % f and b start at 0 so offset seq by one
% %     
% %     % this could be made a whole lot more efficient, but it's better to be
% %     % explicit
% %     
% %     if debugON
% %         L11 = testL1(em, a, timeDiff, nStates, L);
% %     end
% %     
% %     for k = 1:nStates
% %         for l = 1:nStates
% %             for i = 1:L-1
% %                 temp(i) = exp( logf(k,i) + logGTR(k,l)*timeDiff(i+1) + logGE(l,i+1) + logb(l,i+1))./scale(i+1);
% %                 A(k,l) = A(k,l) + temp(i)*timeDiff(i+1);
% %             end
% %         end
% %     end
% %     
% %     for k=1:nStates
% %         for l=1:nStates
% %             a(k,l)=A(k,l)/sum(sum(A(k,:)));
% %         end
% %     end
% %     
% %     if debugON
% %         a
% %         L12 = testL1(em, a, timeDiff, nStates, L);
% %         L12 - L11
% %         if L12 - L11 < 0
% %             disp('L1 is decreasing')
% %         end
% %     end
    
    % We use scale everything to stay out of numerical trouble. See Durbin
    % sect. 3.6
    [f, p1, scale]=n_forward(sa, em, a, timeDiff);
    b=n_backward(em, a, scale, timeDiff);
    
    % calculate improvement
    % if we go backwards, something bad has happened
    if abs((comdiff) + ((sum(p1)-loglik))) < 1e-5
        disp(['Oscillating divergence = ' num2str(comdiff)]);
        loglik = -realmax; 
        break
    end
    if iter > iter_max - 2
        disp(['maximum number of iteration reached. iter = ' num2str(iter)]);
        loglik = -realmax;
        break
    end
    
    comdiff=sum(p1)-loglik;
%     if debugON
%         comdiff
%         a
%         [thetaE1 sigmaE1; thetaE2 sigmaE2]
%     end
    if comdiff < 0
        disp(['Negative diff = ' num2str(comdiff)]);
        loglik = -realmax;
%         break
    end
    
    % if this happens, it's very bad
    if ~isfinite(p1)
        disp('Invalid prob');
        loglik = -realmax;
        break
    end
    
    % probability of each state at each position
    pr = f.*b;
    pr(:,1) = [];
    
    % get expected values
    logf = log(f);
    logb = log(b);
    logGE = log(em);
    logGTR = log(a);
    % f and b start at 0 so offset seq by one
    
    % this could be made a whole lot more efficient, but it's better to be
    % explicit
    
%     if debugON
%         L11 = testL1(sa, em, a, timeDiff, nStates, L);
%     end
    
    tempA = zeros(nStates, nStates, L);
    for k = 1:nStates
        for l = 1:nStates
            for i = 3:L
                tempA(k,l,i) = exp( logf(k,i-1) + logGTR(k,l)*timeDiff(i) + logGE(l,i) + logb(l,i))./scale(i);
                A(k,l) = A(k,l) + tempA(k,l,i)*timeDiff(i);
            end
        end
    end
    
%     if debugON
%         L31 = testL3(a, timeDiff, nStates, L, tempA);
%     end
    
    for k=1:nStates
        for l=1:nStates
            a(k,l)=A(k,l)/sum(sum(A(k,:)));
        end
    end
    sa = pr(:,1);
    
%     if debugON
%         disp([num2str(a)])
%         L12 = testL1(sa, em, a, timeDiff, nStates, L);
%         L32 = testL3(a, timeDiff, nStates, L, tempA);
% %         tempLL1 = L12-L11;       
%         tempLL1 = L32-L31;       
%         if tempLL1 < 0
%             disp(['LL1 = ' num2str(tempLL1) ' is decreasing'])
%         else
%             disp(['LL1 = ' num2str(tempLL1)])
%         end
%     end
    
    % % %     V = zeros(nStates, nvars); % we don't really need V, it's just for debugging
    % % %     for k = 1:nStates
    % % %         m = zeros(1, nvars);
    % % %         v = zeros(1, nvars);
    % % %         for l = 1:nvars
    % % %             m(l) = sum(pr(k,:) .* x(l,:))/sum(pr(k,:));
    % % %             v(l) = sum(pr(k,:) .* x(l, :).^2)/sum(pr(k, :));
    % % %         end
    % % %         M(k, :) = M(k, :) + m;
    % % %         V(k, :)=  V(k, :) + v - m.^2;
    % % %     end
    % % %
    % % %     % estimate covariance - mle estimate
    % % %     for k = 1:nStates
    % % %         for l = 1:nvars
    % % %             for l2 = 1:nvars
    % % %                 s = sum(pr(k, :) .* ((x(l,:) - M(k, l)) .* ((x(l2,:) - M(k, l2)))))/sum(pr(k,:));
    % % %                 if k == 1
    % % %                     COV1(l, l2) = s;
    % % %                 else
    % % %                     COV2(l, l2) = s;
    % % %                 end
    % % %             end
    % % %         end
    % % %     end
    % % %
    % % %     for k = 1:nStates
    % % %         for l = 1:nvars
    % % %             means(k, l) = M(k, l);
    % % %         end
    % % %     end
    % % %
    % % %     cov1 = COV1;
    % % %     cov2 = COV2;
     
%     if debugON
%         L21 = testL2(x, thetaE1, sigmaE1, thetaE2, sigmaE2, a, sa, means);
%         L41 = testL4(x, thetaE1, sigmaE1, thetaE2, sigmaE2, a, pr, means);
%         tic
%     end
    
    for k = 1 : nStates
%         if maxTdiff == 1
%             tempYX = zeros(nvars,nvars+1);
%             tempXX = zeros(nvars+1,nvars+1);
%             for i = 2 : L
%                 tempYX = tempYX + pr(k,i-1)*x(2:end, i)*[1; x(2:end, i-1)]';
%                 tempXX = tempXX + pr(k,i-1)*[1; x(2:end, i-1)]*[1; x(2:end, i-1)]';
%             end
%             cthetaE = tempYX/(tempXX);
%             
%             thetaE = cthetaE(:,2:end);
%             means(:,k) = cthetaE(:,1);
%             
%             sigmaE = comSigmaE(cthetaE,x(2:end,1:end-1),x(2:end,2:end),timeDiff(2:end),pr(k,:));
%             
%         else
            if k == 1
                preThetaE = thetaE1;
                preSigmaE = sigmaE1; 
            else
                preThetaE = thetaE2;
                preSigmaE = sigmaE2;            
            end

%             LB = -ones(nvars,nvars+1) + [zeros(nvars,1) 1e-4*ones(nvars)];
%             UB =  ones(nvars,nvars+1) - [zeros(nvars,1) 1e-4*ones(nvars)];
%             
%             totalUpdate = 10;
%             iterTol = 1e-4;
%             count(k) = 0;
%             exitFlagSigma = 1;
%             exitflag = 1;
%             while totalUpdate > iterTol || exitflag == 0 || exitFlagSigma == 0
                tolEig = 1e-10;
                if abs(eigs(preThetaE,1)) > 1-tolEig
                    preThetaE = preThetaE/norm(preThetaE,1);
                end
                precTheta = [means(:,k) preThetaE];

                %                 [thetaE,fval,exitflag] = fminsearch(@(theta) comErrorVec2(theta,x(2:5,1:end-1),x(2:5,2:end),timeDiff(2:end),sigmaE,pr(k,:)), preThetaE, optimset('Display',displayOption));
                %                 [thetaE,fval,exitflag] = fminsearchbnd(@(theta) comErrorVec2(theta,x(2:5,1:end-1),x(2:5,2:end),timeDiff(2:end),sigmaE,pr(k,:)), preThetaE, LB, UB, optimset('Display',displayOption));
%                 [cthetaE,fval,exitflag] = fminsearchbnd(@(ctheta) comErrorVec2(ctheta,x(2:end,1:end-1),x(2:end,2:end),timeDiff(2:end),preSigmaE,pr(k,:)), precTheta, LB, UB, optimset('Display',displayOption));
%                 [cthetaE,fval,exitflag] = fminsearchcon(@(ctheta) comErrorVec2(ctheta,x(2:end,1:end-1),x(2:end,2:end),timeDiff(2:end),preSigmaE,pr(k,:)), precTheta, [],[],[],[], @(ctheta) norm(ctheta(:,2:end),1)-1, optimset('Display',displayOption));
                [cthetaE,fval,exitflag] = fminsearchcon(@(ctheta) comErrorVec2(ctheta,x(2:end,1:end-1),x(2:end,2:end),timeDiff(2:end),preSigmaE,pr(k,:)), precTheta, [],[],[],[], @(ctheta) abs(eigs(ctheta(:,2:end),1))-(1-tolEig), optimset('Display',displayOption));
                
                thetaE = cthetaE(:,2:end);
                means(:,k) = cthetaE(:,1);
                
                %             sigmaE = comSigmaE(thetaE,x(2:5,1:end-1),x(2:5,2:end),timeDiff(2:end),pr(k,:));
                %             sigmaE = comSigmaE(cthetaE,x(2:end,1:end-1),x(2:end,2:end),timeDiff(2:end),pr(k,:));
                
                prC = pr(k,:);
                diffYX = x(2:end,2:end) - cthetaE*[ones(1,L-1); x(2:end,1:end-1)];
                tempSigma = zeros(nvars,nvars);
                for i = 1 : L-1
                    temp = diffYX(:,i);
                    tempSigma = tempSigma + prC(i)*(temp*temp');
                end
                sigmaE = tempSigma/sum(prC);
                
%                 tTheta = zeros(nvars, nvars, maxTdiff);
%                 tSigma = zeros(nvars, nvars, maxTdiff);
%                 tTheta(:,:,1) = eye(nvars);
%                 tSigma(:,:,1) = sigmaE;
%                 errSigma = zeros(maxTdiff,1);
%                 [~, errSigma(1)] = cholcov(tSigma(:,:,1), 0);
%                 for i = 2 : maxTdiff
%                     tTheta(:,:,i) = tTheta(:,:,i-1) + thetaE^(i-1);
%                     tSigma(:,:,i) = tTheta(:,:,i)*sigmaE*tTheta(:,:,i)';
%                     [~, errSigma(i)] = cholcov(tSigma(:,:,i), 0);
%                     if errSigma(i) ~= 0
%                         exitFlagSigma = 0;
%                     end
%                 end
%                 
%                 diffTheta = thetaE-preThetaE;
%                 diffSigma = sigmaE-preSigmaE;
%                 temp1 = sum(sum(diffTheta.*diffTheta + diffSigma.*diffSigma));
%                 temp2 = sum(sum(thetaE.*thetaE + sigmaE.*sigmaE));
%                 totalUpdate = temp1/temp2*100;
%                 
%                 count(k) = count(k)+ 1;
%                 preThetaE = thetaE;
%                 preSigmaE = sigmaE;                
%             end
     
            
%         end  
        
        if k == 1
            thetaE1 = thetaE; sigmaE1 = sigmaE;
        else
            thetaE2 = thetaE; sigmaE2 = sigmaE;
        end

    end
%     if debugON     
%         L22 = testL2(x, thetaE1, sigmaE1, thetaE2, sigmaE2, a, sa, means);
%         L42 = testL4(x, thetaE1, sigmaE1, thetaE2, sigmaE2, a, pr, means);
%         toc
% %         tempLL2 - L22 - L21;       
%         tempLL2 = L42 - L41;       
%         if tempLL2 < 0
%             disp(['LL2 = ' num2str(tempLL2) ' is decreasing'])
%         else
%             disp(['LL2 = ' num2str(tempLL2)])
%         end
%     end
        
    % make sure the covariances are legit
    [~, err] = cholcov(sigmaE1, 0);
    if err ~= 0
        disp('sigmaE1 is not positive definite')
        loglik = -realmax;
        break;
    end
    
    [~, err] = cholcov(sigmaE2, 0);
    if err ~= 0
        disp('sigmaE2 is not positive definite')
        loglik = -realmax;
        break;
    end
    
    % update log likelihood and do it again
    loglik=sum(p1);
    if ~isfinite(loglik)
        disp('Invalid LL');
        loglik = -realmax;
        break;
    end
    
    iter = iter + 1;
    disp(['LL = ' num2str(loglik) ...
        ' diff = ' num2str(comdiff) ...
        ' iter = ' int2str(iter)])
    
%     if debugON
%     else
%         save(sprintf('%d_iter%d.mat',fnum,iter))
%     end

    x = data_augmentation(orgX, thetaE1, thetaE2, sigmaE1, sigmaE2, means, pr);
end