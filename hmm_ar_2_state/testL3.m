% function p = testL3(a, timeDiff, nStates, L, logf, logb, logGE, scale)
function p = testL3(a, timeDiff, nStates, L, tempA)
logGTR = log(a);

testL = zeros(nStates,nStates);
for k = 1:nStates
    for l = 1:nStates
        for i = 3:L
%             temp = exp( logf(k,i) + logGTR(k,l)*timeDiff(i+1) + logGE(l,i+1) + logb(l,i+1))./scale(i+1);
            testL(k,l) = testL(k,l) + tempA(k,l,i)*logGTR(k,l)*timeDiff(i);
        end
    end
end
p = sum(sum(testL));