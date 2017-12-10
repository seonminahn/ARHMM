function p = testL1(sa, em, a, timeDiff, nStates, L)

[f, p1, scale]=n_forward(sa, em, a, timeDiff);
b=n_backward(em, a, scale, timeDiff);
pr = f.*b;
pr(:,1) = [];
logf = log(f);
logb = log(b);
logGE = log(em);
logGTR = log(a);
testL = zeros(nStates,nStates);
for k = 1:nStates
    for l = 1:nStates
        for i = 3:L
            temp = exp( logf(k,i-1) + logGTR(k,l)*timeDiff(i) + logGE(l,i) + logb(l,i))./scale(i);
            testL(k,l) = testL(k,l) + temp*logGTR(k,l)*timeDiff(i);
        end
    end
end
p = sum(sum(testL));
