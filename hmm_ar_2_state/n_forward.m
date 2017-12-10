function [fs p s]=n_forward(sa, em,a,timeDiff)
% em(state, :) - emission probability of emitting mean, std when in state
% a(state1, state2) = probability of transition from state1 to state 2
% returns:
%   fs - scaled forward matrix
%   p - log P(data| parameters)
%   s - scale factor for each position, 1 X length of data

nStates=size(em,1);
L=size(em, 2);

fs = zeros(nStates,L);
fs(:,1) = 1;
% fs(:,1) = 0;
fs(:,2) = em(:,2).*sa;
s = zeros(1,L);
s(1) = 1;
s(2) = sum(fs(:,2));
fs(:,2) =  fs(:,2)./s(2);
for count = 3:L
    aCurrent = a^timeDiff(count);
    for state = 1:nStates
        fs(state,count) = em(state, count) .* (sum(fs(:,count-1) .*aCurrent(:,state)));
    end
    % scale factor normalizes sum(fs,count) to be 1. 
    s(count) =  sum(fs(:,count));
    fs(:,count) =  fs(:,count)./s(count);
end

p = sum(log(s));
end



% % % function [fs p s]=n_forward(em,a)
% % % % em(state, :) - emission probability of emitting mean, std when in state
% % % % a(state1, state2) = probability of transition from state1 to state 2
% % % % returns:
% % % %   fs - scaled forward matrix
% % % %   p - log P(data| parameters)
% % % %   s - scale factor for each position, 1 X length of data
% % % 
% % % nStates=size(em,1);
% % % L=size(em, 2);
% % % 
% % % fs = zeros(nStates,L);
% % % fs(1,1) = 1;  % assume that we start in state 1.
% % % s = zeros(1,L);
% % % s(1) = 1;
% % % for count = 2:L
% % %     for state = 1:nStates
% % %         fs(state,count) = em(state, count) .* (sum(fs(:,count-1) .*a(:,state)));
% % %     end
% % %     % scale factor normalizes sum(fs,count) to be 1. 
% % %     s(count) =  sum(fs(:,count));
% % %     fs(:,count) =  fs(:,count)./s(count);
% % % end
% % % 
% % % p = sum(log(s));
% % % end


