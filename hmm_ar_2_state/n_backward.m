function bs=n_backward(em, a, scale, timeDiff)
% em(state, :) - emission probability of emitting mean, std when in state
% a(state1, state2) = probability of transition from state1 to state 2
% returns:
%   bs - scaled backward matrix

[nStates, L] = size(em);

bs = ones(nStates,L);
for count = L-1:-1:2
    aCurrent = a^timeDiff(count+1);

    for state = 1:nStates
      bs(state,count) = (1/scale(count+1)) * sum( aCurrent(state,:)'.* bs(:,count+1) .* em(:,count+1));       
    end
end

end

% % % function bs=n_backward(em,a, scale)
% % % % em(state, :) - emission probability of emitting mean, std when in state
% % % % a(state1, state2) = probability of transition from state1 to state 2
% % % % returns:
% % % %   bs - scaled backward matrix
% % % 
% % % nStates=size(em,1);
% % % L=size(em, 2);
% % % 
% % % bs = ones(nStates,L);
% % % for count = L-1:-1:1
% % %     for state = 1:nStates
% % %       bs(state,count) = (1/scale(count+1)) * sum( a(state,:)'.* bs(:,count+1) .* em(:,count+1));       
% % %     end
% % % end
% % % 
% % % end