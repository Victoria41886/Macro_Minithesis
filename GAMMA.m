function [GAMMA]=GAMMA(X,R,PARAM)
%In part 1/ GAMMA=.4

%In part 2/
GAMMA=PARAM(7)+PARAM(8)*((X+R)/PARAM(11))^PARAM(9);
end

