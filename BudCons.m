function [BudCons]=BudCons(C1,R,S,PARAM)

CR=cumprod((1+R-PARAM(5)))/(1+R(1)-PARAM(5));
%Cumulated interest rate net of depreciation
DF=1./CR
BudCons=sum((ConsOpt(C1,R,PARAM)-S).*DF);
%If the agent budget constraint is verified then BudCons=0
end
