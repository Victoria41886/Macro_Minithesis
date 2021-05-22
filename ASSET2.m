function [ASSET]=ASSET2(S,Cpath,R,PARAM,Asset_ss)

ASSET(1)=Asset_ss(21-size(R,1));

for ii=2:size(R,1)
    ASSET(ii)=S(ii-1)-Cpath(ii-1)+ASSET(ii-1).*(1+R(ii-1)-PARAM(5));
    %from the budget constraint
end

end
