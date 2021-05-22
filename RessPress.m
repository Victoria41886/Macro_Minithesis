function [RessPress]=RessPress(E,PARAM,SS)

%AT SS
X_0=SS.x;

X(1)=X_0*(1-PARAM(10))+E(1);
GAMMA(1)= PARAM(7)+PARAM(8)*((X(1)+E(1))/PARAM(11))^PARAM(9);

for ii=2:size(E,1)
    %X(ii)=X(ii-1)-X(ii-1).^PARAM(10)+E(ii-1);
    X(ii)=X(ii-1)*(1-PARAM(10))+E(ii-1);
    %Pressure on ressources at the beginning of each period
    
    GAMMA(ii)=PARAM(7)+PARAM(8)*((X(ii)+E(ii))/PARAM(11))^PARAM(9);
    %Cost of ressource extraction at each period
end

RessPress(:,1)=X;
RessPress(:,2)=GAMMA;

end
