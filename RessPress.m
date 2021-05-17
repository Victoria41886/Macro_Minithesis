function [RessPress]=RessPress(E,PARAM,SS)

%At SS
%X_0=Ess^(1/PARAM(10))

%For now let's just set it to 8
X_0=SS.x;

%X(1)=X_0-X_0.^PARAM(10)+E(1);
X(1)=X_0*(1-PARAM(10))+E(1);

for ii=2:size(E,1)
    %X(ii)=X(ii-1)-X(ii-1).^PARAM(10)+E(ii-1);
    X(ii)=X(ii-1)*(1-PARAM(10))+E(ii-1);
    %Pressure on ressources at the beginning of each period
    
    GAMMA(ii)=PARAM(7)+PARAM(8)*(X(ii)+E(ii))^PARAM(9)/1e7;
    %Cost of ressource extraction at each period
end
GAMMA(1)= PARAM(7)+PARAM(8)*(X(1)+E(1))^PARAM(9)/1e7;

RessPress(:,1)=X;
RessPress(:,2)=GAMMA;

end
