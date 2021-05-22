function [ConsOpt]=ConsOpt(C1,R,PARAM)

CR=cumprod((1+R-PARAM(5)))/(1+R(1)-PARAM(5));
ConsOpt=C1*PARAM(4).^([0:size(R,1)-1].'/PARAM(6)).*CR.^(1/PARAM(6));

%From 3
%If C(1)=C1 then the vector C (based on utility maximization) is ConsOpt
end
