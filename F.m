function [F]=F(K,AN,E,PARAM)
F=K.^PARAM(1).*(AN).^(1-PARAM(1)-PARAM(2)).*(PARAM(3)*E).^PARAM(2);
%BASED ON (1)
end