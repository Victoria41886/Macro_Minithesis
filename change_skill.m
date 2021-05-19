function [A] = change_skill(R,PARAM,skill)
%UNTITLED This has at the moment two columns, but it could also give a
%value depending on skill


%Summary of this function goes here
%   This function depicts a slowly decreasing productivity if Emissiom
%   remain constant. The initial value of A, calibrated by the authors,  however is very high, such
%   that we still have massive growth.

E=zeros(100,1);
m_t= zeros(100,1);
G = zeros(100,1);
O = zeros(100,1);
F = zeros(100,1);
A= zeros(100,1);
A_t=zeros(100,1);

% Calibration values
cc1 = 0.9301;
cc2= 0.9846;
cc3= 1.42/2.980;
cc4= 0.9112;
cc5= 2.980; % I take the value for sc. 2
cc6=  0.5 ;    % missing in paper, this value is just a placeholder
cc7= 0.998;
cc8= 0.0150; % damage from a 3 degree increase in Sc. 2
cc9= 0.028;
A_0= 2.5; % in the paper the initial productivitz is 2.5, which seems pretty high


% now we introduce a skill specific paramter to have a different evolution
% of A for high and low skilled workers

if skill==1
    cc10=1
else
    cc10= 0.5 
end

%Initial Values
m_t(1)= 590; % preindustrial carbon level
G(1)=0; % initial surface temperature change
O(1)=0; % initial ocean temperature change
E(1)= PARAM(3).*R(1);
A_t(1)= A_0 * exp((cc9/cc2)*(1-exp(-cc2*1)));


for i = 2:1:100
    E(i) = PARAM(3).*R(i);
    m_t(i) = m_t(1) + cc1*E(i-1) + cc2*(m_t(i-1) - m_t(1));
    F(i)= log(m_t(i)/m_t(1))/log(2) + cc3;
    G_t(i) = cc4*G(i-1) + cc5*F(i)+ cc6*O(i-1);
    O(i)= cc7*O(i-1)+(1-cc7)*G(i-1);
    A_t(i)= A_0 * exp((cc9/cc2)*(1-exp(-cc2*i)));
    A(i,1)=(1+cc10*cc8*G_t(i).^2).^(-1/(1-PARAM(1)))*A_t(i);
end

A(1,:)=A_0


end

