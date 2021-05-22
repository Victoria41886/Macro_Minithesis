function [A] = change(R,PARAM,Tg)
%UNTITLED Summary of this function goes here
%   This function depicts a slowly decreasing productivity if Emissiom
%   remain constant. The initial value of A, calibrated by the authors,  however is very high, such
%   that we still have massive growth.
%R=5*ones(100,1)
%Each period correspond to 4 years
R=repelem(R,4)/4
Tg=repelem(Tg,4)/4

SIZE=size(R,1);

E=zeros(SIZE,1);
m_t= zeros(SIZE,1);
G = zeros(SIZE,1);
O = zeros(SIZE,1);
F = zeros(SIZE,1);
Al= zeros(SIZE,1);
A_t=zeros(SIZE,1);

% Calibration values
cc1 = 0.9301;
cc2= 0.9846;
cc3= 1.42/2.980;
cc4= 0.9112;
cc5= 2.980; % I take the value for sc. 2
cc6=1-cc4 ;    % missing in paper, I do as for O
cc7= 0.998;
cc8= 0.0150; % damage from a 3 degree increase in Sc. 2
cc9= 0.028;
cc10= 0.7; % high skilled worker factor
A_0= 2.5; % in the paper the initial productivitz is 2.5, which seems pretty high



%Initial Values
m_t(1)= 590; % preindustrial carbon level
G(1)=0; % initial surface temperature change
O(1)=0; % initial ocean temperature change
E(1)= PARAM(3).*R(1);
A_t(1)= A_0 * exp((cc9/cc2)*(1-exp(-cc2*1)));


for i = 2:SIZE
    E(i) = PARAM(3).*R(i);
    m_t(i) = m_t(1) + cc1*E(i-1) + cc2*(m_t(i-1) - m_t(1));
    F(i)= log(m_t(i)/m_t(1))/log(2) + cc3;
    G(i) = cc4*G(i-1) + cc5*F(i)+ cc6*O(i-1);
    O(i)= cc7*O(i-1)+(1-cc7)*G(i-1);
    A_t(i)= A_0 * exp((cc9/cc2)*(1-exp(-cc2*i)));
    Al(i,1)=((Tg(i-1))^(1/8))*(1+cc8/9*G(i).^2).^(-1/(1-PARAM(1)))*A_t(i);%low skilled
    Al(i,2)=((Tg(i-1))^(1/8))*(1+cc10*cc8/9*G(i).^2).^(-1/(1-PARAM(1)))*A_t(i);%high skilled
end

Al(1,:)=A_0

A=Al;
A=Al(1:4:SIZE,:)
end

