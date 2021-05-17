clear

%PARAMETERS
alpha=.33;
PARAM(1)=alpha;

theta=0.4;
PARAM(2)=theta;

eff=.5;
PARAM(3)=eff;

B=.98; % beta
PARAM(4)=B;

dep=0; %tau? 
PARAM(5)=dep;

sigma=2;
PARAM(6)=sigma;

Er_1=.3;
Er_2=1;
Er_3=4;
Er_4=0.01;

PARAM(7)=Er_1;
PARAM(8)=Er_2;
PARAM(9)=Er_3;
PARAM(10)=Er_4;
%Except the last one which has to be below one, these numbers could not be
%more random


L=[ones(15,1);zeros(5,1)];
N=sum(L);


%% 1/ First climate change has no impact and ressources are extracted at a constant cost

%{

r=.1303;
w=1;

AO=AgentOptim(r,w,A,L,PARAM);

figure(1)
plot(1:20,AO(:,1),"b--",1:20,AO(:,2),"r-")
%Looks good
%}

%% Steady State

GUESS=[39,5];
A=1;

S=fsolve(@(guess) SolvEco(guess,L,A,PARAM),GUESS);
S=abs(S);

%Let's do some verif
verif=SolvEco(S,L,A,PARAM);
Economy(S(1),N,S(2),A,PARAM);

%AMAZING

%Let's do a few plots
SS.k=S(1);    
SS.e=S(2);

%At SS
%SS.x=SS.e^(1/PARAM(10));
SS.x=SS.e/Er_4

ECO=Economy(SS.k,N,SS.e,A,PARAM);


%We assume that agents are born with 0 asset
Ini.asset=0

AO=AgentOptim(ECO(2)*ones(20,1),ECO(3),A,L,PARAM,Ini);
SS.asset=AO(:,2)
figure(2)
plot(1:20,AO(:,1),"b--",1:20,AO(:,2),"r-")

Yss=F(SS.k,N,SS.e,A,PARAM)

%Verification that the income shares are right
A*N*ECO(1)/Yss
SS.k*ECO(2)/Yss
SS.e*ECO(3)/Yss
%Niceeee


%%
%Now with vectors 


R=cumprod(1.01*ones(20,1))-1

W=unidrnd(100,20,1)/100
A=unidrnd(100,20,1)/100+.5


clear AgentOptim
AO=AgentOptim(R,W,A,L,PARAM,Ini);

figure(1)
plot(1:20,AO(:,1),"b--",1:20,AO(:,2),"r-")
legend("Consumption","Assets")


%% Now for all the agents

%Wages and interest rates are random
Wtime=unidrnd(100,100,1)/100;
Rtime=unidrnd(100,100,1)/500;

%Wtime=ones(100,1)*.3;
Rtime=ones(100,1)/100;

Atime=cumprod(1.001*ones(100,1));

AllAgents=AGENTS(Rtime,Wtime,Atime,L,PARAM,SS)

Cagents=AllAgents(:,:,1);
Aagents=AllAgents(:,:,2);

figure(2)
for jj=1:4:81
plot(1:100,AllAgents(:,jj,1),"b-",1:100,AllAgents(:,jj,2),"r:")
hold on
end

figure(3)
jj=2
plot(1:100,AllAgents(:,jj,1),"b--",1:100,AllAgents(:,jj,2),"r-")


%% 4/ And now it's linked to the whole economy


%Imagine high skill worker: L changes, w changes, initial asset changes
L(1)=0

Atime=cumprod(1.001*ones(100,1));

GUESS=[ones(100,1)*39,ones(100,1)*5];

SolveTimeEco(GUESS,L,PARAM,SS);

options.MaxFunEvals = 50000 

Solution=abs(fsolve(@(x) SolveTimeEco(x,L,PARAM,SS),GUESS, options))


RessPress(Solution(:,2),PARAM,SS)

%Let's do some verif
verif=SolveTimeEco(Solution,L,PARAM,SS)

ECO=Economy(Solution(:,1),N,Solution(:,2),Atime,PARAM);
AllAgents=AGENTS(ECO(:,2),ECO(:,3),Atime,L,PARAM,SS);

for jj=1:3:100
plot(1:100,AllAgents(:,jj,1),"b-",1:100,AllAgents(:,jj,2),"r:")
hold on
end


%%


%{

%% 3/ Climate Change for 10 periods
E=zeros(10,1)
m_t= zeros(10,1)
G = zeros(10,1)
O = zeros(10,1)
F = zeros(10,1)
A= zeros(10,1)
A_t=zeros(10,1)

% Calibration values
cc1 = 0.9301
cc2= 0.9846
cc3= 1.42/2.980
cc4= 0.9112
cc5= 2.980 % I take the value for sc. 2
cc6=  0.5     % missing in paper, this value is just a placeholder
cc7= 0.998
cc8= 0.0150 % damage from a 3 degree increase in Sc. 2
cc9= 0.028
A_0= 2.5 % not sur, check again

%Initial Values
m_t(1)= 590 % preindustrial carbon level
G(1)=0 % initial surface temperature change
O(1)=0 % initial ocean temperature change
E(1)= PARAM(3)*Ess
A_t(1)= A_0 * exp((cc9/cc2)*(1-exp(-cc2*1)))

for i = 2:1:10
    E(i) = PARAM(3)*Ess
    m_t(i) = m_t(1) + cc1*E(i-1) + cc2*(m_t(i-1) - m_t(1))
    F(i)= log(m_t(i)/m_t(1))/log(2) + cc3
    G_t(i) = cc4*G(i-1) + cc5*F(i)+ cc6*O(i-1)
    O(i)= cc7*O(i-1)+(1-cc7)*G(i-1)
    A_t(i)= A_0 * exp((cc9/cc2)*(1-exp(-cc2*i)))
    A(i)=(1+cc8*G_t(i).^2).^(-1/(1-PARAM(1)))*A_t(i)
end


%% Trying to include Climate chnage into the model


% My main concern right now is this constantly changing A, which leads to
% constantly changing wages, I do not know to handle this, since this far
% we have only handled fixed wages, which were known in advanve or a one
% time shock and convergence to a new equilibirium. But calculating a new
% equilibrium for each time period does not make a lot of sense either.



%%Here I first try to incorporate it into agents optimization,
%%while R is constant at the intitial steady state level -> shoudl we
%%endogeneize Ressource use?



% This section can be basically skipped. I was too sentimental to delete it. The idea was to
% calculate a new optimized consumption for each period according to the
% current wages and interest-rates.
%However, we can probably not continue using the formula for ideal
%consumption any more which I did and moreover, you only compute the
%consumption for one generation for one period, sowe do not get any
%sensible value for K. This is as if only one generation would exist at a
%time.
ECO=repmat(Economy(Kss,N,Ess,1,PARAM),10,1) % This is just some initial guess to establish the ECO matrix, be aware that the algorithm updates only lines 2 to 10
ECO(1,:)=Economy(Kss,N,Ess,1,PARAM) % correcting period 1 
K=zeros(10,1)
K(1,1)=Kss
AgentOptim=zeros(10,2)
x=zeros(10,1)
AO= zeros(10,1)
for i = 2:1:10
    K(i)=abs(K(i-1,1)); % we use last periods K to compute ECO later on, at the end of the step we will replcae k(i) with the actual k(i) of this period, computed as the sum of ai
    R=abs(Ess);
        
    N=sum(L);
    ECO(i,:)=Economy(K(i,1),N,R,A(i,1),PARAM)
    
     x(i,1)=sum(ECO(i,1).*(1+ECO(i,2)-PARAM(5)).^-(i))/sum((1+ECO(i,2)-PARAM(5)).^-(i).*(PARAM(4)*(1+ECO(i,2)-PARAM(5))).^(i/PARAM(6)))
   %optimal consumption for period i 
    AgentOptim(i,1)= x(i,1)*(PARAM(4)*(1+ECO(i,2)-PARAM(5))).^(i/PARAM(6))
    AgentOptim(i,2)= ECO(i-1,1)-AgentOptim(i-1,1)+AgentOptim(i-1,2)*(1+ECO(i,2)-PARAM(5))

    %AO(i,1)=AgentOptim2(ECO(i,2),ECO(i,1),A(i),L,PARAM) 
    % the issue with agent optim here, because agent optim takes a value for 
    % for wage and r and then calculates c and a for all periods, we need
    % to update our wage and r in each period, with AO we get an entire 
    
    K(i)=sum(AgentOptim(:,2))
end


% In this step i do sort of the same as above but do not calculate only one new
% optimal Consumption for each period, but basically an optimal consumption
% for each generation according to the current wage and interest rate. 

ECO=repmat(Economy(Kss,N,Ess,1,PARAM),10,1) % This is just some initial guess to establish the ECO matrix, be aware that the algorithm updates only lines 2 to 10
ECO(1,:)=Economy(Kss,N,Ess,1,PARAM) % correcting period 1 
K=zeros(10,1)
K(1,1)=Kss

for i = 2:1:10
    K(i)=abs(K(i-1,1)); % we use last periods K to compute ECO later on, at the end of the step we will replcae k(i) with the actual k(i) of this period, computed as the sum of ai
    R=abs(Ess);
        
    N=sum(L);
    ECO(i,:)=Economy(K(i,1),N,R,A(i,1),PARAM)
    AO = AgentOptim(ECO(i,2),ECO(i,1),A(i,1),L,PARAM) 
       
    K(i)=sum(AO(:,2))
end





%}
