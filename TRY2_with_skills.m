%To use this File SolvoEco and SolvoTimeEco need to be modified to take N
%as an input argument and not calculate within these function files
%
%Later the fsolve for the whole economy however, runs forever without
%getting to any solution

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



%% Including High skilled workers into the model 


%New steady state for high skilled workers
GUESS=[39,5];
A=1;                          % The A we use here should be the same as the initial A for the climate change process
L=zeros(20,2)
L(1:15,1)=1                     % the first column represents L for Low skills, the 2nd col for high skills
L(2:15,2)=1.1
N= sum(, 'all')

S=fsolve(@(guess) SolvEco(guess,L(:,2),N,A,PARAM),GUESS);
S=abs(S);

%Let's do some verif
verif=SolvEco(S,L(:,2),N,A,PARAM);
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

AO=AgentOptim(ECO(2)*ones(20,1),ECO(3),A,L(:,2),PARAM,Ini);
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

%Dynamic whole economy 
Atime=cumprod(1.001*ones(100,1));

GUESS=[ones(100,1)*39,ones(100,1)*5];

SolveTimeEco(GUESS,L(:,2),N,PARAM,SS);

options.MaxFunEvals = 10000 
Solution=abs(fsolve(@(x) SolveTimeEco(x,L(:,2),N,PARAM,SS),GUESS,options));


RessPress(Solution(:,2),PARAM,SS)

%Let's do some verif
verif=SolveTimeEco(Solution,L,PARAM,SS)

ECO=Economy(Solution(:,1),N,Solution(:,2),Atime,PARAM);
AllAgents=AGENTS(ECO(:,2),ECO(:,3),Atime,L(:,2),PARAM,SS);

for jj=1:3:100
plot(1:100,AllAgents(:,jj,1),"b-",1:100,AllAgents(:,jj,2),"r:")
hold on
end


