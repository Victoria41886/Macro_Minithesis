% to use this file it is necessary to use the skill adjusted function files
% Economy, change_skill, SolvoEco and SolvoTimeEco

clear

%PARAMETERS
alpha=.33; %similar in paper
PARAM(1)=alpha;

theta=0.05; %from paper
PARAM(2)=theta;

eff=.5; % what is this parameter???
PARAM(3)=eff;

B=.98; % beta, similar in paper
PARAM(4)=B;

dep=0.1; %delta 
PARAM(5)=dep;

sigma=0.97;% from paper
PARAM(6)=sigma;

%Er_1=113; These would be the values from the paper
%Er_2=700;
%Er_3=4;

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

GUESS=[39,5];
A=1.5;     % I think this A should be the same as our initial A in the climate change function,
           % there A=2.5 originally, but now it is set to 1.5 
skill=1.1

L_l=0.5*[ones(15,1);zeros(5,1)]; %low skilled workers
L_h=zeros(20,1)                  %high skilled workers
L_h(2:15)=0.5
N=sum(L_h).*skill+sum(L_l);

S=fsolve(@(guess) SolvEco_skill(guess,L_h,N,A,PARAM,skill),GUESS);
S=abs(S);

%Let's do some verif
%verif=SolvEco(S,L,N,A,PARAM,skill);
%Economy(S(1),N,S(2),A,PARAM,skill);

%AMAZING

%Let's do a few plots
SS.k=S(1);    
SS.e=S(2);

%At SS
%SS.x=SS.e^(1/PARAM(10));
SS.x=SS.e/Er_4

ECO=Economy_skill(SS.k,N,SS.e,A,PARAM,skill);


%We assume that agents are born with 0 asset
Ini.asset=0

AO=AgentOptim(ECO(2)*ones(20,1),ECO(1),A,L_h,PARAM,Ini);  %HERE WAS A MISTAKE; ECO1=W 
SS.asset=AO(:,2)
figure(2)
plot(1:20,AO(:,1),"b--",1:20,AO(:,2),"r-")


%% 

%Dynamic whole economy 

Atime=cumprod(1.001*ones(100,1));

GUESS=[ones(100,1)*39,ones(100,1)*5];

SolveTimeEco_skill(GUESS,L_h,N,PARAM,SS,skill);

options.MaxFunEvals = 15000 
Solution=abs(fsolve(@(x) SolveTimeEco_skill(x,L_h,N,PARAM,SS,skill),GUESS, options))


RessPress(Solution(:,2),PARAM,SS)

%Let's do some verif
%verif=SolveTimeEco(Solution,L,PARAM,SS)

ECO=Economy_skill(Solution(:,1),N,Solution(:,2),Atime,PARAM,skill);
AllAgents=AGENTS(ECO(:,2),ECO(:,3),Atime,L_h,PARAM,SS);

for jj=1:3:100
plot(1:100,AllAgents(:,jj,1),"b-",1:100,AllAgents(:,jj,2),"r:")
hold on
end

%% Same for low skilled workers
% always clear before running this whole economy one time

clear

%PARAMETERS
alpha=.33;
PARAM(1)=alpha;

theta=0.05;
PARAM(2)=theta;

eff=.5; % what is this parameter???
PARAM(3)=eff;

B=.98; % beta
PARAM(4)=B;

dep=0.1; %delta 
PARAM(5)=dep;

sigma=0.97;
PARAM(6)=sigma;

%Er_1=113; These would be the values from the paper
%Er_2=700;
%Er_3=4;

Er_1=.3;
Er_2=1;
Er_3=4;
Er_4=0.01;

PARAM(7)=Er_1;
PARAM(8)=Er_2;
PARAM(9)=Er_3;
PARAM(10)=Er_4;

skill=1;
GUESS=[39,5];
A=1.5; 

L_l=0.5*[ones(15,1);zeros(5,1)]; %low skilled workers
L_h=zeros(20,1)                  %high skilled workers
L_h(2:15)=0.5
N=sum(L_h).*skill+sum(L_l);

%Steady state
S=fsolve(@(guess) SolvEco_skill(guess,L_l,N,A,PARAM,skill),GUESS);
S=abs(S);

SS.k=S(1);    
SS.e=S(2);
SS.x=SS.e/Er_4

ECO=Economy_skill(SS.k,N,SS.e,A,PARAM,skill);

Ini.asset=0
AO=AgentOptim(ECO(2)*ones(20,1),ECO(1),A,L_l,PARAM,Ini);
SS.asset=AO(:,2)
figure(2)
plot(1:20,AO(:,1),"b--",1:20,AO(:,2),"r-")


%whole economy
Atime=cumprod(1.001*ones(100,1));

GUESS=[ones(100,1)*39,ones(100,1)*5];

SolveTimeEco_skill(GUESS,L_l,N,PARAM,SS,skill);

options.MaxFunEvals = 15000 

Solution=abs(fsolve(@(x) SolveTimeEco_skill(x,L_l,N,PARAM,SS,skill),GUESS, options))


RessPress(Solution(:,2),PARAM,SS)

%Let's do some verif
%verif=SolveTimeEco(Solution,L,PARAM,SS)

ECO=Economy_skill(Solution(:,1),N,Solution(:,2),Atime,PARAM,skill);
AllAgents=AGENTS(ECO(:,2),ECO(:,3),Atime,L_l,PARAM,SS);

for jj=1:3:100
plot(1:100,AllAgents(:,jj,1),"b-",1:100,AllAgents(:,jj,2),"r:")
hold on
end
