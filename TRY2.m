clear

%PARAMETERS
alpha=.38;
PARAM(1)=alpha;

theta=0.05;
PARAM(2)=theta;

eff=.27;
PARAM(3)=eff;

B=.96; % beta
PARAM(4)=B;

dep=0.045; 
PARAM(5)=dep;

sigma=.97;
PARAM(6)=sigma;

Er_1=113;
Er_2=700;
Er_3=4;
Er_4=.01;

PARAM(7)=Er_1;
PARAM(8)=Er_2;
PARAM(9)=Er_3;
PARAM(10)=Er_4;

PARAM(11)=100

%THIS WORKS BETTER WITH THESE PARAMATERS
PARAM(7)=.3;
PARAM(8)=1;
PARAM(10)=Er_4;




%Except the last one which has to be below one, these numbers could not be
%more random

%TotalPopulation =7.5e9     *7.5e9/20
GenSize=7.5e9/20;
GenSize=1;
L=[ones(15,1)*GenSize;zeros(5,1)];
N=sum(L);

%L as divided in High and low skilled workers
L= ones(20,2);
L(16:20,:)=0;
L(1,2)=0;
L(2:15,2)=1.2;
N=sum(L, 'all');


%% INITIAL Steady State


guess=[100,100];
%If it does not work just increase the guess
A=2.5;

S=fsolve(@(guess) SolvEco(guess,L,A,PARAM),guess)
S=abs(S)


%Let's do some verif
%verif=SolvEco(S,L,A,PARAM);

%AMAZING

%Let's do a few plots
SS.k=S(1);    
SS.e=S(2);

%At SS
%SS.x=SS.e^(1/PARAM(10));
SS.x=SS.e/Er_4

SS.ECO=Economy(SS.k,N,SS.e,A,PARAM);

tau= 0.1 %individual tax rate
SS.T=SS.ECO(1)*A*sum(L,'all')*tau % Steady state Tax income

%We assume that agents are born with 0 asset
Ini.asset=0;

AO_l=AgentOptim(SS.ECO(2)*ones(20,1),SS.ECO(1),A,L(:,1),PARAM,Ini);
AO_h=AgentOptim(SS.ECO(2)*ones(20,1),SS.ECO(1),A,L(:,2),PARAM,Ini);
SS.asset_l=AO_l(:,2)
SS.asset_h=AO_h(:,2)

% Should SS.asset be individual for high and lowskilled workers? Probably
% yes

%sum(SS.asset)

figure(2)
plot(1:20,AO_l(:,1),"b--",1:20,AO_l(:,2),"r-")
figure(3)
plot(1:20,AO_h(:,1),"b--",1:20,AO_h(:,2),"r-")

AN=A*N
SS.y=F(SS.k,AN,SS.e,PARAM);
SS.y=SS.y(1);
%Verification that the income shares are right
A*N*SS.ECO(1)/SS.y;
SS.k*SS.ECO(2)/SS.y;
SS.e*SS.ECO(3)/SS.y;
%Niceeee


%% Let's try stuff








%% 4/ And now it's linked to the whole economy

%Let's model a shock
%PARAM(9)=2;
%Ressources became cheaper!

TIME=100


GUESS=[ones(TIME,1)*SS.k,ones(TIME,1)*SS.e, ones(TIME,1)*SS.T];

SolveTimeEco(GUESS,L,PARAM,SS);

options.MaxFunEvals = 30000 

Solutiont=abs(fsolve(@(x) SolveTimeEco(x,L,PARAM,SS),GUESS, options))


RPt=RessPress(Solutiont(:,2),PARAM,SS);

%Let's do some verif
verift=SolveTimeEco(Solutiont,L,PARAM,SS);

At=change(Solutiont(:,2),PARAM, Solutiont(:,3));
%At=2.5*ones(TIME,1);


ECOt=Economy_extended(Solutiont(:,1),L,Solutiont(:,2),At,PARAM);
AllAgents=AGENTS(ECOt(:,2),ECOt(:,1),At,L,PARAM,SS);


%Graph for agent's consumption choice

%low skilled
figure(3)
for jj=1:1:TIME
plot(1:TIME,AllAgents(:,jj,1),"b-",1:TIME,AllAgents(:,jj,2),"r:")
hold on
end

%High skilled workers
figure(4)
for jj=1:1:TIME
plot(1:TIME,AllAgents(:,jj,3),"b-",1:TIME,AllAgents(:,jj,4),"r:")
hold on
end




%ECOt+SS
WAGE=[SS.ECO(1)*ones(20,1);ECOt(:,1)];
IRATE=[SS.ECO(2)*ones(20,1);ECOt(:,2)];
GAMMA=[SS.ECO(3)*ones(20,1);ECOt(:,3)];
AA=[2.5*ones(20,1);At(:,1)];

%relevant values (marginal products)
figure(5)
plot(1:TIME+20,WAGE,"b-",1:TIME+20,IRATE,"r-",1:TIME+20,AA,"k--")


%relevant values (marginal products)
figure(5)
plot(1:TIME+20,WAGE,"b-",1:TIME+20,IRATE,"r-",1:TIME+20,AA,"k--")


Nt=15

Kt=Solutiont(:,1);
KK=[SS.k*ones(20,1);Kt];

Xt=RPt(:,1);
XX=[SS.x*ones(20,1);Xt];

Et=Solutiont(:,2);
EE=[SS.e*ones(20,1);Et];

Yt=F(Kt,Nt,Et,At,PARAM);
YY=[SS.y*ones(20,1);Yt];

%relevant values (K,E,X,Y)
figure(5)
plot(1:TIME+20,XX,"g--",1:TIME+20,KK,"r-",1:TIME+20,EE,"g-",1:TIME+20,YY,"k-")








Ct=nansum(AllAgents(:,:,1),2)


R=0.0098*ones(50,1)

change(R,PARAM)


