function [SolvEco]=SolvEco_skill(GUESS,L,N,A,PARAM,skill)

Kg=abs(GUESS(1));
Eg=abs(GUESS(2));

%N=sum(L,'all');

ECO=Economy_skill(Kg,N,Eg,A,PARAM,skill);

%We assume that agents are born with 0 asset
Ini.asset=0

AO=AgentOptim(ones(20,1).*ECO(2),ECO(1),A,L,PARAM,Ini);
K=sum(AO(:,2));
%For Kg we have an interest rate which in turn implies K (agent optim)
SolvEco(1)=K-Kg;

%At steady state we know that X=E^(1/Er_4)
%SolvEco(2)=ECO(3)-GAMMA(Eg^(1/PARAM(10)),Eg,PARAM);

%Other Version
%At steady state X=E/Er_4
SolvEco(2)=ECO(3)-GAMMA(Eg/PARAM(10),Eg,PARAM);
end




