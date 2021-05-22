function [SolvEco]=SolvEco(GUESS,L,A,PARAM)

Kg=abs(GUESS(1));
Eg=abs(GUESS(2));

N=sum(L,'all');

ECO=Economy(Kg,N,Eg,A,PARAM);

%We assume that agents are born with 0 asset
Ini.asset=0;
tau=0.1
AO_l=AgentOptim(ones(20,1).*ECO(2),ECO(1),A,L(:,1),PARAM,Ini);
AO_h=AgentOptim(ones(20,1).*ECO(2),ECO(1),A,L(:,2),PARAM,Ini);
K=sum(AO_l(:,2))+sum(AO_h(:,2));
%For Kg we have an interest rate which in turn implies K (agent optim)
SolvEco(1)=K-Kg;


%Other Version
%At steady state X=E/Er_4
SolvEco(2)=ECO(3)-GAMMA(Eg/PARAM(10),Eg,PARAM);

end




