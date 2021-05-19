function [SolveTimeEco_result]=SolveTimeEco_skill(GUESS,L,N,PARAM,SS,skill)

%Extracting our guess
Kg=abs(GUESS(:,1));
Eg=abs(GUESS(:,2));

% A modelled as a result of emissions, like in the paper
Atime=change_skill(Eg,PARAM,skill)

%N=sum(L);

%Simulating the economy
ECO=Economy_skill(Kg,N,Eg,Atime,PARAM,skill);


%The issue now is that we need model the new AND the old agents in the
%economy

AllAg=AGENTS(ECO(:,2),ECO(:,1),Atime,L,PARAM,SS);

Ktime=nansum(AllAg(:,:,2),1).';

RP=RessPress(Eg,PARAM,SS)

SolveTimeEco_result(:,1)=Ktime-Kg;
%SolveTimeEco(:,2)=ECO(:,3)-RP(:,2);
SolveTimeEco_result(:,2)=ECO(:,3)- RP(:,2) ;            % I changed the name of the output because this gave me also error messages
%sum(abs(SolveTEco),'all')                       % this command gave me error messages 
end






