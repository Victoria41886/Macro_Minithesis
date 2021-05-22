function [SolveTimeEco]=SolveTimeEco(GUESS,L,PARAM,SS)
NN=size(GUESS,1);


Kg=abs(GUESS(:,1));
Eg=abs(GUESS(:,2));
Tg=abs(GUESS(:,3));

%Let's set A
Atime=change(Eg,PARAM,Tg);
%Atime=2.5*ones(NN,1);

N=sum(L, 'all');


%Simulating the economy
ECO=Economy_extended(Kg,L,Eg,Atime,PARAM);

AllAg=AGENTS(ECO(:,2),ECO(:,1),Atime,L,PARAM,SS);

Ktime=sum((nansum(AllAg(:,:,2),2)+nansum(AllAg(:,:,4),2)), 2);

Ttime=Atime(:,1).*ECO(:,1).*sum(L(:,1))+ Atime(:,2).*ECO(:,1).*sum(L(:,2));

RP=RessPress(Eg,PARAM,SS);

SolveTimeEco(:,1)=Ktime-Kg;
SolveTimeEco(:,2)=ECO(:,3)-RP(:,2);
SolveTimeEco(:,3)=Ttime-Tg;


sum(abs(SolveTimeEco),'all')
end


%Don't pay attention to what's below, I was just trying stuff


%Extracting our guess, K and E don't change in the last 40 periods
%Kg=[abs(GUESS(:,1));GUESS(NN,1)*ones(20,1)];
%Eg=[abs(GUESS(:,2));GUESS(NN,1)*ones(20,1)];



%SolveTimeEco(:,1)=[Ktime-Kg;Ktime(NN-19:NN)-Ktime(NN-19)];
%SolveTimeEco(:,2)=[ECO(:,3)-RP(:,2);ECO(NN-19:NN,3)-ECO(NN-19,3)];


%SolveTimeEco(NN:NN+20,:)=0;
%We don't need to check if the sum of Assets=Kapital in the last periods
%since we are not modelling all the agents


%We want to be at steady state in the last two "real" periods and the ones
%after
%SolveTimeEco(1:21,3)=Ktime(NN-20:NN)-Ktime(NN-21);
%SolveTimeEco(22:42,3)=Eg(NN-20:NN)-Eg(NN-21);
%SolveTimeEco(43:63,3)=Atime(NN-20:NN)-Atime(NN-21);





