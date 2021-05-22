function [AgentOptim]=AgentOptim(R,W,A,L,PARAM,SS)

tau=0.1
S=(1-tau)*W.*A.*L;

AsIni=SS.asset(21-size(R,1));

CR=cumprod((1+R-PARAM(5)))/(1+R(1)-PARAM(5));
DF=1./CR;

Cinitial=(sum(S.*DF)+AsIni*(1+R(1)-PARAM(5)))/sum(DF.^(1-1/PARAM(6)).*PARAM(4).^([0:size(R,1)-1].'/PARAM(6)));
%The AsIni*(1+R(1)-PARAM(5)) is to allow us to use the function for agent
%with non-zero intial asset i.e. the agents born before period 1 but still
%alive


AgentOptim(:,1)=ConsOpt(Cinitial,R,PARAM);
AgentOptim(:,2)=ASSET(S,AgentOptim(:,1),R,PARAM,SS);

%Give a colum with consumption and a column with assets based on the agent
%maximization prog
end
