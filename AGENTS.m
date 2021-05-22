function [AGENTS]=AGENTS(Rtime,Wtime,Atime,L,PARAM,SS)

NN=100;

Cagents=NaN(NN,NN);
Aagents=NaN(NN,NN);


Atime_l=Atime(:,1)
Atime_h=Atime(:,2)

%For the old agents
for jj=1:19
    R=Rtime(1:jj);
    W=Wtime(1:jj);
    A=Atime_l(1:jj);
    Lold=L(21-jj:20,1);
    Asset_ss=SS.asset_l;
    Agent=AgentOptim2(R,W,A,Lold,PARAM,Asset_ss);
    Cagents_l(1:jj,jj)=Agent(:,1);
    Aagents_l(1:jj,jj)=Agent(:,2);
end

%For the new agents
for ii=1:NN-19
    R=Rtime(ii:ii+19);
    W=Wtime(ii:ii+19);
    A=Atime_l(ii:ii+19);
    Asset_ss=SS.asset_l;
    Agent=AgentOptim2(R,W,A,L(:,1),PARAM,Asset_ss);
    Cagents_l(ii:ii+19,ii+19)=Agent(:,1);
    Aagents_l(ii:ii+19,ii+19)=Agent(:,2);
end





% high skilled workers
%old agents
for jj=1:19
    R=Rtime(1:jj);
    W=Wtime(1:jj);
    A=Atime_h(1:jj);
    Lold=L(21-jj:20,2);
    Asset_ss=SS.asset_h;
    Agent=AgentOptim2(R,W,A,Lold,PARAM,Asset_ss);
    Cagents_h(1:jj,jj)=Agent(:,1);
    Aagents_h(1:jj,jj)=Agent(:,2);
end


%For the new agents
for ii=1:NN-19
    R=Rtime(ii:ii+19);
    W=Wtime(ii:ii+19);
    A=Atime_h(ii:ii+19);
    Asset_ss=SS.asset_h;
    Agent=AgentOptim2(R,W,A,L(:,2),PARAM,Asset_ss);
    Cagents_h(ii:ii+19,ii+19)=Agent(:,1);
    Aagents_h(ii:ii+19,ii+19)=Agent(:,2);
end

%results low skilled workers
AGENTS(:,:,1)=Cagents_l;
AGENTS(:,:,2)=Aagents_l;

%results high skilled workers
AGENTS(:,:,3)=Cagents_h
AGENTS(:,:,4)=Aagents_h

end


