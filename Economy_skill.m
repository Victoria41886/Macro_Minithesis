function [Economy]=Economy_skill(K,N,E,A,PARAM,skill)

Economy(:,1)=skill*((1-PARAM(1)-PARAM(2))*F(K,N,E,A,PARAM)./(A.*N));
Economy(:,2)=(PARAM(1)*F(K,N,E,A,PARAM)./K);
Economy(:,3)=(PARAM(2)*F(K,N,E,A,PARAM)./E);

end

