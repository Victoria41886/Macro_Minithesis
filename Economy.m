function [Economy]=Economy(K,N,E,A,PARAM)
AN=N.*A
Economy(:,1)=((1-PARAM(1)-PARAM(2))*F(K,AN,E,PARAM)./(AN));
Economy(:,2)=(PARAM(1)*F(K,AN,E,PARAM)./K);
Economy(:,3)=(PARAM(2)*F(K,AN,E,PARAM)./E);


end

