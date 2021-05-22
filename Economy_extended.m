function [Economy]=Economy_extended(K,L,E,A,PARAM)
AN= sum(L(:,1)).*A(:,1)+sum(L(:,2)).*A(:,2);

Economy(:,1)=((1-PARAM(1)-PARAM(2))*F(K,AN,E,PARAM)./(AN));
Economy(:,2)=(PARAM(1)*F(K,AN,E,PARAM)./K);
Economy(:,3)=(PARAM(2)*F(K,AN,E,PARAM)./E);
end