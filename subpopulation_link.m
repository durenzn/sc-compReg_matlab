function Match=subpopulation_link(E_mean_healthy,E_mean_CLL,O_mean_healthy,O_mean_CLL)
r1=corr(E_mean_healthy,E_mean_CLL);
r2=corr(O_mean_healthy,O_mean_CLL);
rr1=r1-sum(r1).*sum(r1')'/sum(sum(r1));
rr2=r2-sum(r2).*sum(r2')'/sum(sum(r2));
rr=rr1+rr2;
clear a
[a(:,1),a(:,2)]=find(rr>0);
b=find(rr>0);
a(:,3)=r1(b);
a(:,4)=r2(b);
a(:,5)=rr1(b);
a(:,6)=rr2(b);
a(:,7)=rr(b);
[d f]=sort(a(:,7),'descend');
a=a(f,:);
a=a((a(:,5)>0).*(a(:,6)>0)==1,:);

S1=[];
S2=[];
Match=[];
for i=1:size(a,1)
    idx=ismember(a(i,1),S1)+ismember(a(i,2),S2);
    S1=[S1;a(i,1)];
    S2=[S2;a(i,2)];
    if idx<2
        Match=[Match;a(i,:)];
    end
end