function o1=region_readCount(Peak_class,S,healthy_opn,PeakName,isAll)
K=max(S);
u=unique(cell2mat(Peak_class(:,2)));
PeakO=1*(healthy_opn>0);
if nargin == 4
for i=1:length(u)
    s=Peak_class(cell2mat(Peak_class(:,2))==u(i),1);
    [d f]=ismember(s,PeakName);
	o=sum(PeakO(f,:))./sum(PeakO);
    for j=1:K
        o1(i,j)=mean(o(S==j));
    end
end
else
    for i=1:length(u)
        s=Peak_class(cell2mat(Peak_class(:,2))==u(i),1);
        [d f]=ismember(s,PeakName);
        o=sum(PeakO(f,:))./sum(PeakO);
        o1(i,:)=o;
    end
end
o1=full(o1);