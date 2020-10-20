function [E1_mean,E2_mean,Symbol,O1_mean,O2_mean,Element_name]=cluster_profile(O1,E1,O1_idx,E1_idx,O2,E2,O2_idx,E2_idx,Symbol1,Symbol2,PeakName1,PeakName2,PeakName_intersect)
K1=max(max(O1_idx),max(E1_idx));
K2=max(max(O2_idx),max(E2_idx));
%%%Expression profile
Symbol=intersect(Symbol1,Symbol2);
[~,f1]=ismember(Symbol,Symbol1);
[~,f2]=ismember(Symbol,Symbol2);
for i=1:K1
    E1_mean(:,i)=mean(E1(f1,E1_idx==i),2);
end
for i=1:K2
    E2_mean(:,i)=mean(E2(f2,E2_idx==i),2);
end
%%%Accessibility profile
[~,f1]=ismember(PeakName_intersect(:,1),PeakName1);
[~,f2]=ismember(PeakName_intersect(:,2),PeakName2);
d1=ismember(PeakName1,PeakName_intersect(:,1));
d2=ismember(PeakName2,PeakName_intersect(:,2));
Element_name=[PeakName_intersect(:,1);PeakName1(d1==0);PeakName2(d2==0)];
m=size(PeakName_intersect,1);
O1_mean=zeros(length(Element_name),K1);
O2_mean=zeros(length(Element_name),K2);
for i=1:K1
    O1_mean(1:m,i)=mean(O1(f1,O1_idx==i)>0,2);
    O1_mean(1+m:m+sum(d1==0),i)=mean(O1(d1==0,O1_idx==i)>0,2);
end
for i=1:K2
    O2_mean(1:m,i)=mean(O2(f2,O2_idx==i)>0,2);
    O2_mean(m+sum(d1==0)+1:end,i)=mean(O2(d2==0,O2_idx==i)>0,2);
end
O1_mean=O1_mean./mean(O1_mean);
O2_mean=O2_mean./mean(O2_mean);