function [A Feature]=getmatrix(filename,feature,fileHead)
D=dlmread(filename,' ',fileHead,0);
D_num=dlmread(filename,' ',[fileHead-1,0,fileHead-1,1]);
A=sparse(D(:,1),D(:,2),D(:,3),D_num(1),D_num(2));
fileID = fopen(feature);
C = textscan(fileID,'%s %s %s %*[^\n]');
fclose(fileID);
Feature=[C{1,1} C{1,2} C{1,3}];
