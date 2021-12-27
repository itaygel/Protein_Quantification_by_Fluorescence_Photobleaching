function [b] = max_hist(v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


a=histogram(v,150);
y=a.Values;
x=a.BinEdges;
BinWidth=a.BinWidth;


[m,mm]=max(y);
b=x(mm)+BinWidth/2;









end

