function [ Translation_Scale,Tpoint ] = MyNormalizing( point )
%Normalization Operating for Coordinates of 2D image points
%Input:point 2*n matrix for n points
%Output:Translation_Scale,the normalizing matrix; Tpoint, normalized points  
%
%Reference Multiple View Geometry in Computer Vision, 2nd edition,P109
%Author: Peike Zhang, Yuanxin Wu, Qi Cai and Danping Zou
%2016
%Version 2
Xcentroid=0;
Ycentroid=0;
%%
for i=1:size(point,2)
    Xcentroid=Xcentroid+point(1,i);
    Ycentroid=Ycentroid+point(2,i);
end
Xcentroid=Xcentroid/size(point,2);
Ycentroid=Ycentroid/size(point,2);
%%
Length=0;
for j=1:size(point,2)
    Length=Length+sqrt(point(1,j)^2+point(2,j)^2);
    %Length=Length+norm(matches(1:2,j));
end
Length=Length/size(point,2);
Translation=[1,0,-Xcentroid;
             0,1,-Ycentroid;
             0,0,1];   
scale=sqrt(2)/Length;

Scale=[scale,0,0;
       0,scale,0;
       0,0    ,1];

%Scale=diag([scale,scale,1]);
Translation_Scale=Scale*Translation;
%%
    Newpoint=Translation_Scale*[point;ones(1,size(point,2))];
    Tpoint=Newpoint(1:2,:);
end 