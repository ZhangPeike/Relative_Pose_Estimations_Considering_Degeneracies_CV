function [ T,Tpoint ] = NormalizingV2( point )
%对图像点的归一化处理
%Normalization
%point是非齐次坐标，Newpoint临时是齐次坐标，后转为非齐次
%过程，首先使图像点的中心平移到坐标原点，
%平移和缩放用3*3矩阵表示
%参考MVG P109
%先求中心，再求平均距离
%T=[*,*,*;*,*,*;0,0,1]
%作者：张培科
%2016年2月6日10:32:25 
%Version 2
%Newpoint，是其次坐标（*；*；1），再取非齐次形式
%输入输出的是非齐次坐标
%%
%点云中心
centroid=mean(point);
%%
%将点云中心移动到远点
newpoint(1,:)=point(1,:)-centroid(1);
newpoint(2,:)=point(2,:)-centroid(2);
dist=sqrt(newpoint(1,:).^2+newpoint(2,:).^2);
meandist = mean(dist(:));  % Ensure dist is a column vector for Octave 3.0.1
scale = sqrt(2)/meandist;
T = [scale   0   -scale*centroid(1)
     0     scale -scale*centroid(2)
     0       0      1      ];
Newpoint=T*[point;ones(1,size(point,2))];
%归一化后的点的非齐次坐标
Tpoint=Newpoint(1:2,:);
end 