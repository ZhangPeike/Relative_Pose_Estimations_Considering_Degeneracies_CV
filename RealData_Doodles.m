clear;
close all;
warning off all;
disp('For Planar Structure,H F E');
load('DoodleMatches.mat');
tF=0.001;
tE=1/(3520*3520);
%tE=0.00000001;
K =[3520,0,1760;
    0,3520,1760;
    0,0,1];
frames.imsize=[4160,3120];
index=1;
%{
frames.images{1}='G:\ComputerVision2016\3D Program\摄像机标定\棋盘格小黄鸭手机低分辨率\20160817标定（涂鸦和展板）\IMG_20160812_161145_HDR.jpg';
frames.images{2}='G:\ComputerVision2016\3D Program\摄像机标定\棋盘格小黄鸭手机低分辨率\20160817标定（涂鸦和展板）\IMG_20160812_161219_HDR.jpg';
frames.imsize = size(imresize(imread(frames.images{1}),1/7));
K=[500,0,250;0,500,250;0,0,1];
pair=match2viewSIFT(frames,1,2);
%}
%showMatchesTopLeft(pair,frames);
pair.Fmatches=[];
pair.Fmatches(:,1)=pair.matches(:,1);
for j=2:size(pair.matches,2)
    if pair.Fmatches(1,index)==pair.matches(1,j)&&pair.Fmatches(2,index)==pair.matches(2,j)||pair.Fmatches(3,index)==pair.matches(3,j)&&pair.Fmatches(4,index)==pair.matches(4,j)
        fprintf('第%d个Removed! \n',j);  
    else
        index=index+1;
        pair.Fmatches(:,index)=pair.matches(:,j);
    end
end
[F,inliersF]=ransacfitfundmatrix(pair.Fmatches(1:2,:),pair.Fmatches(3:4,:),tF,0);
%{
inliersF(:,46)=[];
inliersF(:,137)=[];
inliersF(:,148)=[];
%}
E=K'*F*K;
RtF=EtoRtV3(E,K,K,pair.Fmatches(:,inliersF));
fprintf('F Eul:%.3f %.3f %.3f\n',dcm2eul(RtF(:,1:3))*180/pi);
%H = H_DLT( [pair.Fmatches(1:2,inliersF);pair.Fmatches(3:4,inliersF)] );
H = vgg_H_from_x_lin(pair.Fmatches(1:2,inliersF),pair.Fmatches(3:4,inliersF));
[ RtH,isPureRotation ]=HtoRtV5( H,K,K,pair.Fmatches(:,inliersF) );
fprintf('H Eul:%.3f %.3f %.3f\n',dcm2eul(RtH(:,1:3))*180/pi);
[ Right_Rt, Qt, D, M, P_E, P_H]=Universial_Relative_PoseV3_Out( K,K,pair.matches,0.1,1.5);
