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
[E_once,inliersE]=ransac5point_Num_of_Inliers_EveryE(pair.Fmatches(1:2,:),pair.Fmatches(3:4,:),tE,K,0);
%[E_once,inliersE]=ransac5point_Num_of_Inliers(pair.Fmatches(1:2,:),pair.Fmatches(3:4,:),tE,K,0);
RtE=EtoRtV3(E_once,K,K,pair.Fmatches(:,inliersE));
fprintf('E Eul:%.3f %.3f %.3f\n',dcm2eul(RtE(:,1:3))*180/pi);
%H = H_DLT( [pair.Fmatches(1:2,inliersF);pair.Fmatches(3:4,inliersF)] );
H = vgg_H_from_x_lin(pair.Fmatches(1:2,inliersF),pair.Fmatches(3:4,inliersF));
[ RtH,isPureRotation ]=HtoRtV5( H,K,K,pair.Fmatches(:,inliersF) );
fprintf('H Eul:%.3f %.3f %.3f\n',dcm2eul(RtH(:,1:3))*180/pi);
[ Right_Rt, Qt, D, M, P_E, P_H]=Universial_Relative_PoseV3_Out( K,K,pair.matches,0.1,1.5);
%%
%{
X1=vgg_get_homg(pair.Fmatches(1:2,inliersF));
X2=vgg_get_homg(pair.Fmatches(3:4,inliersF));
d = vgg_H_sampson_distance_sqr(H,X1,X2);
H_d=sum(d);
[ SampsonDist_row,SD_avg,sigma] = Fdistance( F,pair.Fmatches(:,inliersF) );
F_d=sum(SampsonDist_row);
fprintf('F Error:%.3f\n',F_d);
fprintf('H Error:%.3f\n',H_d);
disp('***');
%%

X_Triangualed_F = Triangulation(K*[eye(3) zeros(3,1)],K*RtF,pair.Fmatches(:,inliersF));
X_Triangualed_F_Inhomo = X_Triangualed_F(1:3,:)./X_Triangualed_F([4 4 4],:);
[ ~, ReprojectionPointInHomoF ] = MyReprojection( K*[eye(3) zeros(3,1)],K*RtF,X_Triangualed_F );
figure; 
plot3(X_Triangualed_F_Inhomo(1,:),X_Triangualed_F_Inhomo(2,:),X_Triangualed_F_Inhomo(3,:),'.g');
grid on;
title('F');
[ Error_F,ErrorSum_F ] = MyReErrorIndi( pair.Fmatches(:,inliersF),ReprojectionPointInHomoF );
%[ ErrorSum_F ] = MyReError( pair.Fmatches(:,inliersF),ReprojectionPointInHomoF );
X_Triangualed_H = Triangulation(K*[eye(3) zeros(3,1)],K*RtH,pair.Fmatches(:,inliersF));
X_Triangualed_H_Inhomo = X_Triangualed_H(1:3,:)./X_Triangualed_H([4 4 4],:);
[ ~, ReprojectionPointInHomoH ] = MyReprojection( K*[eye(3) zeros(3,1)],K*RtH,X_Triangualed_H );
figure; 
plot3(X_Triangualed_H_Inhomo(1,:),X_Triangualed_H_Inhomo(2,:),X_Triangualed_H_Inhomo(3,:),'.g');
grid on;
title('H');
[ Error_H,ErrorSum_H ] = MyReErrorIndi( pair.Fmatches(:,inliersF),ReprojectionPointInHomoH );
%[ ErrorSum_H ] = MyReError( pair.Fmatches(:,inliersF),ReprojectionPointInHomoH );
X_Triangualed_E = Triangulation(K*[eye(3) zeros(3,1)],K*RtE,pair.Fmatches(:,inliersE));
X_Triangualed_E_Inhomo = X_Triangualed_E(1:3,:)./X_Triangualed_E([4 4 4],:);
[ ~, ReprojectionPointInHomoE ] = MyReprojection( K*[eye(3) zeros(3,1)],K*RtE,X_Triangualed_E );
figure; 
plot3(X_Triangualed_E_Inhomo(1,:),X_Triangualed_E_Inhomo(2,:),X_Triangualed_E_Inhomo(3,:),'.g');
grid on;
title('E');
[ Error_E,ErrorSum_E ] = MyReErrorIndi( pair.Fmatches(:,inliersE),ReprojectionPointInHomoE );
figure
plot(Error_H,'.g');
title('contributing H');
figure
plot(Error_F,'.g');
title('contributing F');
figure
plot(Error_E,'.g');
title('contributing E');
fprintf('F Error:%.3f\n',ErrorSum_F);
fprintf('H Error:%.3f\n',ErrorSum_H);
fprintf('E Error:%.3f\n',ErrorSum_E);
%}
