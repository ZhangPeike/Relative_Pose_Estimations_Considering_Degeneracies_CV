% RANSACFITHOMOGRAPHY - fits 2D homography using RANSAC
%
% Usage:   [H, inliers] = ransacfithomography_vgg(x1, x2, t)
%
% Arguments:
%          x1  - 2xN or 3xN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that point coordinates are not normalised to that their
%                mean distance from the origin is sqrt(2) any more.
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          H       - The 3x3 homography such that x2 = H*x1.
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.
%
% See Also: ransac, homography2d, homography1d

% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% February 2004 - original version
% July     2004 - error in denormalising corrected (thanks to Andrew Stein)

% Adapted to use vgg functions by Peter Kovesi and Andrew Zisserman
% 2012 Updated by Alexander Khanin
% 2017-01-21 10:25:44 upgraded by Peike Zhang
% 1.Estimating model with normalization but reshape the model into original
% one corresponding to raw data input, then to judge inliers and outliers.
% 2.Algebraic Distance is replaced by sampson distance, which is
% approximate to geometric distance and easy.
function [H, inliers] = ransacfithomography_Est_with_Normalizing_Err_original(x1, x2, t)

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [rows,npts] = size(x1);
    if rows == 2
        %error('x1 and x2 must have 3 rows');
        %disp('Input for estimating H is 2 rows matches');
        [~,npts] = size(x1);
        x1=[x1;ones(1,npts)];
        x2=[x2;ones(1,npts)];
    end
    
    if npts < 4
        error('Must have at least 4 points to fit homography');
    end
    %global x;
    %x=[x1;x2];
    % Normalise each set of points so that the origin is at centroid and
    % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    % scale parameter is 1.  Note that 'homography2d' will also call
    % 'normalise2dpts' but the code in 'ransac' that calls the distance
    % function will not - so it is best that we normalise beforehand.
    %[x1, T1] = normalise2dpts(x1);
    %[x2, T2] = normalise2dpts(x2);
    
    s = 4;  % Minimum No of points needed to fit a homography.
    
    fittingfn = @wrap_vgg_homography2d;
    distfn    = @homogdist2d;
    degenfn   = @isdegenerate;
    
    % x1 and x2 are 'stacked' to create a 6xN array for ransac
    [~, inliers] = ransac([x1; x2], fittingfn, distfn, degenfn, s, t);

    % Now do a final least squares fit on the data points considered to
    % be inliers.
    H = vgg_H_from_x_lin(x1(:,inliers), x2(:,inliers));
    %H = vgg_H_from_x_nonlin(Hlin,x1(:,inliers), x2(:,inliers));
    % Denormalise
    %H = T2\H*T1;    

%----------------------------------------------------------------------
% Function to evaluate the symmetric transfer error of a homography with
% respect to a set of matched points as needed by RANSAC.
function [inliers, H] = homogdist2d(H, x, t)
    %Distance is squared!
    SampsonDist_row=vgg_H_sampson_distance_sqr(H,x(1:3,:),x(4:6,:));
    inliers=find(SampsonDist_row<t^2);   
end
%{
function [inliers, H] = homogdist2d(H, x, t)
    
    x1 = x(1:3,:);   % Extract x1 and x2 from x
    x2 = x(4:6,:);    
    
    % Calculate, in both directions, the transfered points    
    Hx1    = H*x1;
    invHx2 = H\x2;
    x1=x1./x1([3 3 3],:);
    x2=x2./x2([3 3 3],:);
    Hx1=Hx1./Hx1([3 3 3],:);
    invHx2=invHx2./invHx2([3 3 3],:);
    d2 = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
    inliers = find(abs(d2) < t); 
end
%}

%----------------------------------------------------------------------
% Function to determine if a set of 4 pairs of matched  points give rise
% to a degeneracy in the calculation of a homography as needed by RANSAC.
% This involves testing whether any 3 of the 4 points in each set is
% colinear. 
%{    
function r = isdegenerate(x)

    x1 = x(1:3,:);    % Extract x1 and x2 from x
    x2 = x(4:6,:);    
    r = ...
    iscolinear(x1(:,1),x1(:,2),x1(:,3)) | ...
    iscolinear(x1(:,1),x1(:,2),x1(:,4)) | ...
    iscolinear(x1(:,1),x1(:,3),x1(:,4)) | ...
    iscolinear(x1(:,2),x1(:,3),x1(:,4)) | ...
    iscolinear(x2(:,1),x2(:,2),x2(:,3)) | ...
    iscolinear(x2(:,1),x2(:,2),x2(:,4)) | ...
    iscolinear(x2(:,1),x2(:,3),x2(:,4)) | ...
    iscolinear(x2(:,2),x2(:,3),x2(:,4));
%}
function r = isdegenerate(x)
    r = 0;    
end
function H = wrap_vgg_homography2d(x)
    H = vgg_H_from_x_lin(x(1:3,:),x(4:6,:));
end
end
