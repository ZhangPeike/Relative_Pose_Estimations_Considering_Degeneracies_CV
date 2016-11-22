function [misalign f fobj]= magimu_align_calibration_int4(wib, x, sample_t, C0)
% S1: Initial value of misalignment by analysis
% S2: Refined estimation by Newton method

deg = pi/180;
%% S1: Initial value of misalignment by analysis
len = length(x);
T = 5; % incremental time interval: 1s
dN = floor(T/sample_t);

A = zeros(3*floor(len/dN),9);
b = zeros(3*floor(len/dN),1);

global w_theld;
w_theld = 0.5;

ind = 1;
for k = 2:dN:len-dN
    if max(vmmag(wib(:,k-1:k+dN-1))) < w_theld % remove segement without significant angular velocity (< w_theld rad/s)
        continue;
    end

    tmp = zeros(3,9,dN+1);
    for j = k-1:k+dN-1 % 1 + (2:101), using last epoch k-1 as the start value
        m = x(:,j); w = wib(:,j);
        tmp(:,:,j-k+2) = kron(w',antisymm(m));
    end
    A(3*(ind-1)+1:3*ind,:) = sample_t/2*sum(tmp(:,:,1:end-1)+tmp(:,:,2:end),3);
    b(3*(ind-1)+1:3*ind,:) = x(:,k+dN-1)-x(:,k-1);

    ind = ind + 1;
end
c = A(1:3*(ind-1),:)\b(1:3*(ind-1));
C1 = reshape(c,3,3);
q = dcm2qua(dcmOrth(C1));
disp(['Initial Misalignment angle from IMU body to Mag:(deg) ' num2str(qua2eul(q)/deg)]);

%% S2: Refined solution by Newton method (MLE)
q0 = q;
% Added initial angle error to examine algorithm's robustness
angerr = 0*rand(3,1);
q0 = qmul(q0,eul2qua(angerr*deg));
disp(['Initial Misalignment angle from IMU body to Mag:(deg) ' num2str(qua2eul(q0)/deg) '  Added angle error (deg):' num2str(angerr')]);

lamda = 0;
f = zeros(3,1);
z = [q0; f; lamda];

clear obj;
disp('obj(z) = ');
for nIter = 1:5
    q = z(1:4); f = z(5:7); lamda = z(end);

    obj(nIter) = objective4(x,wib,q,f,T,sample_t,lamda);
    disp(['          ' num2str(obj(nIter)) ': ' num2str(z')]);

    J = zeros(8,1); H = zeros(8);
    Q0 = [1 0 0 0 1 0 0 0 1; 0 0 0 0 0 -1 0 1 0; 0 0 1 0 0 0 -1 0 0; 0 -1 0 1 0 0 0 0 0]';
    Q1 = [0 0 0 0 0 -1 0 1 0; 1 0 0 0 -1 0 0 0 -1; 0 1 0 1 0 0 0 0 0; 0 0 1 0 0 0 1 0 0]';
    Q2 = [0 0 1 0 0 0 -1 0 0; 0 1 0 1 0 0 0 0 0; -1 0 0 0 1 0 0 0 -1; 0 0 0 0 0 1 0 1 0]';
    Q3 = [0 -1 0 1 0 0 0 0 0; 0 0 1 0 0 0 1 0 0; 0 0 0 0 0 1 0 1 0; -1 0 0 0 -1 0 0 0 1]';

    for k = 2:dN:len-dN
        tmp = zeros(3,9,dN+1);
        tmp_f = zeros(3,3,dN+1);
        for j = k-1:k+dN-1 % 1 + (2:101), using last epoch k-1 as the start value
            m = x(:,j); w = wib(:,j);
            tmp(:,:,j-k+2) = kron(w',antisymm(m));
            tmp_f(:,:,j-k+2) = antisymm(m);
        end
        W = sample_t/2*sum(tmp(:,:,1:end-1)+tmp(:,:,2:end),3);
        M = sample_t/2*sum(tmp_f(:,:,1:end-1)+tmp_f(:,:,2:end),3);
        Cv = reshape(DCM(q),9,1);
        alfa = W'*(W*Cv - M*f - (x(:,k+dN-1)-x(:,k-1)));
        dCv = 2*[q(1)  q(2) -q(3) -q(4);...
            -q(4)  q(3)  q(2) -q(1);...
            q(3)  q(4)  q(1)  q(2);...
            q(4)  q(3)  q(2)  q(1);...
            q(1) -q(2)  q(3) -q(4);...
            -q(2) -q(1)  q(4)  q(3);...
            -q(3)  q(4) -q(1)  q(2);...
            q(2)  q(1)  q(4)  q(3);...
            q(1) -q(2) -q(3)  q(4)];

        J(1:4) = J(1:4) + 2*dCv'*alfa;  % Jq
        J(5:7) = J(5:7) - 2*M'*(W*Cv - M*f - (x(:,k+dN-1)-x(:,k-1))); %Jf
        H(1:4,1:4) = H(1:4,1:4) + 2*(dCv'*W'*W*dCv + [alfa'*Q0; alfa'*Q1; alfa'*Q2; alfa'*Q3]); % Jqq
        H(1:4,5:7) = H(1:4,5:7) - 2*dCv'*W'*M; % Jqf
        H(5:7,5:7) = H(5:7,5:7) + 2*M'*M; % Jff
    end
    H(1:4, end) = 2*q; % Jq_lamda
    H(5:7,1:4) = H(1:4,5:7)'; %Jfq;
    H(end,1:4) = H(1:4, end)';
    J(1:4) = 2*lamda*q + J(1:4); % Jq
    J(end) = q'*q-1; % J_lamda
    H(1:4,1:4) = 2*lamda*eye(4) + H(1:4,1:4); % Jqq

    dz = -H\J;%-inv(H)*J;
    z = z + dz;
end

q = z(1:4); f = z(5:7); lamda = z(end);
obj(end+1) = objective4(x,wib,q,f,T,sample_t,lamda);
disp(['          ' num2str(obj(nIter)) ': ' num2str(z')]);

disp(['MLE Misalignment angle from IMU body to Mag:(deg) ' num2str(qua2eul(q)/deg)]);
disp(['MLE Gyro bias:(deg/s) ' num2str((DCM(q)'*f)'/deg)]);

misalign = q;

%% for debug
ind = 1;
alfa = zeros(3,floor(len/dN));
beta = zeros(3,floor(len/dN));

for k = 2:dN:len-dN
    tmp = zeros(3,9,dN+1);
    tmp_f = zeros(3,3,dN+1);
    for j = k-1:k+dN-1 % 1 + (2:101), using last epoch k-1 as the start value
        m = x(:,j); w = wib(:,j);
        tmp(:,:,j-k+2) = kron(w',antisymm(m));
        tmp_f(:,:,j-k+2) = antisymm(m);
    end
    W = sample_t/2*sum(tmp(:,:,1:end-1)+tmp(:,:,2:end),3);
    M = sample_t/2*sum(tmp_f(:,:,1:end-1)+tmp_f(:,:,2:end),3);
    alfa(:,ind) = W*reshape(DCM(q),9,1) - M*f;
    beta(:,ind) = (x(:,k+dN-1)-x(:,k-1));
    ind = ind + 1;
end

fobj = [alfa; beta];

