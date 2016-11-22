function product = qmul(u1,u2)
% multiplication of quaternions represented by 1*4 vector

% forcing to be row vector
if size(u1,1) ~= 1
    u1 = u1';
end
if size(u2,1) ~= 1
    u2 = u2';
end

% output column vector
product = [u1(1)*u2(1)-dot(u1(2:4),u2(2:4)) u1(1)*u2(2:4)+u2(1)*u1(2:4)+cross(u1(2:4),u2(2:4))]';