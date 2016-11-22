function y = dcm2qua(dcm)

angle = [atan2(-dcm(3,2),dcm(2,2)) atan2(-dcm(1,3),dcm(1,1)) asin(dcm(1,2))];

q1 = [cos(angle(2)/2) sin(angle(2)/2)*[0 1 0]]';
q2 = [cos(angle(3)/2) sin(angle(3)/2)*[0 0 1]]';
q3 = [cos(angle(1)/2) sin(angle(1)/2)*[1 0 0]]';
y = qmul(qmul(q1,q2),q3);