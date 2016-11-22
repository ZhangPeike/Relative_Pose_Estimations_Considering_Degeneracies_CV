function y = dcm2eul(dcm)

y = [atan2(-dcm(3,2),dcm(2,2)) atan2(-dcm(1,3),dcm(1,1)) asin(dcm(1,2))];