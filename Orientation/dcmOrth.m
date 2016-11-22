function y = dcmOrth(x)

while 1
    x = (inv(x')+x)/2;
    if sum(abs(diag(x*x'-eye(3))))<1e-10
        break;
    end
end
y = x;