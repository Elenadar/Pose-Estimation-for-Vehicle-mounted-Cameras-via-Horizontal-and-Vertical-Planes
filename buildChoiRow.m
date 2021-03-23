function [A_r,B_r] = buildChoiRow(p,q)
    x = p(1);
    y = p(2);
    z = p(3);
    x_ = q(1);
    y_ = q(2);
    z_ = q(3);
    A_r = [x*y_ -z*y_];
    B_r = [x_*y z_*y];
end