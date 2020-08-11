function[grads] = gradient_compute(x_tilde,knots,order)
L = size(x_tilde,2);
for i = 1:L
    grads{i} = gradient_bspline(x_tilde(1:2,i),knots,order);
end