function [g] = grad(x)
    g = zeros(2,1);
    g(1) = 2*x(1)-2-40*x(1)*(x(2)-x(1).^2);
    g(2) = 20*(x(2)-x(1).^2);
end