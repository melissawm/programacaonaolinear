function [H] = hessian(x)
    H = zeros(2,2);
    H(1,1) = 2-40*(x(2)-x(1).^2)+80*x(1).^2;
    H(1,2) = -40*x(1);
    H(2,1) = H(1,2);
    H(2,2) = 20;
end