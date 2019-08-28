function [f] = fun(x)
    
    % x deve ser um vetor com 2 componentes
    
    f = (1-x(1)).^2+10*(x(2)-x(1).^2).^2;
end