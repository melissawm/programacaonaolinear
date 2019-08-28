%%%%%%%%%%%%%%%%%%%%%
%% Método de Newton 
%%%%%%%%%%%%%%%%%%%%%

% A função a ser minimizada deve estar num arquivo chamado
%                        fun.m
% e o diretório pode ser alterado usando-se o comando abaixo:

addpath("problemas/rosenbrock");

% Por exemplo, caso o problema estivesse em outro diretório, usaríamos
% addtopath("problemas/outro_problema")

% Não se esqueça de calcular o gradiente e a 
% Hessiana (necessária para o método de Newton)

% Escolha um ponto inicial:
x = [-1,2]';

% Escolha uma tolerância a ser aplicada ao gradiente para 
% detectar a convergência:
eps = 1e-5;

% Se quiser ver um gráfico do caminho das iteradas ao final da execução do método, selecione 
%                grafico = 1

grafico = 1;

%% Método de Newton

itnum = 0;

% Vamos guardar aqui as coordenadas dos pontos para fazer
% o gráfico do caminho das iteradas
if grafico
    % Vamos fazer um gráfico com as curvas
    % de nível da função e verificar o 
    % caminho das iteradas
    X = linspace(-2, 2);
    Y = linspace(-5, 5);
    [XX, YY] = meshgrid(X, Y);
    %Z = fun([XX, YY]);
    Z = (1-XX).^2+10*(YY-XX.^2).^2;
    levels = -1:2:80;
    contour(XX, YY, Z, levels, 'linewidth', 2);
    hold on
end

% Calculamos o gradiente e começamos o laço principal:
g = grad(x);

lista_x = x;

while (norm(g) > eps) && (itnum < 2000)
    % Direção de Newton:
    H = hessian(x);
    x = x - H\g;
    
    % Acrescentamos x na posição itnum+1 na lista_x pois itnum inicialmente
    % é zero, e os índices de vetores no matlab começam em 1.
    lista_x(:, end+1) = x;

    disp(['Iteração ', num2str(itnum)])
    x
    
    itnum = itnum + 1;
    g = grad(x);
end

if grafico
    plot(lista_x(1,:), lista_x(2,:), 'mo-', 'linewidth', 2);
end

% Vamos verificar se obtivemos convergência, ou se 
% o método parou por ter atingido o número máximo de iterações:
if norm(g) <= eps
    disp("========================================")
    disp("           Convergência!                ")
    disp("  x final: ")
    x
    disp("  f final: ")
    fun(x)
    disp("  g final: ")
    g
    disp("========================================")
else
    disp("++++++++++++++++++++++++++++++++++++++++")
    disp("           Não convergimos!             ")
    disp("++++++++++++++++++++++++++++++++++++++++")
end
hold off;
