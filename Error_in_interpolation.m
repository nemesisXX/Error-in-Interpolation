%%
% [-1,1] using chebyshev points
error_max = zeros(8,1);
error_2norm = zeros(8,1);
a = 1;
for n = [2,4,8,16,32,64,128,256]
    i = linspace(0,n,n+1);
    % Chebyshev points
    x = cos((i+0.5)*pi/(n+1));
    f = @(x) (x>=0); % Step function using logic
    y = f(x);
    x_0 = [-1:1/(10*n):1];
    y_0 = lagrange_interpolant(x,y,x_0);
    error_max(a) = max(abs(y_0-f(x_0)));
    sum_series = sum((y_0-f(x_0)).^2);
    error_2norm(a) = sqrt((2/(10*n))*sum_series);
    a = a+1;
end
error_max
error_2norm
function y0 = lagrange_interpolant(x,y,x0)
% x is the vector of abscissas
% y is the matching vector of ordinates
% x0 represents the target to be interpolated
% y0 represents the solution from the lagrange interpolation
y0 = 0;
n = length(x);
for j = 1:n
    t = 1;
    for i =1:n
        if i~=j
            t = t.*(x0-x(i))/(x(j)-x(i));
        end
    end
    y0 = y0+t*y(j);
end
end
