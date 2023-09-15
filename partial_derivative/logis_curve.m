%% Extra: Functions that are for defining the partial derivatives of each parameters 
function f = logis_curve(xi, theta)
  deno = (1+theta(2)*exp(-theta(3)*xi));
  f = [1/deno, -(theta(1)*exp(-theta(3)*xi))/deno^2,  ...
    theta(1)*theta(2)*xi*exp(-theta(3)*xi)/deno^2]';
end