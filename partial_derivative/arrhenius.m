function f = arrhenius(xi, theta)
  f = exp(-theta(2)/xi) .* [1, -theta(1)/xi]';
end
