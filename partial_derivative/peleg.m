function f = peleg(xi, theta)
  deno = (theta(1) + xi * theta(2))^2;
  f = [-xi / deno; -xi^2 / deno];
end
