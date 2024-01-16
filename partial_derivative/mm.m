function f = mm(xi, theta)
  f = [xi/(theta(2)+xi), -theta(1)*xi/(theta(2)+xi)^2]';
end
