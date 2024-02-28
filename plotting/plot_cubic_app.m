%%Compute D-optimal designs for polynomial regression with exact n points
 
clear;
criterion = "A";
tol = 1E-4; % for finding and filtering out the points 
N = [5, 11,21,51,101,201, 501 1001,2001]';
% N = 5;
loss = zeros(2, size(N,1));
% runningtime = cputime;  %record computation time

 % number of design points for initial design
  
a =  -1;   %[a, b] is the design space
b =   1;  
p = 3;            % degree of polynomial regression model  
q = p+1; % how many beta's (degree + 1 intercept term)



% The following vectors and matrices are used in the information
% matrices below.
       

for i = 1:size(N,1)
  u = linspace(a, b, N(i)); %equally spaced N points in [a,b]
  f = power(u(:), 0:p);   
  %% 1. Compute the initial proxy approximate design 
  cvx_begin
    cvx_precision high
    variable w(1,N(i));
    expression A(q,q); 
    
    % here we compute the information matrix at
    for j=1:N(i)
      f1 = f(j,:)';      
      A = A + (f1 * f1') * w(j);
    end
      
    if criterion == "D"
      minimize (-log(det_rootn(A)));
    elseif criterion == "A"
      minimize( trace_inv(A) );   %A-opt
    else
      fprintf('Does not run.');
    end
    0 <= w <= 1;
    sum(w)==1;
  cvx_end

  % organizing the result
  design_app = [u(find(w(1,:)>tol)); w(find(w(1,:)>tol))]; %optimal design 
  d00 = design_app(1,:); % support points
  w00 = design_app(2,:); % optimal weight
  L00 = cvx_optval; % optimal objective value
  loss(1, i) = L00;
  loss(2, i) = N(i);
end


%% the optimal values
design_true_D = [-1, -0.447, 0.447, 1; 0.25, 0.25, 0.25,0.25];
design_true_A = [-1, -0.464, 0.464, 1; 0.151, 0.349, 0.349, 0.151];
B = zeros(q,q);
for j=1:size(design_true_D,2)
  x = design_true_D(1, j);
  f = power(x, 0:p)';
  w = design_true_D(2, j);
  B = B + (f * f') * w;
end
opt_D = -log(det(B)^(1/q));
C = zeros(q,q);
for j=1:size(design_true_A,2)
  x = design_true_A(1, j);
  f = power(x, 0:p)';
  w = design_true_A(2, j);
  C = C + (f * f') * w;
end

opt_A = trace(inv(C));   %A-opt
