%%Compute D-optimal designs for polynomial regression with exact n points
 
clear;
criterion = "A";
tol = 1E-4; % for finding and filtering out the points 
% N = 21;
Nsim = 500;
% n = [10, 20 , 30, 40, 50]';
n = [10, 20,30,40]';
tol_annealing = 1E-40;

% number of design points for initial design
  
a =  -1;   %[a, b] is the design space
x0 = 0;
b =   1;  
p = 3;            % degree of polynomial regression model  
q = p+1; % how many beta's (degree + 1 intercept term)



% The following vectors and matrices are used in the information
% matrices below.


if criterion == "D"
    design_true = [-1, -0.447, 0.447, 1; 0.25, 0.25, 0.25,0.25];
elseif criterion == "A"
  design_true = [-1, -0.464, 0.464, 1; 0.151, 0.349, 0.349, 0.151];
else
  fprintf('Does not run.');
end


w00 = design_true(2,:);
d00 = design_true(1,:);
n_pt = length(d00);
n_half = n_pt/2;
d00_half = d00((n_half+1):end);
w00_half = w00((n_half+1):end);
design_true_half = design_true(:,(n_half+1):end);


c0 = 1; % max number of points to be changed in the annealing algorithm
Nt = 200; % number of iterations per temperature change
T0 = 5; % initial temperature
M0 = 2000; % number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0 * alpha^M0; %minimum temperature

delta = 2*(b-a)/500; % neighbourhood size, in this setting, it is 0.2
w00 = design_true(2,:);
d00 = design_true(1,:);

my_loss = zeros(5, size(n,1));
for pig = 1:size(n,1)
  disp(n(pig))
  n_kai = n(pig)/2;
  LOSS = zeros(Nsim, 4);
  for ell=1:Nsim 
    disp(ell)
    rng(ell);  %random seed number
    % rng(466)

    w01_half = initializeExact(w00_half, n_kai)/2'; %convert approximate design lazily to an exact design
    d0_half = design_true_half(1, :);
    w0_half = w01_half;
    index_zero = (d0_half == 0)';

    w0_half(index_zero) = w0_half(index_zero)/2;
    FIM = zeros(q, q);
    d0_full = sort([-d0_half, d0_half]);
    w0_full = [flip(w0_half), w0_half];
    
    k = length(w0_half); 
    
    for j=1:size(design_true,2)
        fx = power(d0_full(j), 0:p)';
        FIM = FIM + (fx * fx') * w0_full(j);
   end

  L0 = -log(det(FIM)^(1/q));
    if criterion == "D"
      L0 = -log(det(FIM)^(1/q));
    elseif criterion == "A"
      L0 = trace(inv(FIM));   %A-opt
    else
      fprintf('Does not run.');
    end
  
  
  % store loss at each iteration for plotting
  
  loss = zeros(M0*Nt,1);
  loss(1) = L0;

  %% 3. ANNEALING ALGORITTHM
  
  % This algorithm is an implementation of simulated annealing with a constraint
  % on the generation of points. New generated points for candidate designs
  % are within a predefinied neighbourhood of some random subset of points from the
  % previous design. This is to test if an exact design can be further
  % optimized by replacing some of the support points with nearby points.



  num_iters = 1;
  T = T0;
  L_prev = 0;
  while(T > Tmin && abs(L_prev - L0) > tol_annealing )
  % while(T > Tmin)
    L_prev = L0;
    % GENERATING RANDOM CANDIDATE DESIGN
    for h = 1:Nt
      k = length(w0_half); % current number of support points
      ci =  min([k, randi([1 c0])]); % number of points to be replaced
      toRemove = randperm(k, ci); % randomly select ci points to remove, return an index
      
      % creating new candidate design
    wi_half = zeros(k+ci, 1)'; 
    di_half = zeros(k+ci, 1)';
    
    
    % copy previous design
    wi_half(1:k) = w0_half;
    di_half(1:k) = d0_half;
      
      % remove selected points, here, there may be multiple points
      % $ here we essentially decrease the weight by 1/n (remove one point)
       for j = toRemove
          wi_half(j) = wi_half(j) - 1/n(pig); 
      end
    
      % add weight for new points
      % $ here we add the new (equal) weight to the new points
      % $ as we know that in D-opt, the weight are equal
    for j = 1:ci
      wi_half(k+j) = wi_half(k+j) + 1/n(pig); 
    end
      
      % generate new points in neighbourhood of those previously removed
      i = 1;
    for j = toRemove
      % $ delta is the neighbour
      di_nb = di_half(j) + unifrnd(-1,1)*delta; % ! check this neighbour definition
      
      % di_nb = di(j) + (2*rand(1)-1)*delta; % ! check this neighbour definition
      di_nb = min([max([x0, di_nb]),b]); % this needs more work for high-D 
      % ! actually here we can use the rejection proposal kind of approach
      di_half(k+i) = di_nb;
      i = i + 1; 
    end
      % points are generated in a circular, spherical, or hyperspherical
      % region (depending on the dimension of the design) of radius
      % delta
      
      % remove support points with zero weight
     di_half = di_half(wi_half>tol);
    wi_half = wi_half(wi_half>tol);
      
      index_zero = (di_half == 0)';
  
      wi_half(index_zero) = wi_half(index_zero)/2;
      
      di_full = sort([-di_half, di_half]);
      wi_full = [flip(wi_half), wi_half];
      % compute loss of candidate design
       FIMi = zeros(q, q);
    for j=1:length(wi_full)   
      fx = power(di_full(j), 0:p)';
      FIMi = FIMi + (fx * fx') * wi_full(j);
    end
    
   
       if criterion == "D"
         Li = -log(det(FIMi)^(1/q));
      elseif criterion == "A"
        Li = trace(inv(FIMi));   %A-opt
      else
        fprintf('Does not run.');
      end   
      % PROCEED WITH ANNEALING STEP
      prob = exp(-(Li-L0)/T); % acceptance probability
      if prob > rand(1) % criterion for accepting random design
        L0 = Li;
        d0_half = di_half;
        w0_half = wi_half; 
        num_iters = num_iters + 1;
        loss(num_iters) = L0;
      end
    end
    T = alpha*T;
  end
  loss1 = loss(2:num_iters);
  di_final = sort([-di_half, di_half]);
  wi_final = [flip(wi_half), wi_half];
    %% 4. PLOTTING RESULTS
    % here, we group the values that are the same together
  design_ex_temp = round(sortrows([di_final; wi_final]),4);
  n_count = groupcounts(design_ex_temp(1,:)'); % this 
  sum(wi_final);
  LOSS(ell,:) = [ell, min(loss1(end)), sum(wi_final), sum(n_count)];
  end
  LOSS_filter = LOSS(round(LOSS(:,3), 3) == 1 & LOSS(:,4) == n(pig), :);
  [M, I] = min(LOSS_filter(:,2));
  my_loss(1, pig) = M;
  my_loss(2, pig) = n(pig);
  my_loss(3, pig) = LOSS_filter(I,3);
  my_loss(4, pig) = LOSS_filter(I,4);
  my_loss(5, pig) = LOSS_filter(I,1);
end


% here, we group the values that are the same together
design_ex_temp = round(sortrows([di_final; wi_final]),4);
val = unique(design_ex_temp(1,:))';
n_count = groupcounts(design_ex_temp(1,:)'); % this 
sum(wi_final);
design_ex = [val, n_count]'

sum(design_ex(2,:))
w_final = design_ex(2, :)/n(pig);
d_final = design_ex(1, :);
FIM_final = zeros(q,q);
for j=1:size(design_ex,2)   
  fx = power(d_final(j), 0:p)';
  FIM_final = FIM_final + (fx * fx') * w_final(j);
end
if criterion == "D"
  Loss_ex = exp(-log(det(FIM_final)^(1/q)));
elseif criterion == "A"
  Loss_ex = trace(inv(FIM_final));
else
  disp("invalid input")
end


my_table = array2table(my_loss,  'RowNames', {'loss', 'n', 'sum of weight', 'n_design', 'seed'});
my_table

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
L_val = [exp(opt_D), opt_A, Loss_ex].';
table(L_val, 'RowNames', {'opt D',  'opt A', 'Exact'})
