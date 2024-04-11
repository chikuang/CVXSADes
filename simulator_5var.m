%% Compute optimal designs with the model
% y = x1 + x2 + x3 + x4 +x5
clear; clc;
criterion = "D";
S = [ -1, 1];
p = 5; % Dimension
Nsim = 10;
N1 = 5; % Number of design points of each dimension
N = N1^p;
n = [9]'; %[8, 9,10]'; 
delta = 1E-2;
%% Generate multi-dimensional space (using some Matlab tricks)
X = cell(1, p);
X(:) = {linspace(S(1), S(2), N1)};
u = sortrows(combvec(X{:}).') ; 

%% 0. Initialization
tol = 10^(-4); % for finding and filtering out the points 
tol_annealing = 1E-40;
% q = 5;
q = 6; % with intercept

%The following vectors and matrices are used in the information
%matrices below.

%% Step 2, calculate ethe optimal design of in S_t
cvx_begin quiet
cvx_precision best
  % variables wk(Np, 1) del(1) 
  variable w(N)
  expression M(q,q); 
  % minimize del(1)
  % subject to
  for j = 1:N
      xx = [1, u(j, :)]';
      M = M + w(j) * (xx * xx');
  end
  
   if criterion == "D"
    minimize (-log(det_rootn(M)));
  elseif criterion == "A"
    minimize( trace_inv(M) );   %A-opt
  else
    fprintf('Does not run.');
  end
  ones(N, 1)' * w == 1;
  -w <= zeros(N, 1);
cvx_end
        

% organizing the result
kk = find(w > tol); % because the computer does not have exact zero
design_app = [u(kk,:), w(kk)]; %optimal design 
d00 = design_app(:, 1:end-1); % support points
w00 = design_app(:, end); % optimal weight
L00 = cvx_optval; % optimal objective value


%Find n exact design points using an annealing algorithm with the
%following setting
c0 = 1; %max number of points to be changed in the annealing algorithm
Nt = 500; %number of iterations per temperature change
T0 = 500; %initial temperature
M0 = 1E4; %number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0*alpha^M0; %minimum temperature



my_loss = zeros(5, size(n,1));

for pig = 1:size(n,1)
  LOSS = zeros(Nsim, 4);
  n_i = n(pig);
  disp(n_i)
 
  for ell=1:Nsim 
    % disp(ell)
    rng(ell);  %random seed number
    loss = zeros(1, M0*Nt);
    FIM = zeros(q, q);
    for j = 1:length(d00)  
      xx = [1, d00(j, :)]';
      FIM = FIM + w00(j) * (xx * xx');
    end
    if criterion == "D"
      L0 = -log(det(FIM)^(1/q));
    elseif criterion == "A"
      L0 = trace(inv(FIM));   %A-opt
    else
      fprintf('Does not run.');
    end

    loss(1) = L0;
    w01 = initializeExact3(w00, n_i); %convert approximate design lazily to an exact design
    
    if length(w01) < size(design_app,1)
      d0 = design_app(1:length(w01), 1:(end-1));
    else 
      d0 = design_app(:,1:(end-1));
    end
    w0 = w01;
    k = length(w0); 
    
% calculate the FIM and objective function value for this exact design
  FIM1 = zeros(q, q);
  for j = 1:length(w0)
    xx = [1, d0(j, :)]';
    FIM1 = FIM1 + w0(j) * (xx * xx');
  end
  if criterion == "D"
    L0 = -log(det(FIM1)^(1/q));
  elseif criterion == "A"
    L0 = trace(inv(FIM1));   %A-opt
  else
    fprintf('Does not run.');
  end
  loss(2) = L0;


    num_iters = 2;
    T = T0;
    L_prev = 0;
  while(T > Tmin && abs(L_prev - L0) > tol_annealing )
    L_prev = L0;
    % GENERATING RANDOM CANDIDATE DESIGN
    for h = 1:Nt
      k = length(w0); % current number of support point
      ci =  min([k, randi([1 c0])]); % number of points to be replaced
      toRemove = randperm(k, ci); % randomly select ci points to remove, return an index
      
      % creating new candidate design
      wi = zeros(k+ci, 1); 
      di = zeros(size(d0,1) + ci, size(d0, 2)); % ! This will be depends on the dimension
      
      % copy previous design
      wi(1:k) = w0;
      di(1:k,:) = d0;
      
      % remove selected points, here, there may be multiple points
      % $ here we essentially decrease the weight by 1/n (remove one point)
      for j = toRemove
        wi(j) = wi(j) - 1/n_i; 
      end
    
      % add weight for new points
      % $ here we add the new (equal) weight to the new points
      % $ as we know that in D-opt, the weight are equal
      for j = 1:ci
        wi(k+j) = wi(k+j) + 1/n_i; 
      end
      
      % generate new points in neighbourhood of those previously removed
      i = 1;
      for j = toRemove
        % $ delta is the neighbour
        di_nb = zeros(1,p);
        for dog = 1:p
            d_d = unifrnd(-1,1) * delta;
            di_nb(dog) = di(j, dog) + d_d;
            di_nb(dog) = min([max([S(1), di_nb(dog)]),S(2)]); 
        end
   
        % ! actually here we can use the rejection proposal kind of approach
        di(k+i, :) = di_nb;
        i = i + 1; 
      end
      % points are generated in a circular, spherical, or hyperspherical
      % region (depending on the dimension of the design) of radius
      % delta
      
      % remove support points with zero weight
      di = di(wi>tol, :);
      wi = wi(wi>tol);
      % compute loss of candidate design
  
      FIMi = zeros(q, q);
      for j = 1:length(wi)
        xx = [1, di(j, :)]';
        FIMi = FIMi + wi(j) * (xx * xx');
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
        d0 = di;
        w0 = wi;
        num_iters = num_iters + 1;
        loss(num_iters) = L0;
      end
    end
    T = alpha*T;
  end
    loss1 = loss(2:num_iters);
    design_ex_temp = round(sortrows([d0, w0]),4);
    % val = unique(design_ex_temp(:,1:(end-1)),'rows');
    n_count = groupcounts(design_ex_temp(:,1:(end-1)));
    sum(w0);
    % if size(loss1, 2) == 0
    %   loss1 = Inf;
    % end
    LOSS(ell,:) = [ell, loss1(end), sum(w0), sum(n_count)];
  end
  LOSS_filter = LOSS(round(LOSS(:,3), 3) == 1 & LOSS(:,4) == n_i, :);
  [M, I] = min(LOSS_filter(:,2));
  my_loss(1, pig) = M;
  my_loss(2, pig) = n_i;
  my_loss(3, pig) = LOSS_filter(I,3);
  my_loss(4, pig) = LOSS_filter(I,4);
  my_loss(5, pig) = LOSS_filter(I,1);
end

%% 4. PLOTTING RESULTS

% here, we group the values that are the same together
design_ex_temp = round(sortrows([d0, w0]),4);
val = unique(design_ex_temp(:,1:(end-1)),'rows');
n_count = groupcounts(design_ex_temp(:, 1:(end-1)));
sum(w0)
design_ex = [val, n_count];
plot(loss1)
design_app 
design_ex
sum(design_ex(:, end)) % number of points
L_val = [L00, min(loss1), loss1(end)].';
table(L_val, 'RowNames', {'Approx',  'Exact (min)', 'Exact (end)'})

my_table = array2table(my_loss, ...
  'RowNames', {'loss', 'n', 'sum of weight', 'n_design', 'seed'});
my_table