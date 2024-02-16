clear;
criterion = "D";

%% 0. Initialization
tol = 1E-4; % for finding and filtering out the points 
tol_annealing = 1E-40;
Nsim = 20;

% S1 = [0, 5]; % Design space in each dimension
% S1 = [175, 275];
% S2 = [10, 30];
S1 = [-2, +2];
S2 = [-2, +2];
p = 2; % Dimension
N1 = 21; % Number of design points of each dimension
N2 = 21;
N = N1*N2;
%% Generate multi-dimensional space (using some Matlab tricks)
X = cell(1, p);
X(1) = {linspace(S1(1), S1(2), N1)};
X(2) = {linspace(S2(1), S2(2), N2)};
u = sortrows(combvec(X{:}).') ; 
q = 6;
%% 1. Compute the initial proxy approximate design 
cvx_begin
  cvx_precision high
  variable w(1,N);
  expression A(q,q); 
  
  % here we compute the information matrix at
  for j=1:N
    xj = u(j,:);
    fj = second_order(xj);
    A = A + (fj * fj') * w(j);
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
kk = find(w > tol); % because the computer does not have exact zero
design_app = [u(kk,:), w(kk)']; %optimal design 
d00 = design_app(:,1:end-1); % support points
w00 = design_app(:,end); % optimal weight
L00 = cvx_optval; % optimal objective value

%% 2. Find n exact design points using an annealing algorithm with the
% following setting
%% 2. Find n exact design points using an annealing algorithm with the
% following setting
n = 9; % number of support point in the exact design
c0 = 1; % max number of points to be changed in the annealing algorithm
Nt = 200; % number of iterations per temperature change
T0 = 0.1; % initial temperature
M0 = 500; % number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0 * alpha^M0; %minimum temperature

delta = 2*(S1(2)-S1(1))/(N-1); %neighbourhood size




LOSS = zeros(Nsim, 2);

for ell=1:Nsim 
  disp(ell)
  % rng(ell);  %random seed number
  % rng(202)
  % w01 = initializeExact2(w00, n); % convert approximate design, but
  % retain the number of design points
  w01 = initializeExact(w00, n); %convert approximate design lazily to an exact design
  d0 = design_app(:,1:(end-1));
  w0 = w01;
  
  k = length(w0); 
  
  FIM = zeros(q, q);
  for j=1:size(design_app,1)    
    xj = d0(j,:);
    fj = second_order(xj);
    FIM = FIM + (fj * fj') *  w0(j);
  end
  
  if criterion == "D"
    L0 =  -log(det(FIM)^(1/q));
  elseif criterion == "A"
    L0 = trace(inv(FIM));   %A-opt
  else
    fprintf('Does not run.');
  end
  
  
  % L0 = -log(det(FIM)^(1/q));  %D-optimality
  
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
      k = length(w0); % current number of support points
      ci =  min([k, randi([1 c0])]); % number of points to be replaced
      toRemove = randperm(k, ci); % randomly select ci points to remove, return an index
      
      % creating new candidate design
      wi = zeros(k+ci, 1); 
      di = zeros(size(d0,1)+ci, size(d0,2)); % ! This will be depends on the dimension
    
      
      
      % copy previous design
      wi(1:k) = w0;
      di(1:k,:) = d0;
      
      % remove selected points, here, there may be multiple points
      % $ here we essentially decrease the weight by 1/n (remove one point)
      for j = toRemove
        wi(j) = wi(j) - 1/n; 
      end
    
      % add weight for new points
      % $ here we add the new (equal) weight to the new points
      % $ as we know that in D-opt, the weight are equal
      for j = 1:ci
        wi(k+j) = wi(k+j) + 1/n; 
      end
      
      % generate new points in neighbourhood of those previously removed
      i = 1;
      for j = toRemove
        % $ delta is the neighbour
        dx = unifrnd(-1,1)*delta;
        dy = unifrnd(-1,1)*delta;
        
        di_nb = zeros(2,1);
        di_nb(1) = di(j,1) + dx;
        di_nb(1) = min([max([S1(1), di_nb(1)]),S1(2)]);
        di_nb(2) = di(j,2) + dy;
        di_nb(2) = min([max([S1(1), di_nb(2)]),S1(2)]);
        % di_nb = di(j) + unifrnd(-1,1)*delta; % ! check this neighbour definition
  
        % ! actually here we can use the rejection proposal kind of approach
        di(k+i,:) = di_nb;
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
      for j=1:length(wi)   
        xj = di(j,:);
        fj = second_order(xj);
        FIMi = FIMi + (fj * fj') * wi(j);
      end
      
      if criterion == "D"
        Li =  -log(det(FIMi)^(1/q));
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
  loss1 = loss(1:num_iters);
  
  %% 4. PLOTTING RESULTS
  % here, we group the values that are the same together
  design_ex_temp = round(sortrows([d0, w0]),4);
  % table(design_ex)
  val = unique(design_ex_temp(:,1:(end-1)),'rows');
  n_count = groupcounts(design_ex_temp(:,1:(end-1)));
  % sum(w0)
  design_ex = [val, n_count];
  LOSS(ell,:) = [ell, loss1(end)];
end

sortrows(LOSS, 2, "descend");
min(LOSS(:,2))
L00