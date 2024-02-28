clear;
criterion = "D";

n = 10; % number of support point in the exact design
N = 1001;  % number of design points for initial design
Nsim = 500;
LOSS = zeros(Nsim, 2);

beta = [0.9506, 2.0434]';
a =  0;  
b =  4;

%% 0. Initialization
tol = 1E-4; % for finding and filtering out the points 
tol_annealing = 1E-40;

 
q = length(beta);  
% q = p+1; % how many beta's (degree + 1 intercept term)
u = linspace(a, b, N); %equally spaced N points in [a,b]

%% 1. Compute the initial proxy approximate design 
cvx_begin
  cvx_precision best
  variable w(1,N);
  expression M(q,q); 
  
  % here we compute the information matrix at
   for i = 1:N
      xx = u(i);
      rx = [1, xx]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      M = M + w(i) * Gamma * (rx * rx');
   end
  
  minimize (-log(det_rootn(M)));
  0 <= w <= 1;
  sum(w)==1;
cvx_end

% organizing the result
design_app = [u(find(w(1,:)>tol)); w(find(w(1,:)>tol))]; %optimal design 
d00 = design_app(1,:); % support points
w00 = design_app(2,:); % optimal weight
L00 = cvx_optval; % optimal objective value

%% 2. Find n exact design points using an annealing algorithm with the
% following setting

c0 = 1; % max number of points to be changed in the annealing algorithm
Nt = 200; % number of iterations per temperature change
T0 = 0.1; % initial temperature
M0 = 500; % number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0 * alpha^M0; %minimum temperature

delta = 2*(b-a)/(N-1); % neighbourhood size, in this setting, it is 0.2

for ell=1:Nsim 
  disp(ell)
  rng(ell);  %random seed number
  % rng(202)
  % w01 = initializeExact2(w00, n); %convert approximate design lazily to an exact design
  w01 = initializeExact(w00, n); %convert approximate design lazily to an exact design
  d0 = design_app(1,:);
  w0 = w01;
  
  
  k = length(w0); 

  FIM = zeros(q, q);
  for j = 1:size(design_app, 1)   
    xx = d0(j);
    rx = [1, xx]';
    Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
    wj = w0(j);
    FIM = FIM + wj * Gamma * (rx * rx');
  end
  
 
    L0 = -log(det(FIM)^(1/q));  %D-optimality
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
      di = zeros(k+ci, 1);
      
      % copy previous design
      wi(1:k) = w0;
      di(1:k) = d0;
      
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
        di_nb = di(j) + unifrnd(-1,1)*delta; % ! check this neighbour definition
  
        % di_nb = di(j) + (2*rand(1)-1)*delta; % ! check this neighbour definition
        di_nb = min([max([a, di_nb]),b]); % this needs more work for high-D 
        % ! actually here we can use the rejection proposal kind of approach
        di(k+i,:) = di_nb;
        i = i + 1; 
      end
      % points are generated in a circular, spherical, or hyperspherical
      % region (depending on the dimension of the design) of radius
      % delta
      
      % remove support points with zero weight
      di = di(wi>tol);
      wi = wi(wi>tol);
        
      FIMi = zeros(q, q);
      for j=1:length(wi)
        xx = di(j);

        rx = [1, xx]';
        Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
        wj = wi(j);
        FIMi = FIMi + wj * Gamma * (rx * rx');
      end
      
      Li = -log(det(FIMi)^(1/q));  %D-optimality
        
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
  val = unique(design_ex_temp(:,1));
  n_count = groupcounts(design_ex_temp(:,1));
  sum(w0);
  
  design_ex = [val, n_count]';
  LOSS(ell,:) = [ell, loss1(end)];
end

design_app
design_ex
L_val = [L00, loss1(1), loss1(end)].';
rowname = {'Approx', 'Exact (initial)', 'Exact (Final)'}.';
table(L_val, 'RowNames', rowname)
