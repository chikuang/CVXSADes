%%Compute D-optimal designs for polynomial regression with exact n points
 
clear;
criterion = "A";
tol = 1E-4; % for finding and filtering out the points 
N = 501;
n = [10, 20 , 30, 40, 50]';
my_loss = zeros(2, size(n,1));
tol_annealing = 1E-40;
% runningtime = cputime;  %record computation time

 % number of design points for initial design
  
a =  -1;   %[a, b] is the design space
b =   1;  
p = 3;            % degree of polynomial regression model  
q = p+1; % how many beta's (degree + 1 intercept term)



% The following vectors and matrices are used in the information
% matrices below.
Nsim = 500;
 if criterion == "D"
    design_true = [-1, -0.447, 0.447, 1; 0.25, 0.25, 0.25,0.25];
  elseif criterion == "A"
    design_true = [-1, -0.464, 0.464, 1; 0.151, 0.349, 0.349, 0.151];
  else
    fprintf('Does not run.');
  end



c0 = 1; % max number of points to be changed in the annealing algorithm
Nt = 200; % number of iterations per temperature change
T0 = 0.1; % initial temperature
M0 = 500; % number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0 * alpha^M0; %minimum temperature

delta = 2*(b-a)/(N-1); % neighbourhood size, in this setting, it is 0.2
w00 = design_true(2,:);
d00 = design_true(1,:);


for pig = 1:size(n,1)
  n(pig)
  LOSS = zeros(Nsim, 3);
  for ell=1:Nsim 
    % disp(ell)
    rng(ell);  %random seed number
    w01 = initializeExact(w00, n(pig)); %convert approximate design lazily to an exact design
    d0 = design_true(1,:);
    w0 = w01;
    
    k = length(w0); 
    
    FIM_temp = FIM_polyP(d0, p);
  
    FIM = zeros(q, q);
    for j=1:size(design_true,2)    
      FIM = FIM + FIM_temp(:,:,j) * w0(j);
    end
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
        wi(j) = wi(j) - 1/n(pig); 
      end
    
      % add weight for new points
      % $ here we add the new (equal) weight to the new points
      for j = 1:ci
        wi(k+j) = wi(k+j) + 1/n(pig); 
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
      
      % compute loss of candidate design
      % FIMi = FIM_polyP(di, p);
      % FIMi = sum(FIMi.*reshape(wi,1,1,[]),3);
      FIMi_temp = FIM_polyP(di, p);
  
      FIMi = zeros(q, q);
      for j=1:length(wi)   
        FIMi = FIMi + FIMi_temp(:,:,j) * wi(j);
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
  loss1 = loss(1:num_iters);
  
  %% 4. PLOTTING RESULTS
  % here, we group the values that are the same together
  design_ex_temp = round(sortrows([d0, w0]),4);
  val = unique(design_ex_temp(:,1));
  n_count = groupcounts(design_ex_temp(:,1));
  sum(w0);
  LOSS(ell,:) = [ell, loss1(end), sum(w0)];
  end
  LOSS_filter = LOSS(round(LOSS(:,3), 3) == 1, :);
  my_loss(1, pig) = min(LOSS_filter(:,2));
  my_loss(2, pig) = n(pig);
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
for j=1:size(design_true,2)
  x = design_true(1, j);
  f = power(x, 0:p)';
  w = design_true(2, j);
  C = C + (f * f') * w;
end

opt_A = trace(inv(C));   %A-opt
