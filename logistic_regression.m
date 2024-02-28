%% Compute D-optimal designs for logistic regression
%%with exact n points
 
criterion = "A";

% beta = [-0.5, 0.7, 0.38]';
% beta = [-0.4926, -0.6280, -0.3283]';
% beta = [3, 0.5, 0.75]';
beta = [-0.550, -0.157, 0.194]'; % values used in p186, Berger and Wong.
S1 = [-40, +40]; 
S2 = [-40, +40]; 
p = 2; % Dimension
N1 = 201; % Number of design points of each dimension
N = N1^p;
n = 15; 

%% Generate multi-dimensional space (using some Matlab tricks)
X = cell(1, p);
X(:) = {linspace(S1(1), S1(2), N1)};
u = sortrows(combvec(X{:}).') ; 

%% 0. Initialization
tol = 10^(-4); % for finding and filtering out the points 
tol_annealing = 1E-40;
q = length(beta);

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
  for i = 1:N
      xx = u(i,:);
      rx = [1, xx]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      M = M + w(i) * Gamma * (rx * rx');
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
cvx_status


%Find n exact design points using an annealing algorithm with the
%following setting
c0 = 1; %max number of points to be changed in the annealing algorithm
Nt = 200; %number of iterations per temperature change
T0 = 500; %initial temperature
M0 = 10000; %number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0*alpha^M0; %minimum temperature

%% ! This part needs to be changed
delta = (S1(2)-S1(1))/(N1-1); %neighbourhood size

w01 = initializeExact(w00, n); %convert approximate design lazily to an exact design
d0 = design_app(:,1:(end-1));
w0 = w01;


k = length(w0); 

% calculate the FIM and objective function value for this exact design
FIM = zeros(q, q);
for j = 1:size(design_app, 1)   
  xx = d0(j, :);
  rx = [1, xx]';
  Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
  FIM = FIM + w0(j) * Gamma * (rx * rx');
end
if criterion == "D"
  L0 = -log(det(FIM)^(1/q));
elseif criterion == "A"
  L0 = trace(inv(FIM));   %A-opt
else
  fprintf('Does not run.');
end

  %D-optimality

% store loss at each iteration for plotting
loss = zeros(1, M0*Nt);
loss(1) = L0;

%% ANNEALING ALGORITTHM

% This algorithm is an implementation of simulated annealing with a constraint
% on the generation of points. New generated points for candidate designs
% are within a predefinied neighbourhood of some random subset of points from the
% previous design. This is to test if an exact design can be further
% optimized by replacing some of the support points with nearby points.

% rng(523803);  %random seed number
num_iters = 1;
T = T0;
L_prev = 0;

while(T > Tmin && abs(L_prev - L0) > tol_annealing )
  L_prev = L0;
  % GENERATING RANDOM CANDIDATE DESIGN
  for h = 1:Nt
    k = length(w0); % current number of support points
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
      d_d1 = unifrnd(-1,1) * delta;
      d_d2 = unifrnd(-1,1) * delta;
      
      di_nb = zeros(1,2);
      di_nb(1) = di(j,1) + d_d1;
      di_nb(1) = min([max([S1(1), di_nb(1)]),S1(2)]); 
      di_nb(2) = di(j,2) + d_d2;
      di_nb(2) = min([max([S2(1), di_nb(2)]),S2(2)]);
      % di_nb = di(j) + unifrnd(-1,1)*delta; % ! check this neighbour definition

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
      xx = di(j, :);
      rx = [1, xx]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      FIMi = FIMi + Gamma * (rx * rx') * wi(j);
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


%% 4. PLOTTING RESULTS

% Plot loss
% figure;
% plot(1:num_iters,loss(1:num_iters));
% xlabel("Iteration");
% ylabel("Loss");
% title("Annealing Schedule");

loss1 = loss(2:num_iters);

% here, we group the values that are the same together
design_ex_temp = round(sortrows([d0, w0]),4);
% table(design_ex)
val = unique(design_ex_temp(:,1:(end-1)),'rows');
n_count = groupcounts(design_ex_temp(:, 1:(end-1)));
sum(w0)
design_ex = [val, n_count];
plot(loss1)
design_app 
design_ex
sum(design_ex(:, 3))
L_val = [L00, min(loss1), loss1(end)].';
table(L_val, 'RowNames', {'Approx',  'Exact (min)', 'Exact (end)'})