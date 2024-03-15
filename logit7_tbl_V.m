%% Compute D-optimal designs for logistic regression
%%with exact n points
 
criterion = "D";
beta = [-0.4926, -0.6280, -0.3283, 0.4378, 0.5283, -0.6120, -0.6837, -0.2061]'; 
S1 = [0, 3];
p = 7; % Dimension
N1 = 4; % Number of design points of each dimension
N = N1^p;
n = 31; 

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
Nt = 500; %number of iterations per temperature change
T0 = 500; %initial temperature
M0 = 15*10000; %number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0*alpha^M0; %minimum temperature

%% ! This part needs to be changed
delta = (S1(2)-S1(1))/(floor(N/4)-1); %neighbourhood size

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
      
      di_nb = zeros(1,p);
      for dog = 1:p
          d_d = unifrnd(-1,1) * delta;
          di_nb(dog) = di(j, dog) + d_d;
          di_nb(dog) = min([max([S1(1), di_nb(dog)]),S1(2)]); 
      end

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



%% Optimal value 
% Define the data
% Define the data
val_paper = [
    1 3 0 3 3 3 0 0 0.0232;
    2 0 3 3 3 3 0 0 0.0262;
    3 0 3 3 0 0 0 3 0.0151;
    4 0 3 3 3 0 3 0 0.0270;
    5 3 3 3 3 0 0 0 0.0191;
    6 3 0 3 0 0 0 0 0.0113;
    7 0 0 3 3 0 3 0 0.0103;
    8 3 0 3 0 0 0 0 0.0215;
    9 0 3 3 0 0 0 3 0.0220;
    10 3 0 0 3 0 0 0 0.0106;
    11 3 0 3 3 0 0 3 0.0366;
    12 0 0 0 3 3 0 3 0.0389;
    13 0 3 3 3 0 0 0 0.0318;
    14 0 3 0 3 0 0 3 0.0309;
    15 0 0 3 3 0 3 3 0.0171;
    16 0 0 3 3 3 0 3 0.0315;
    17 0 3 3 0 0 0 3 0.0385;
    18 0 0 0 0 0 0 0 0.0367;
    19 0 0 3 0 3 0 0 0.0319;
    20 0 0 3 3 3 0 0 0.0323;
    21 0 0 0 3 3 0 0 0.0128;
    22 0 0 0 3 0 0 3 0.0217;
    23 3 3 3 3 0 0 0 0.0333;
    24 0 0 3 3 3 0 3 0.0114;
    25 0 0 0 3 0 3 0 0.0389;
    26 0 0 0 0 0 0 0 0.0107;
    27 0 0 3 0 0 3 0 0.0271;
    28 0 0 3 0 0 0 0 0.0380;
    29 0 3 3 3 3 0 0 0.0263;
    30 3 0 3 3 0 0 3 0.0298;
    31 0 0 3 3 0 3 0 0.0346;
    32 3 0 0 3 0 0 0 0.0316;
    33 0 3 0 3 0 0 0 0.0312;
    34 0 0 3 0 0 0 0 0.0131;
    35 0 0 0 3 0 0 3 0.0225;
    36 0 3 0 3 0 0 0 0.0138;
    37 0 3 3 3 0 3 3 0.0300;
    38 3 0 3 3 0 0 0 0.0100;
    39 0 0 3 3 3 3 0 0.0118;
    40 0 0 3 0 0 0 3 0.0146;
    41 0 0 3 3 0 3 3 0.0242
];



disp(val_paper);
d_paper = val_paper(:,2:8);
w_paper = val_paper(:, 9);
FIM_paper = zeros(q, q);
for j = 1:length(w_paper)
  xx = d_paper(j, :);
  rx = [1, xx]';
  Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
  FIM_paper = FIM_paper + Gamma * (rx * rx') * w_paper(j);
end
L_paper = -log(det(FIM_paper)^(1/q));

L_val = [L00, min(loss1), loss1(end)].';
table(L_val, 'RowNames', {'Approx',  'Exact (min)', 'Exact (end)'})


d_ex = design_ex(:, 1:7);
w_ex = design_ex(:, end)/n;
FIM_ex = zeros(q, q);
for j = 1:length(w_ex)
  xx = d_ex(j, :);
  rx = [1, xx]';
  Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
  FIM_ex = FIM_ex + Gamma * (rx * rx') * w_ex(j);
end
L_ex = -log(det(FIM_ex)^(1/q))