%% Compute D-optimal designs for logistic regression
%%with exact n points
 
clear;
runningtime=cputime;  %record computation time

beta = [-0.2, 0.3, -0.15]';
% beta = [-0.4926, -0.6280, -0.3283]';

rng('default')
S1 = [-1, +1]; % Design space in each dimension
p = 2; % Dimension
N1 = 5; % Number of design points of each dimension
N = N1^p;

%% Generate multi-dimensional space (using some Matlab tricks)
X = cell(1, p);
X(:) = {linspace(S1(1), S1(2), N1)};
u = sortrows(combvec(X{:}).') ; 

%% 0. Initialization
tol = 10^(-4); % for finding and filtering out the points 
tol_annealing = 1E-40;
% N = 21;  % number of design points for initial design
q = length(beta);

%The following vectors and matrices are used in the information
%matrices below.

%% Step 2, calculate ethe optimal design of in S_t
cvx_begin
cvx_precision best
  % variables wk(Np, 1) del(1) 
  variable w(N)
  expression I_design(q,q); 
  % minimize del(1)
  % subject to
  I_design = zeros(q, q) ;
  for i = 1:N
      xx = u(i,:);
      rx = [1, xx]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      I_design = I_design + w(i) * Gamma * (rx * rx');
  end
  
  % -log_det(I_design) <= del;
  minimize(-log_det(I_design))
  ones(N, 1)' * w == 1;
  -w <= zeros(N, 1);
cvx_end
        

% organizing the result
kk = find(w > tol); % because the computer does not have exact zero
design_app = [u(kk,:), w(kk)]; %optimal design 
d00 = design_app(:,1:end-1); % support points
w00 = design_app(:,end); % optimal weight
L00 = cvx_optval; % optimal objective value


%Find n exact design points using an annealing algorithm with the
%following setting
n = 15; 
c0 = 1; %max number of points to be changed in the annealing algorithm
Nt = 200; %number of iterations per temperature change
T0 = 55; %initial temperature
M0 = 1000; %number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0*alpha^M0; %minimum temperature

%% ! This part needs to be changed
delta = 2*(S1(2)-S1(1))/(N-1); %neighbourhood size

w01 = initializeExact(w00, n); %convert approximate design lazily to an exact design
d0 = design_app(1,:);
w0 = w01;


k=length(w0); 

% calculate the FIM and objective function value for this exact design
FIM = zeros(q, q);
for i = 1:size(design_app,2)   
  xx = d0(j);
  rx = [1, xx]';
  Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
  FIM = FIM + w0(i) * Gamma * (rx * rx');
end
L0 = -log(det(FIM)^(1/q));  %D-optimality

% store loss at each iteration for plotting
loss = zeros(1,M0*Nt);
loss(1) = L0;

%% ANNEALING ALGORITTHM

% This algorithm is an implementation of simulated annealing with a constraint
% on the generation of points. New generated points for candidate designs
% are within a predefinied neighbourhood of some random subset of points from the
% previous design. This is to test if an exact design can be further
% optimized by replacing some of the support points with nearby points.

rng(523803);  %random seed number
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
    
    % compute loss of candidate design

    FIMi = zeros(q, q);
    for j=1:size(design_app,2)   
      xx = di(j);
      rx = [1, xx]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      FIMi = FIMi + Gamma * (rx * rx') * w0(j);
    end
    
    Li = -log(det(FIMi)^(1/q));
 
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


%% PLOTTING RESULTS

% Plot loss
figure;
plot(1:num_iters,loss(1:num_iters));
xlabel("Iteration");
ylabel("Loss");
title("Annealing Schedule");
textbox_str = {sprintf('n = %d', n), sprintf('T0 =  %.d', T0),sprintf('M0 = %d', M0),sprintf('Nt = %d', Nt),sprintf('c0 = %d', c0),sprintf('delta = %.2f', delta)};
textbox = annotation('textbox', [0.75, 0.7, 0.2, 0.1], 'String', textbox_str, 'FitBoxToText', 'on', 'EdgeColor', 'none', 'BackgroundColor', 'white');

loss1=loss(1:num_iters);
% here, we group the values that are the same together
design_ex_temp = round(sortrows([d0, w0]),4);
val = unique(design_ex_temp(:,1));
weight = groupcounts(design_ex_temp(:,1));
sum(w0)
design_ex = [val, weight]';

figure; 
% scatter(d0,n*w0,"blue");
histogram(d0, 20, "EdgeColor", "red")
xlabel("Design space");
ylabel("Frequency");
title("Exact design distribution (n = " + sum(design_ex(2,:)) + ")")
ax = gca;
ax.YTick = unique( round(ax.YTick) );

resulttime = cputime-runningtime  %computation time
 

design_app
design_ex
[L00, loss1(1), min(loss1)]    
    