%% Compute D-optimal designs for logistic regression with interaction

%% Reference paper: 
% Haines, L. M., & Kabera, G. M. (2018). D-optimal designs 
% for the two-variable binary logistic regression model with interaction. Journal of Statistical Planning and Inference, 193, 136-150.
 

criterion = "D";
% beta = [1, 2, 2, 0.2]';
%beta = [-3, 4, 6, 1]' ; %Example 4.2 (b)
beta = [-2.2054, 13.5803, 2.2547, 1.6262]' ; %A real example
S1 = [0, 2]; 
S2 = [0, 2]; 
p = 2; % Dimension
Nsim = 20;
N1 = 51; % Number of design p       oints of each dimension
N = N1^p;
% n = [5, 10, 15, 20]'; 
n = 10

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
      rx = [1, xx, xx(1) * xx(2)]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      M = M + w(i) * Gamma * (rx * rx');
  end
  
   if criterion == "D"
    minimize (-log(det_rootn(M)));
  elseif criterion == "A"
    minimize( trace_inv(M) );   %A-opt
   elseif criterion == "E"
     minimize(-lambda_min(M))
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
T0 = 5; %initial temperature
M0 = 10000; %number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0*alpha^M0; %minimum temperature

%% ! This part needs to be changed
delta = (S1(2)-S1(1))/(N1-1); %neighbourhood size



FIM = zeros(q, q);
for j = 1:size(design_app, 1)   
  xx = d00(j, :);
   rx = [1, xx, xx(1) * xx(2)]';
  Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
  FIM = FIM + w00(j) * Gamma * (rx * rx');
end
if criterion == "D"
  L0 = -log(det(FIM)^(1/q));
elseif criterion == "A"
  L0 = trace(inv(FIM));   %A-opt
elseif criterion == "E"
  L0 = -lambda_min(FIM);
else
  fprintf('Does not run.');
end


my_loss = zeros(5, size(N,1));
for pig = 1:size(n,1)
  LOSS = zeros(Nsim, 4);
  n_i = n(pig);
  disp(n_i)
 
  for ell=1:Nsim 
    disp(ell)
    rng(ell);  %random seed number
    rng(2);
    loss = zeros(1, M0*Nt);
    loss(1) = L0;
    w01 = initializeExact(w00, n_i); %convert approximate design lazily to an exact design
    d0 = design_app(:,1:(end-1));
    w0 = w01;
    k = length(w0); 
    
% calculate the FIM and objective function value for this exact design
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
        rx = [1, xx, xx(1) * xx(2)]';
        Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
        FIMi = FIMi + Gamma * (rx * rx') * wi(j);
      end
  
  
       if criterion == "D"
        Li = -log(det(FIMi)^(1/q));
      elseif criterion == "A"
        Li = trace(inv(FIMi));   %A-opt
      elseif criterion == "E"
        Li = -lambda_min(FIMi);
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
    val = unique(design_ex_temp(:,1:(end-1)),'rows');
    n_count = groupcounts(design_ex_temp(:,1:(end-1)));
    sum(w0);
    LOSS(ell,:) = [ell, loss1(end),sum(w0), sum(n_count)];
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
plot(exp(loss1))
design_app 
design_ex
LOSS
Eff=exp(L00)/exp(loss1(end))
exp(L00)
exp(loss1(end))

sum(design_ex(:, 3)) % number of points
L_val = [L00, min(loss1), loss1(end)].';
table(L_val, 'RowNames', {'Approx',  'Exact (min)', 'Exact (end)'})