% Athuor: Chi-Kuang Yeh
% July 31, 2018
% Reference: https://arxiv.org/pdf/1701.00888.pdf
% Description: 
%   Computing D-optimal design under OLSE

%% 1.  Input
criterion = 'c';
if criterion == "c"
  cVec = [1, 0, 0]';
end
range = [1, 61]';
a = range(1);
b = range(2);
theta = [0.07, 0.93, 0.96]';
p0 = theta(1); p1 = theta(2); p2 = theta(3);
N = 61;
Nsim = 10;
n = [14]';
tol_annealing = 1E-40;
tol = 1E-4; % for finding and filtering out the points 
%% 2. Initilization
  % discretized equally spaced space
 
u =  linspace(a, b, N);
w = zeros(N+1, 1); del = 0 ; 
one_vec = ones(N,1); zero_vec = zeros(N, 1); % +1 as we want to include 2 boundary points
M = zeros(3,3);
%%% cvx part
cvx_begin quiet
  cvx_precision best
  variables w(N,1) del(1)
  minimize del(1)
  subject to
  % constructing the B matrix
   for i = 1:N
     pi = p1-(p1+p2-1)*(1-p0)^u(i);
     lam = inv(pi*(1-pi));
     % now we have to find f
     f0 = u(i)*(p1+p2-1)*(1-p0)^(u(i)-1);
     f1 = 1 - (1-p0)^u(i);
     f2 = -(1-p0)^u(i);
     f = [f0,f1,f2]';
     M = M + w(i)*lam*f*f';
   end
    % constrains
    if criterion == "D"
      -log(det_rootn(M)) <= del;
    elseif criterion == "A"
       trace_inv(M) <= del;
    elseif criterion == "c"
       matrix_frac(cVec,M) <= del;
    else
      disp("Error");
    end
    -w <= zero_vec;
    one_vec' * w == 1;
cvx_end
  
%%% manage the outputs
kk = find(w>tol); % 
design_app = [u(kk);w(kk)'];
d00 = design_app(1,:); % support points
w00 = design_app(2,:); % optimal weight
L00 = cvx_optval; % optimal objective value
cvx_status


%% 2. Find n exact design points using an annealing algorithm with the
% following setting

c0 = 1; % max number of points to be changed in the annealing algorithm
Nt = 600; % number of iterations per temperature change
T0 = 0.1; % initial temperature
M0 = 1E4; % number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0 * alpha^M0; %minimum temperature

delta = 2*(b-a)/500; % neighbourhood size, in this setting, it is 0.2
q = 3;
my_loss = zeros(5, size(n,1));
for pig = 1:size(n,1)
  LOSS = zeros(Nsim, 4);
  n_i = n(pig);
  disp(n_i)
  for ell=1:Nsim 
    disp(ell)
    rng(ell);  %random seed number
    rng(43);
    % w01 = initializeExact2(w00, n); %convert approximate design lazily to an exact design
    w01 = initializeExact(w00, n_i); %convert approximate design lazily to an exact design
    w0 = w01;
    d0 = design_app(1,:);
    k = length(w0); 
  
    FIM = zeros(q, q);
    for i=1:length(w0)  
      x = design_app(1,i);
        pi = p1-(p1+p2-1)*(1-p0)^x;
       lam = inv(pi*(1-pi));
       % now we have to find f
       f0 = x*(p1+p2-1)*(1-p0)^(x-1);
       f1 = 1 - (1-p0)^x;
       f2 = -(1-p0)^x;
       f = [f0,f1,f2]';
       FIM = FIM + w0(i)*lam*(f * f' );
    end
    
    if criterion == "D"
      L0 = -log(det(FIM)^(1/q));  %D-optimality
    elseif criterion == "A"
      L0 = trace(inv(FIM));  %D-optimalityminimize( trace_inv(A) )   %A-opt
    elseif criterion == "c"
      L0 = cVec' * FIM * cVec;      
    else 
      fprintf('Does not run.')
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
          % di_nb = di(j) + unifrnd(-1,1)*delta; % ! check this neighbour definition

          cat = binornd(1,0.5);
          if cat == 1
            cat = cat;
          else 
            cat = -1;
          end
            di_nb = di(j) + cat;
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
           x = di(j);
           pi = p1-(p1+p2-1)*(1-p0)^x;
           lam = inv(pi*(1-pi));
           % now we have to find f
           f0 = x*(p1+p2-1)*(1-p0)^(x-1);
           f1 = 1 - (1-p0)^x;
           f2 = -(1-p0)^x;
           f = [f0,f1,f2]';
           FIMi = FIMi + wi(j)*lam*(f * f' );
        end
    
      if criterion == "D"
        Li = -log(det(FIMi)^(1/q));  %D-optimality
      elseif criterion == "A"
        Li = trace(inv(FIMi));  %D-optimalityminimize( trace_inv(A) )   %A-opt
      elseif criterion == "c"
        %Li = cVec' * FIMi * cVec;   
        Li=matrix_frac(cVec,FIMi);
      else 
        fprintf('Does not run.')
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
    
    %% 4. PLOTTING RESULTS
    % here, we group the values that are the same together
    design_ex_temp = round(sortrows([d0, w0]),4);
    val = unique(design_ex_temp(:,1));
    n_count = groupcounts(design_ex_temp(:,1));
    sum(w0);
    design_ex = [val, n_count]';
    LOSS(ell,:) = [ell, min(loss1(end)), sum(wi), sum(n_count)];
  end
  LOSS_filter = LOSS(round(LOSS(:,3), 3) == 1 & LOSS(:,4) == n_i, :);
  [M, I] = min(LOSS_filter(:,2));
  my_loss(1, pig) = M;
  my_loss(2, pig) = n_i;
  my_loss(3, pig) = LOSS_filter(I,3); 
  my_loss(4, pig) = LOSS_filter(I,4);
  my_loss(5, pig) = LOSS_filter(I,1);
end
my_table = array2table(my_loss, ...
  'RowNames', {'loss', 'n', 'sum of weight', 'n_design', 'seed'});
my_table
sortrows(LOSS, 2, "descend");
L00
min(LOSS(:,2))

design_app
design_ex
loss1(end)
L00
Effc=L00/loss1(end)
 
LOSS

plot(1:length(loss1), loss1) %annealing loss function plot
xlabel('iteration')
ylabel('Loss function')
%ylim([exp(-2) exp(-1.0)])