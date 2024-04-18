%Compute maximin D-optimal design in Application 3
%%Four dose-response models 

clear;
criterion = "A";
tol = 10^(-4);
tol_annealing = 1E-40;
Nsim = 1;
N = 201;         %number of design points   
n = 20;
a = 0;  b = 500; %[a, b] is the design space
u = linspace(a,b,N); %equally spaced N points in [a,b]
v0 = 60; v1 = 294; v2=25;   %true parameter values for Emax I model 
v02 = 60;  v12 = 340; v22 = 107.14; %true parameter values for Emax II model 
v03 = 49.62; v13 = 290.51; v23 = 150; v33 = 45.51; %parameter values for logistic model
  
%Vectors and matrices are used in the information matrices below.
p1=2; p2=3; p3=p2; p4=4; %# of parameters in the four models 
F1i=zeros(p1,p1,N); F2i=zeros(p2,p2,N); F3i=zeros(p3,p3,N); F4i=zeros(p4,p4,N);  
for j=1:N 
  d=u(:,j);
  f1=[1 d];
  f2=[1 d/(v2+d) -v1*d/(v2+d)^2];
  f3=[1 d/(v22+d) -v12*d/(v22+d)^2];
  g=exp((v23-d)/v33);
  f4=[1 1/(1+g) -v13*g/v33/((1+g)^2)  v13*g*(v23-d)/v33^2/((1+g)^2)];
  F1i(:,:,j)=(f1'*f1);
  F2i(:,:,j)=(f2'*f2);
  F3i(:,:,j)=(f3'*f3);
  F4i(:,:,j)=(f4'*f4);
end
  
%Compute the D-optimal design for linear model
cvx_begin
    cvx_precision high
    variable w(1,N);
    expression A1(p1,p1); 
    for j=1:N
        A1 = A1+F1i(:,:,j)*w(j);
    end        

     if criterion == "D"
      minimize( - det_rootn(A1) )  
    elseif criterion == "A"
      minimize( trace_inv(A1) );
    else
      fprintf('Does not run.');
     end
    0 <= w <= 1;
    sum(w)==1;
cvx_end

  if criterion == "D"
      Loss1 = 1/(det(A1))^(1/p1);
    elseif criterion == "A"
      Loss1 = trace(inv(A1));   %A-opt
    else
      fprintf('Does not run.');
  end
%design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   
    
%Compute the D-optimal design for Emax I model
cvx_begin
    cvx_precision high
    variable w(1,N);
    expression A2(p2,p2); 
    for j=1:N
        A2 = A2+F2i(:,:,j)*w(j);
    end
     if criterion == "D"
      minimize( - det_rootn(A2) )  
    elseif criterion == "A"
      minimize( trace_inv(A2) );   %A-opt
    else
      fprintf('Does not run.');
     end
    0 <= w <= 1;
    sum(w)==1;
cvx_end
if criterion == "D"
    Loss2=1/(det(A2))^(1/p2);
  elseif criterion == "A"
    Loss2 = trace(inv(A2));   %A-opt
  else
    fprintf('Does not run.');
end

%design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   

%Compute the D-optimal design for Emax II model
cvx_begin
    cvx_precision high
    variable w(1,N);
    expression A3(p3,p3); 

    for j=1:N
        A3 = A3+F3i(:,:,j)*w(j);
    end               
     if criterion == "D"
      minimize( - det_rootn(A3) );
    elseif criterion == "A"
      minimize( trace_inv(A3) );
    else
      fprintf('Does not run.');
     end 
    0 <= w <= 1;
    sum(w)==1;
cvx_end
if criterion == "D"
    Loss3=1/(det(A3))^(1/p3);
  elseif criterion == "A"
    Loss3 = trace(inv(A3));   %A-opt
  else
    fprintf('Does not run.');
end
%design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   

%Compute the D-optimal design for logistic model
cvx_begin
    cvx_precision high
    variable w(1,N);
    expression A4(p4,p4); 
    for j=1:N
        A4 = A4+F4i(:,:,j)*w(j);
    end
    if criterion == "D"
      minimize( - det_rootn(A4) );
    elseif criterion == "A"
      minimize( trace_inv(A4) );;   %A-opt
    else
      fprintf('Does not run.');
     end 
    0 <= w <= 1;
    sum(w)==1;
cvx_end
if criterion == "D"
    Loss4=1/(det(A4))^(1/p4);
  elseif criterion == "A"
    Loss4 = trace(inv(A4));   %A-opt
  else
    fprintf('Does not run.');
end
%design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']   

%Compute the Maximin D-efficiency design for the 4 models
cvx_begin
    cvx_precision high
    variable w(1,N+1);
    expression B1(p1,p1); 
    expression B2(p2,p2);
    expression B3(p3,p3); 
    expression B4(p4,p4); 
    for j=1:N
        B1 = B1+F1i(:,:,j)*w(j);
        B2 = B2+F2i(:,:,j)*w(j);
        B3 = B3+F3i(:,:,j)*w(j);
        B4 = B4+F4i(:,:,j)*w(j);
    end         
    

    if criterion == "D"
      maximize w(N+1)
      det_rootn(B1)- w(N+1)/Loss1 >=0; 
      det_rootn(B2)- w(N+1)/Loss2 >=0; 
      det_rootn(B3)- w(N+1)/Loss3 >=0; 
      det_rootn(B4)- w(N+1)/Loss4 >=0; 
    elseif criterion == "A"
      minimize w(N+1)
      trace_inv(B1) - w(N+1)*Loss1 <=0
      trace_inv(B2) - w(N+1)*Loss2 <=0
      trace_inv(B3) - w(N+1)*Loss3 <=0
      trace_inv(B4) - w(N+1)*Loss4 <=0
    else
      fprintf('Does not run.');
     end 
    0 <= w;
    sum(w(1:N))==1;
cvx_end

if criterion == "D"
    Loss1d=1/(det(B1))^(1/p1);
    Loss2d=1/(det(B2))^(1/p2);
    Loss3d=1/(det(B3))^(1/p3);
    Loss4d=1/(det(B4))^(1/p4);
  elseif criterion == "A"
    Loss1d = trace(inv(B1));
    Loss2d = trace(inv(B2));
    Loss3d = trace(inv(B3));
    Loss4d = trace(inv(B4));
  else
    fprintf('Does not run.');
end


eff1=Loss1/Loss1d;  %Efficiency at the maximin design
eff2=Loss2/Loss2d;
eff3=Loss3/Loss3d;
eff4=Loss4/Loss4d;
w=w(1:N);
design_app =[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))']'; %Maximin design   
    

L00 = min([eff1, eff2, eff3,eff4]) ;
eff_app = [eff1, eff2, eff3,eff4] ;
L0 = -L00;


%% Start the annealing part
c0 = 1; % max number of points to be changed in the annealing algorithm
Nt = 200; % number of iterations per temperature change
T0 = 0.1; % initial temperature
M0 = 500; % number of temperature changes before algorithm stops
alpha = 0.9;  %temperature cooling rate

Tmin = T0 * alpha^M0; %minimum temperature

delta = 2*(b-a)/(N-1); % neighbourhood size, in this setting, it is 0.2
my_loss = zeros(5, size(N,1));
LOSS = zeros(Nsim, 4);

for pig = 1:size(n,1)
  LOSS = zeros(Nsim, 4);
  n_i = n(pig);
  disp(n_i)
  
  for ell=1:Nsim 
    % disp(ell)
    ell = 94
    rng(ell);  %random seed number
    
    loss = zeros(1, M0*Nt);
    loss(1) = L0;
    %% Annealing
    c0 = 1; % max number of points to be changed in the annealing algorithm
    Nt = 200; % number of iterations per temperature change
    T0 = 0.1; % initial temperature
    M0 = 500; % number of temperature changes before algorithm stops
    alpha = 0.9;  %temperature cooling rate
    
    Tmin = T0 * alpha^M0; %minimum temperature
    
    delta = 2*(b-a)/(N-1); % neighbourhood size, in this setting, it is 0.2
    w00 = design_app(2,:);
    d00 = design_app(1,:);
    
    
    w01 = initializeExact(w00, n_i); %convert approximate design lazily to an exact design
    d0 = design_app(1,:);
    w0 = w01;
    
    loss = zeros(M0*Nt,1);
    loss(1) = L0;
    num_iters = 1;
    T = T0;
    
    L_prev = 0;
    
    while(T > Tmin && abs(L_prev -L0) > tol_annealing )
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
          di_nb = di(j) + unifrnd(-1,1)*delta; % ! check this neighbour definition
          di_nb = min([max([a, di_nb]),b]); % this needs more work for high-D 
          % ! actually here we can use the rejection proposal kind of approach
          di(k+i) = di_nb;
          i = i + 1; 
        end
        
        % remove support points with zero weight
        di = di(wi>tol);
        wi = wi(wi>tol);
        
        kk = length(wi);
        % compute loss of candidate design
        F1i=zeros(p1,p1,kk); F2i=zeros(p2,p2,kk); F3i=zeros(p3,p3,kk); F4i=zeros(p4,p4,kk);  
        for j=1:kk
          d=di(j);
          f1=[1 d];
          f2=[1 d/(v2+d) -v1*d/(v2+d)^2];
          f3=[1 d/(v22+d) -v12*d/(v22+d)^2];
          g=exp((v23-d)/v33);
          f4=[1 1/(1+g) -v13*g/v33/((1+g)^2)  v13*g*(v23-d)/v33^2/((1+g)^2)];
          F1i(:,:,j)=(f1'*f1);
          F2i(:,:,j)=(f2'*f2);
          F3i(:,:,j)=(f3'*f3);
          F4i(:,:,j)=(f4'*f4);
        end
        
    
        B1=zeros(p1,p1); B2=zeros(p2,p2); B3=zeros(p3,p3); B4=zeros(p4,p4);
        for j=1:kk
          B1 = B1+F1i(:,:,j)*wi(j);
          B2 = B2+F2i(:,:,j)*wi(j);
          B3 = B3+F3i(:,:,j)*wi(j);
          B4 = B4+F4i(:,:,j)*wi(j);
        end  
        if criterion == "D"
          Loss1d=1/(det(B1))^(1/p1);
          Loss2d=1/(det(B2))^(1/p2);
          Loss3d=1/(det(B3))^(1/p3);
          Loss4d=1/(det(B4))^(1/p4);
        elseif criterion == "A"
          Loss1d = trace(inv(B1));
          Loss2d = trace(inv(B2));
          Loss3d = trace(inv(B3));
          Loss4d = trace(inv(B4));
        else
          fprintf('Does not run.');
        end
        eff1=Loss1/Loss1d;  %Efficiency at the maximin design
        eff2=Loss2/Loss2d;
        eff3=Loss3/Loss3d;
        eff4=Loss4/Loss4d;
        Li = -min([eff1, eff2, eff3,eff4]) ;
    
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
    val = unique(design_ex_temp(:,1));
    n_count = groupcounts(design_ex_temp(:,1));
    sum(w0);
    
    design_ex = [val, n_count];
    design_ex;
    
    %% Calculate the final loss
    
    d_final = design_ex(:,1);
    n_final = design_ex(:,2);
    n_total = sum(n_final);
    n_pt = size(design_ex,1);
    w_final = n_final/n_total;
    F1i=zeros(p1,p1,n_pt); F2i=zeros(p2,p2,n_pt); F3i=zeros(p3,p3,n_pt); F4i=zeros(p4,p4,n_pt);  
    for j=1:n_pt
      d=d_final(j);
      f1=[1 d];
      f2=[1 d/(v2+d) -v1*d/(v2+d)^2];
      f3=[1 d/(v22+d) -v12*d/(v22+d)^2];
      g=exp((v23-d)/v33);
      f4=[1 1/(1+g) -v13*g/v33/((1+g)^2)  v13*g*(v23-d)/v33^2/((1+g)^2)];
      F1i(:,:,j) = (f1'*f1);
      F2i(:,:,j) = (f2'*f2);
      F3i(:,:,j) = (f3'*f3);
      F4i(:,:,j) = (f4'*f4);
    end
    
    
    B1=zeros(p1,p1); B2=zeros(p2,p2); B3=zeros(p3,p3); B4=zeros(p4,p4);
    for j=1:n_pt
      B1 = B1 + F1i(:,:,j)*w_final(j);
      B2 = B2 + F2i(:,:,j)*w_final(j);
      B3 = B3 + F3i(:,:,j)*w_final(j);
      B4 = B4 + F4i(:,:,j)*w_final(j);
    end  

    if criterion == "D"
        Loss1d=1/(det(B1))^(1/p1);
        Loss2d=1/(det(B2))^(1/p2);
        Loss3d=1/(det(B3))^(1/p3);
        Loss4d=1/(det(B4))^(1/p4);
    elseif criterion == "A"
        Loss1d = trace(inv(B1));
        Loss2d = trace(inv(B2));
        Loss3d = trace(inv(B3));
        Loss4d = trace(inv(B4));
    else
        fprintf('Does not run.');
    end
    eff1=Loss1/Loss1d;  %Efficiency at the maximin design
    eff2=Loss2/Loss2d;
    eff3=Loss3/Loss3d;
    eff4=Loss4/Loss4d;
    eff_ex = [eff1, eff2, eff3,eff4];
    n_total;
    eff_app;
    Li = min([eff1, eff2, eff3,eff4]);
  
    LOSS(ell,:) = [ell, Li, sum(w_final), n_total];
  end
  LOSS_filter = LOSS(round(LOSS(:,3), 3) == 1 & LOSS(:,4) == n_i, :);
  [M, I] = max(LOSS_filter(:,2));
  my_loss(1, pig) = M;
  my_loss(2, pig) = n_i;
  my_loss(3, pig) = LOSS_filter(I,3);
  my_loss(4, pig) = LOSS_filter(I,4);
  my_loss(5, pig) = LOSS_filter(I,1);
end
my_table = array2table(my_loss, ...
  'RowNames', {'loss', 'n', 'sum of weight', 'n_design', 'seed'});
my_table


my_loss(1,:)/L00
