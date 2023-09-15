%%Compute D-optimal designs for polynomial regression 
%%with exact n points
 

clear;
runningtime=cputime;  %record computation time
    tol = 10^(-4);
    N=21;         %number of design points   
     
    a=  -1;   %[a, b] is the design space
    b=   1;   %
     
     p=3;            %degree of polynomial regression model 
     
    u=linspace(a,b,N); %equally spaced N points in [a,b]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               MODEL SET UP                       %        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    
    %The following vectors and matrices are used in the information
    %matrices below.
    f = power(u(:), 0:p);          
    
    Fi=zeros(p+1,p+1,N);
    
    for j=1:N 
        f1=[f(j,:) ];
        Fi(:,:,j)=f1'*f1;
    end
   
%Compute the D-optimal design   
        cvx_begin
            cvx_precision high
            variable w(1,N);
            expression A(p+1,p+1); 
            
            for j=1:N
                A = A+Fi(:,:,j)*w(j);
            end
              
            %minimize( trace_inv(A) )   %A-opt
            minimize (-log(det_rootn(A)))
            0 <= w <= 1;
            sum(w)==1;
        cvx_end
        
       

    design=[u(find(w(1,:)>tol))' w(1,find(w(1,:)>tol))'] %optimal design 
    d00=u(find(w(1,:)>tol)) %support points
    w00=w(w(1,:) > tol) %optimal weight
    L00=cvx_optval

%Find n exact design points using an annealing algorithm with the
%following setting
n=19;
c0=1; %max number of points to be changed in the annealing algorithm
Nt=200; %number of iterations per temperature change
T0=0.1; %initial temperature
M0=500; %number of temperature changes before algorithm stops
alpha=0.9;  %temperature cooling rate

Tmin=T0*alpha^M0; %minimum temperature

delta=2*(b-a)/(N-1); %neighbourhood size

w00 = initializeExact(w00,n) %convert approximate design lazily to an exact design
d0=d00;
w0=w00;


k=length(w0); 

FIM=FIM_1dregPp(d00,3);
FIM=sum(FIM.*reshape(w0,1,1,[]),3);
q=p+1;

L0=-log(det(FIM)^(1/q));  %D-optimality

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
while(T > Tmin)

    % GENERATING RANDOM CANDIDATE DESIGN
    for h = 1:Nt
        
        k = length(w0); % current number of support points
        ci =  min([k,randi([1 c0])]); % number of points to be replaced
        toRemove = randperm(k,ci); % randomly select ci points to remove
        
        % creating new candidate design
        wi = zeros(k+ci,1); 
        di = zeros(k+ci,1);
        
        % copy previous design
        wi(1:k) = w0(:);
        di(1:k,:) = d0(:,:);
        
       
        % remove selected points
        for j = toRemove
            wi(j) = wi(j) - 1/n; 
        end
    
        % add weight for new points
        for j = 1:ci
            wi(k+j) = wi(k+j) + 1/n; 
        end
        
        % generate new points in neighbourhood of those previously removed
        i = 1;
        for j = toRemove
                di_nb =di(j,:)+ (2*rand(1)-1)*delta;
                di_nb=min([max([a, di_nb]),b]); 
                di(k+i,:) = di_nb;
                i = i + 1; 
        end
        % points are generated in a circular, spherical, or hyperspherical
        % region (depending on the dimension of the design) of radius
        % delta
        
        % remove support points with zero weight
        di = di(wi>0.0001,:);
        wi = wi(wi>0.0001);
        
        % compute loss of candidate design
        FIMi = FIM_1dregPp(di,3);
        FIMi = sum(FIMi.*reshape(wi,1,1,[]),3);
        Li = -log(det(FIMi)^(1/q));
   
   % PROCEED WITH ANNEALING STEP

    prob = exp(-(Li-L0)/T); % acceptance probability

    if prob > rand(1) % criterion for accepting random design
        L0 = Li;
        d0 = di;
        w0 = wi;
    
     num_iters = num_iters+1;
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
min(loss1)
d0
w0
sum(w0)
design=sort(d0)

%Need to work on this plot below!!!!!
figure;
scatter(d0,n*w0,"blue");
xlabel("support points");
ylabel("weights");
title("Exact design distribution")

% Plot initial and final design
%figure;
%plt1 = scatter3(d0(:,1),d0(:,2),d0(:,3),"filled", "red",'MarkerFaceAlpha', 0.5,'MarkerEdgeColor', 'black');
%title("Exact Design");
%xlabel("x1");
%ylabel("x2");
%zlabel("x3");
%hold on
%plt2 = scatter3(d(:,1),d(:,2),d(:,3),"filled", "blue",'MarkerFaceAlpha', 0.5);
%legend([plt1,plt2], ["Final Design", "Initial Design"])
%grid on

resulttime=cputime-runningtime  %computation time
 
     
[L00, loss1(1), min(loss1)]