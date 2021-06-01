% RUN SECTION BY SECTION : variable names might repeat
close all
%% Given Variables: Asset returns, target portfolio return, assets covariance
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
R = 12; 
covariance = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  

%% Setting up variables according to quadprog() documentation
H = covariance;
f = []; %f is not used, hence [] the empty brace

% This constraint allows us to specify our target return
A1 = returns;
b1 = R;

% This constraints makes portfolio weights (x) equal to 1
A2 = [1,1,1,1,1];
b2 = 1;

% Combining the above 2 equality constraints into 1; to make quadprog happy 
Aeq = [A1; A2];
beq = [b1; b2];

% Inequality constraint to disallow short-selling
Aineq = -eye(5);
bineq = zeros(5,1); 


%% Solving
%optimal port. weights found by quadprog
x = quadprog( H, f, Aineq, bineq, Aeq, beq); 

%Calculating portfolio variance (i.e., inherent risk of portfolio)
port_risk = x'*covariance*x;

%% Computing the efficient frontier

Rs = min(returns):0.01:max(returns);
port_risks = zeros(length(Rs),1);
options = optimset('Display', 'off');

for i = 1:length(Rs)
beq = [Rs(i); b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq,[],[],[],options); 
port_risks(i) = x'*covariance*x;
end

%% Plotting the efficient frontier
figure;
hold on
opt_return1 = Rs(port_risks==min(port_risks));
opt_risk1 = min(port_risks);
plot(Rs,port_risks,'b', returns, diag(covariance),'ro')
plot(opt_return1,opt_risk1 ,'ko','MarkerFaceColor', 'g')
title('The Efficient Frontier')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off


%% Bi-criterion Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%
H = covariance;
f = -returns; %f is not used, hence [] the empty brace

% This constraints makes portfolio weights (x) equal to 1
Aeq = [1,1,1,1,1];
beq = 1;

% Inequality constraint to disallow short-selling
Aineq = -eye(5);
bineq = zeros(5,1); 


%% Solving for different alpha with and without short-selling
trials = 100;
x_short = zeros(trials,5);
x_nonshort = zeros(trials,5);
x_shortOwn = zeros(trials,5);
x_nonshortOwn = zeros(trials,5);
alphas = 1:1:trials;
alphas = alphas./(trials+1);
port_risk_short = zeros(trials,1);
port_risk_nonshort = zeros(trials,1);
port_return_short = zeros(trials,1);
port_return_nonshort = zeros(trials,1);

x0 = zeros(5,1);
s0 = ones(2*5,1);
y0 = ones(length(beq),1);
z0 = ones(2*5,1);

%optimal port. weights found by quadprog
for i = 1:trials
    x_nonshort(i,:) = quadprog( alphas(i).*H, (1-alphas(i)).*f', Aineq, bineq, Aeq, beq,[],[],[],options); 
    x_nonshortOwn(i,:) = primalDualInteriorMethod_box(alphas(i).*H, (1-alphas(i)).*f',Aeq',beq,zeros(5,1),ones(5,1),x0,y0,z0,s0);
    x_short(i,:) = quadprog( alphas(i).*H, (1-alphas(i)).*f', [], [], Aeq, beq,[],[],[],options); 
    x_shortOwn(i,:) = EqualityQPSolver(alphas(i).*H, (1-alphas(i)).*f',Aeq', beq, "rangespace");
    
    port_risk_short(i,:) = x_short(i,:)*covariance*x_short(i,:)';
    port_risk_nonshort(i,:) = x_nonshort(i,:)*covariance*x_nonshort(i,:)';
    
    port_return_short(i,:) = -f*x_short(i,:)';
    port_return_nonshort(i,:) = -f*x_nonshort(i,:)';
end



%% Plot the portfolios
figure
hold on
plot(port_return_short,port_risk_short,'b', returns, diag(covariance),'ro')
title('The Bi-Criterion when shorting is allowed')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off

figure
hold on
plot(port_return_nonshort,port_risk_nonshort,'b', returns, diag(covariance),'ro')
title('The Bi-Criterion when shorting is not allowed')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off

%% Introducing a risk free security with return 0 %%%%%%%%%%%%%%%%%%%%%%%%%%
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
covariance = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  
          
returns = [returns, 0];

covariance = [covariance, zeros(5,1)];
covariance = [covariance; zeros(1,6)];

H = covariance;
f = []; %f is not used, hence [] the empty brace

% This constraint allows us to specify our target return
A1 = returns;
b1 = min(A1);

% This constraints makes portfolio weights (x) equal to 1
A2 = [1,1,1,1,1,1];
b2 = 1;

% Combining the above 2 equality constraints into 1; to make quadprog happy 
Aeq = [A1; A2];
beq = [b1; b2];

% Inequality constraint to disallow short-selling
Aineq = -eye(6);
bineq = zeros(6,1); 
%% Computing the efficient frontier with a risk-free asset 

Rs_free = min(returns):0.01:max(returns);
port_risks_free = zeros(length(Rs_free),1);

for i = 1:length(Rs_free)
beq = [Rs_free(i); b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq,[],[],[],options); 
port_risks_free(i) = x'*covariance*x;
end

%% Plotting the efficient frontier with a risk-free asset 
figure
hold on
plot(Rs_free,port_risks_free,'b', returns, diag(covariance),'ro')
title('The Efficient Frontier')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off

%% Finding the optimum for R=14

beq = [14; b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq); 
opt_risk2 = x'*covariance*x;


%% Ploting
figure
hold on
opt_return2 = 14;
plot(Rs_free,port_risks_free,'b', returns, diag(covariance),'mo')
plot(Rs,port_risks,'r')
p1 = plot(opt_return2, opt_risk2,'ko','MarkerFaceColor', 'g');
p2 = plot(opt_return1, opt_risk1,'ko','MarkerFaceColor', 'k');
title('Comparison of Efficient Frontiers')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([p1,p2],{'E[R]=14 with risk free asset','Optimal Var[R] from ex1'},'Location','northwest') 
hold off



