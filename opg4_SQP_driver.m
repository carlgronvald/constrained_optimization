close all
%%
%{
The function SQPSolver is used through out this driver. It has an option
argument which functions as follows:
    Options:
    log - Do you want to log the process
        values: true/false

    method - Which SQP method
        values: {'SQP'        (Plain vanilla SQP)
                ,'SQP_ls'     (Line search SQP)
                ,'SQP_trust'  (Trust region SQP)}

    infesibility_handling - For infesibility handling SQP and SQP_ls 
        values: true/false

    precision - spans from 10^-1 to 10^-9
        values: Integer between 1 and 9

    subsolver - For SQP and SQP_ls without infesibility handling there is 
                a self made interior point algorithm. 
        values: 'own solver' or 'quadprog'
%}

%% Good starting points and constraints for the himmelblau problem

%x0 = [3.584428;-1.848126];
%x0 = [0;3.8416];
%x0 = [-3.779310;-3.283186];
%x0 = [-2.805118; 3.131312]; %Svær for augmented
%x0 = [-0.270845;-0.923039];
%x0 = [-5;-4];
%x0 = [-0.5498;-0.2193];
%x0 = [0.2;3.1];
%x0 = [-200;190000];
%x0 = [15000000;4];
objective = @p2obj;
constraints = @p2cons;
plotter = @problem2Plot;
tracer = @problem2Trace;
[f,df] = objective(x0);
[c,dc] = constraints(x0);

l = [-5;-5];
u = [5;5];
cl = [-10;-10];
cu = [47;70];
x0 = [0.0;0.0];

 
%% SQP BFGS

% Comparison of own solver and the quadprog
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP', 'subsolver', 'own solver');
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
time_own = Hist.timePerformence;
plotter(x0,x,xHist, 'SQP', 'own solver')


options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP', 'subsolver', 'quadprog');
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
time_quad = Hist.timePerformence;
tracer(x0,x,xHist, 'r', 'quadprog')

legend show

figure
hold on
h1 = plot(time_own,'r');
h2 = plot(time_quad,'b');
legend([h1,h2],{'Own solver', 'quadprog'})
title('Comparison of underlying subsolvers')
ylabel('time[s]')
xlabel('Iterations')

% Several traces
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP', 'subsolver', 'own solver');
x0 = [3.584428;-1.848126];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
plotter(x0,x,xHist, 'SQP', 'Trace 1')


x0 = [0;3.8416];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
tracer(x0,x,xHist, 'r', 'Trace 2')


x0 = [-4;4];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
tracer(x0,x,xHist, 'k', 'Trace 3')

legend show

%% SQP BFGS Line Search
% Several traces
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP_ls', 'subsolver', 'own solver');
x0 = [3.584428;-1.848126];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
plotter(x0,x,xHist, 'SQP linesearch', 'Trace 1')


x0 = [0;3.8416];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
tracer(x0,x,xHist, 'r', 'Trace 2')


x0 = [-4;4];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
tracer(x0,x,xHist, 'k', 'Trace 3')

legend show

%% SQP BFGS Line Search with infeasiblity handling
% When infeasiblity handling is turned on, quadprog is the only solver for
%  the sub problem which is avalible.

x0 = [-9;-3];
try
    options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP_ls');
    [x,z, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
catch
    disp('We see that the subproblem is infeasible due to an infeasible linearization');
end
options = struct('log',true, 'infesibility_handling', true, 'method', 'SQP_ls');
[x,z, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
plotter(x0,x,xHist, 'SQP linesearch with infeasiblity handling', 'Trace 1',9)
legend show

%% SQP BFGS Trust Region
% When one uses the SQP based on a trust region, quadprog is the only solver for
%  the sub problem which is avalible.
options = struct('log',true, 'method', 'SQP_trust', 'trust_region', 0.5);
x0 = [0;0];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
plotter(x0,x,xHist, 'SQP trust region', 'Trace 1')


x0 = [0;3.8416];
[x,~, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
tracer(x0,x,xHist, 'r', 'Trace 2')


x0 = [-4;4];
[x,z, Hist] = SQPSolver(x0,objective,constraints,l,u,cl,cu,options);
xHist = Hist.xHist;
tracer(x0,x,xHist, 'k', 'Trace 3')

legend show

%% objective, constraints and derivitives
function [f,df] = obj(x)

x1 = x(1);
x2 = x(2);

temp1 = x1^2+x2-11;
temp2 = x1+x2^2-7;

f = temp1^2+temp2^2;

df = zeros(2,1);
df(1) = 4*x1*temp1+2*temp2;
df(2) = 2*temp1+4*x2*temp2;

end

function [c,dc] = cons(x,full,d)

if nargin<2
    full = false;
    d = 0;
end

if nargin<3 && full
    error('One also needs to input d to det the full dc')
end  


x1 = x(1);
x2 = x(2);

c = zeros(2,1);
c(1) = (x1+2)^2-x2;
c(2) = -4*x1+10*x2;




dc = zeros(2,2);
dc(1,1) = 2*x1+4;
dc(1,2) = -1;
dc(2,1) = -4;
dc(2,2) = 10;

dc = dc';

if full
    c = [x; -x; c; -c];     
    c = c-d;

    dc = [eye(2) -eye(2) dc -dc];
end

end


function [] = himmelTrace(x0,x,xHist, tracecolor, legendname)
hold on 
h0 = plot(x0(1,:),x0(2,:),'k.','markersize',40);
h1 = plot(x(1,:),x(2,:),'r.','markersize',40);
plot(xHist(1,:),xHist(2,:),[tracecolor,'-.'],'linewidth',2,'DisplayName',legendname);
set( get( get( h0, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
%legend(h2, legendname)
hold off
end


function [] = himmelPlot(x0,x,xHist, titlename,legendname, edges)
if nargin <6
    edges = 5;
end
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours(edges)
h1 = plot(x0(1,:),x0(2,:),'k.','markersize',40);
h2 = plot(x(1,:),x(2,:),'r.','markersize',40);
plot(xHist(1,:),xHist(2,:),'b-.','linewidth',2,'DisplayName',legendname);
%legend(h1, legendname)
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% show constraints
x1lim = edges; x1 = -x1lim:0.01:x1lim;
x2lim = edges; x2 = -x1lim:0.01:x2lim;
x1f = -edges:0.01:edges;
h3 = fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
h4 = patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

h5 = patch('Faces',[1 2 3 4], 'Vertices', [-edges -5; edges -5; edges -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h6 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; -5 edges; -5 -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h7 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; edges edges; edges 5; -edges 5], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h8 = patch('Faces',[1 2 3 4], 'Vertices', [5 edges; edges edges; edges -edges; 5 -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h5, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h6, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h7, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h8, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% set axis as needed
axis([-edges edges -edges edges])
title(titlename)
hold off
end



function [] = problem2Trace(x0,x,xHist, tracecolor, legendname)
hold on 
h0 = plot(x0(1,:),x0(2,:),'k.','markersize',40);
h1 = plot(x(1,:),x(2,:),'r.','markersize',40);
plot(xHist(1,:),xHist(2,:),[tracecolor,'-.'],'linewidth',2,'DisplayName',legendname);
set( get( get( h0, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
%legend(h2, legendname)
hold off
end


function [] = problem2Plot(x0,x,xHist, titlename,legendname, edges)
if nargin <6
    edges = 5;
end
fig = figure('Position', [100, 0, 1000,1000]);
hold on
problem2Contours(edges)
h1 = plot(x0(1,:),x0(2,:),'k.','markersize',40);
h2 = plot(x(1,:),x(2,:),'r.','markersize',40);
plot(xHist(1,:),xHist(2,:),'b-.','linewidth',2,'DisplayName',legendname);
%legend(h1, legendname)
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% show constraints
x1lim = edges; x1 = -x1lim:0.01:x1lim;
x2lim = edges; x2 = -x1lim:0.01:x2lim;
x1f = -edges:0.01:edges;
h3 = fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
h4 = patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

h5 = patch('Faces',[1 2 3 4], 'Vertices', [-edges -5; edges -5; edges -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h6 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; -5 edges; -5 -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h7 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; edges edges; edges 5; -edges 5], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h8 = patch('Faces',[1 2 3 4], 'Vertices', [5 edges; edges edges; edges -edges; 5 -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h5, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h6, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h7, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h8, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% set axis as needed
axis([-edges edges -edges edges])
title(titlename)
hold off
end

function [f,df] = p2obj(x)
    x1 = x(1);
    x2 = x(2);
    f = (x1-2)^2+x2^2+sin(x1/2);
    %f = 3*x(1)+5*x(2);
    df = zeros(2,1);
    %df(1) = 3;
    %df(2) = 5;
    df(1) = 2*x1-4+1/2*cos(x1);
    df(2) = 2*x2;
end

function [c, dc] = p2cons(x)
    c = zeros(2,1);
    c(1) = -5*x(1)*x(2);
    c(2) = -6*x(2);
    
    dc = zeros(2,2);
    dc(1,1) = 5*x(2);
    dc(1,2) = 5*x(1);
    dc(2,1) = 0;
    dc(2,2) = 6;
    dc = -dc;
end