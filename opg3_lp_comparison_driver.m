%{
This is the driver for exercise 3. 
This file contains:
    - A test for efficiency and correctness under a growing problem using
      a random LP generator with bounds. The size of every increment is 
      controlled by n_large and the number of iterations is controlled
      by iter_test.
%}
%%
%Configure how many tests and how large an LP to test on
iter_test = 20;
n_large = 10;
sizes = n_large:n_large:(iter_test*n_large);
times = zeros(iter_test, 4);
iterations = zeros(iter_test, 4);
solution = zeros(iter_test,2);

%Plot settings
plotIterations = false;
plotTimes = true;
plotLogTimes = false;
plotSolution = false;
plotcvx = false;
plotsimplex = true;

if plotcvx
    plotsimplex = true;
end


for i = 1:iter_test

    %Create a quadrastart = cputime; program and solve it using cvx.
    n = n_large*i;
    [H, g, A, b] = generateRandomEQP(n,n/2);
    l = zeros(n,1);
    u = ones(n,1);

    
    
    if plotcvx
        start = cputime;
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize(g'*x)
            subject to
                A' * x == b
                l <= x
                x <= u
        cvx_end

        times(i,1) = cputime-start;
        iterations(i,1) = cvx_slvitr;
    end
    
    start = cputime;
    options = optimset('Display', 'off');
    options = optimset(options, 'Algorithm', 'dual-simplex');
    [x2, optval, exitflag,output] = linprog(g, [],[], A', b,l,u, options);
    times(i,2) = cputime-start;
    iterations(i,2) = output.iterations;
    
    start = cputime;
    options = optimset('Display', 'off');
    options = optimset(options, 'Algorithm', 'interior-point');
    [x3, optval, exitflag,output] = linprog(g, [],[], A', b,l,u, options);
    times(i,3) = cputime-start;
    iterations(i,3) = output.iterations;
    
    
    x0 = zeros(n,1);
    s0 = ones(2*n,1);
    y0 = ones(length(b),1);
    z0 = ones(2*n,1);
    
    start = cputime;
    [x,y,z,s, iter] = LinearPDIM_box(g,A,b,l,u,x0,y0,z0,s0);
    times(i, 4) = cputime-start;
    iterations(i, 4) = iter;
    
    disp("Iteration " +i+"/"+iter_test);
    disp("Mean error wrt simplex: " + mean(sqrt((x-x2).^2)))
    disp("Mean error wrt interior: " + mean(sqrt((x-x3).^2)))
    
    solution(i,1) = mean(sqrt((x-x2).^2));
    solution(i,2) = mean(sqrt((x-x3).^2));
    
    if exitflag ~= 1
        disp("Iteration "+i+" is infeasible!!")
    end
end


% PLOT ITERATIONS
if plotIterations
    disp("Plotting # of iterations!")
    figure;
    if(plotcvx)
        plot(sizes, iterations(:,1))
        hold on
    end
    if(plotsimplex)
        plot(sizes, iterations(:,2))
        hold on
    end
    plot(sizes, iterations(:,3))
    hold on
    plot(sizes, iterations(:,4))
    hold off
    if(plotcvx)
        legend(["cvx iterations", "linprog dual-simplex iterations", "linprog interior-point iterations", "our solver"])
    elseif(plotsimplex)
        legend(["linprog dual-simplex iterations", "linprog interior-point iterations", "our solver"])
    else
        legend([ "linprog interior-point iterations", "our solver"])
    end
    ylabel("iterations")
    %ylim([0 max(max(iterations))+2])
    xlabel("n")
    title("iterations vs n")

end

if plotTimes
    disp("Plotting times!")
    figure;
    %Note: Time plotted is total wall time between starting and ending
    %solver
    % Not CPU time (since CPU time counts the number of cores used)
    if(plotcvx)
        plot(sizes, times(:,1))
        hold on
    end
    plot(sizes, times(:,2))
    hold on
    plot(sizes, times(:,3))
    plot(sizes, times(:,4))
    hold off
    if(plotcvx)
        legend(["cvx time", "linprog dual-simplex time", "linprog interior-point time", "our solver"])
    else
        legend(["linprog dual-simplex time", "linprog interior-point time", "our solver"])
    end
    ylabel("t [s]")
    xlabel("n")
    title("time vs n")
end

if plotLogTimes
    disp("Plotting times!")
    clear log
    figure;
    %Note: Time plotted is total wall time between starting and ending
    %solver
    % Not CPU time (since CPU time counts the number of cores used)
    if(plotcvx)
        loglog(sizes, times(:,1))
        hold on
    end
    loglog(sizes, times(:,2))
    hold on
    loglog(sizes, times(:,3))
    loglog(sizes, times(:,4))
    hold off
    if(plotcvx)
        legend(["cvx time", "linprog dual-simplex time", "linprog interior-point time", "our solver"])
    else
        legend(["linprog dual-simplex time", "linprog interior-point time", "our solver"])
    end
    ylabel("t [s]")
    xlabel("n")
    title("time vs n")
end

if plotSolution
    disp("Plotting solution!")
    figure;
    plot(sizes, solution)
    title("Correctness")
    ylabel("Mean Squared Error")
    xlabel("n")
    legend(["Relative to Simplex","Relative to interior point"])
end