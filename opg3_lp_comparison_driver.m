%Configure how many tests and how large
largeProb = true;
n_large = 10;
iter_test = 20;
sizes = n_large:n_large:(iter_test*n_large);
times = zeros(iter_test, 4);
iterations = zeros(iter_test, 4);
i = 1;
solution = zeros(iter_test,2);
for i = 1:iter_test

    %Create a quadratic program and solve it using cvx.
    n = n_large*i;
    [H, g, A, b] = generateRandomEQP(n,n/2);
    l = zeros(n,1);
    u = ones(n,1);

    
    tic

    cvx_begin quiet
        %cvx_precision low
       variable x(n)
        minimize(g'*x)
        subject to
            A' * x == b
            l <= x
            x <= u
    cvx_end

    times(i,1) = toc;
    iterations(i,1) = cvx_slvitr;
    
    tic
    options = optimset('Display', 'off');
    options = optimset(options, 'Algorithm', 'dual-simplex');
    [x2, optval, exitflag,output] = linprog(g, [],[], A', b,l,u, options);
    times(i,2) = toc;
    iterations(i,2) = output.iterations;
    
    tic
    options = optimset('Display', 'off');
    options = optimset(options, 'Algorithm', 'interior-point');
    [x3, optval, exitflag,output] = linprog(g, [],[], A', b,l,u, options);
    times(i,3) = toc;
    iterations(i,3) = output.iterations;
    
    
    x0 = zeros(n,1);
    s0 = ones(2*n,1);
    y0 = ones(length(b),1);
    z0 = ones(2*n,1);
    
    tic
    [x,y,z,s, iter] = LinearPDIM_box_normal_equation(g,A,b,l,u,x0,y0,z0,s0);
    times(i, 4) = toc;
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

plotIterations = true;
plotTimes = true;
plotSolution = true;
plotcvx = false;

% PLOT ITERATIONS
if plotIterations
    disp("Plotting # of iterations!")
    figure;
    if(plotcvx)
        plot(sizes, iterations(:,1))
        hold on
    end
    plot(sizes, iterations(:,2))
    hold on
    plot(sizes, iterations(:,3))
    plot(sizes, iterations(:,4))
    hold off
    if(plotcvx)
        legend(["cvx iterations", "linprog dual-simplex iterations", "linprog interior-point iterations", "our solver"])
    else
        legend(["linprog dual-simplex iterations", "linprog interior-point iterations", "our solver"])
    end
    ylabel("iterations")
    ylim([0 max(max(iterations))+2])
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
        legend(["cvx iterations", "linprog dual-simplex iterations", "linprog interior-point iterations", "our solver"])
    else
        legend(["linprog dual-simplex iterations", "linprog interior-point iterations", "our solver"])
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