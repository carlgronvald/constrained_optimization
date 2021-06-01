names = ["LDLdense", "LDLsparse", "LUdense", "LUsparse", "rangespace", "nullspace", "quadprog"];

% What to test
bdependent = true;
ndependent = false;

% Given problem
if bdependent
    tests = 20;
    times = zeros(7,tests);
    answers = zeros(5,tests,7);
    bs = [];
    i=1;
    for b1 = linspace(8.5, 18.68, tests)
        bs = [bs;b1];
        n = 5;
        H = [5.0 1.86 1.24 1.48 -0.46; ...
            1.86 3.0 0.44 1.12 0.52; ...
            1.24 0.44 3.8 1.56 -0.54; ...
            1.48 1.12 1.56 7.2 -1.12; ...
            -0.46 0.52 -0.54 -1.12 7.8];
        g = [-16.1; -8.5; -15.7; -10.02; -18.68];
        A = [16.1 8.5 15.7 10.02 18.68; 1.0 1.0 1.0 1.0 1.0]';
        b = [b1;1];

        tic
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLdense");
        answers(:, i,1) = x;
        times(1, i) = toc;
        tic
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLsparse");
        answers(:, i,2) = x;
        times(2, i) = toc;
        tic
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LUdense");
        answers(:, i,3) = x;
        times(3, i) = toc;
        tic
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LUsparse");
        answers(:, i,4) = x;
        times(4, i) = toc;
        tic
        [x, lambda] = EqualityQPSolver(H,g,A,b, "rangespace");
        answers(:, i,5) = x;
        times(5, i) = toc;
        tic
        [x, lambda] = EqualityQPSolver(H,g,A,b, "nullspace");
        answers(:, i,6) = x;
        times(6, i) = toc;
        options = optimset('Display', 'off');
        tic
        [x, lambda] = quadprog(H,g,[],[],A',b,[],[],[],options);
        answers(:, i,7) = x;
        times(7, i) = toc;


        i = i+1;
    end
    
    answer_diff = zeros(6,20);
    for i=1:6
        answer_diff(i,:) = mean(sqrt((answers(:,:,i)-answers(:,:,7)).^2));
    end
    
    

    hold off
    for i=1:size(times,1)
        plot(bs, times(i,:))
        hold on
    end
    legend(names)
    xlabel("b(1)")
    ylabel("time [s]")
    
    figure
    hold off
    for i=1:6
        plot(bs, (answer_diff(i,:)))
        hold on
    end
    legend(names(1:6))
    xlabel("b(1)")
    ylabel("Error relative to quadprog")
    
end

if ndependent
    tests = 20;
    times = ones(7,tests)*tests;
    ns = zeros(tests,1);
    for i = 1:tests
        n = i*(2000/tests);
        ns(i) = n;
        
        [H, g, A, b] = ProblemEQPRecycling(n, 0.2, 1);
        for k =1:5 % Take the minimum runtime of running it five times, to avoid 
            j=1;
            
            tic
            [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLdense");
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
            tic
            EqualityQPSolver(H,g,A,b, "LDLsparse");
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
            tic
            EqualityQPSolver(H,g,A,b, "LUdense");
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
            tic
            EqualityQPSolver(H,g,A,b, "LUsparse");
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
            
            tic
            EqualityQPSolver(H,g,A,b, "rangespace");
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
            tic
            EqualityQPSolver(H,g,A,b, "nullspace");
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
            options = optimset('Display', 'off');
            tic
            quadprog(H,g,[],[],A',b,[],[],[],options);
            times(j, i) = min(toc, times(j,i));
            j = j+1;
            
        end
        
        disp("n="+n)
    end


    hold off
    for i=1:size(times,1)
        plot(ns, times(i,:))
        hold on
    end
    legend(names, 'Location','northwest')
    xlabel("n")
    ylabel("time [s]")
    title("Growing size problem")
end





