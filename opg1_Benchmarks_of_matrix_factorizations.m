tests = 20;
times = ones(5,tests)*tests;
ns = zeros(tests,1);
for i = 1:tests
    n = i*(1000/tests);
    ns(i) = n;

    [H,g,A,b,x,lambda] = generateRandomEQP(n,n);
    KKT = [H -A;-A', zeros(size(A,2), size(A,2))];
    for k =1:5
        j=1;
        
        tic
        x = lu(KKT);
        times(j, i) = min(toc, times(j,i));
        j = j+1;

        tic
        x = lu(H);
        times(j, i) = min(toc, times(j,i));
        j = j+1;

        tic
        x = ldl(KKT);
        times(j, i) = min(toc, times(j,i));
        j = j+1;

        tic
        x = ldl(H);
        times(j, i) = min(toc, times(j,i));
        j = j+1;

        tic
        x = chol(H);
        times(j, i) = min(toc, times(j,i));
        j = j+1;


    end

    disp("n="+n)
end


hold off
for i=1:5%size(times,1)
    plot(ns, times(i,:))
    hold on
end
legend(["LU, KKT", "LU, H", "LDL, KKT", "LDL, H", "Cholesky, H"])
xlabel("n")
ylabel("time [s]")
title("Growing size problem")



    
