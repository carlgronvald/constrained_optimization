function [x,z,Hist] = SQP_trust(x0,obj,con,l,u,cl,cu,log,precision)
    max_iter = 100;
    n = length(x0);
    epsilon = 10^(-precision);
    penalty = 100;
    dk0 = 0.5;
    dk = dk0;
    time = 0;

    % Define relevant functions call
    x = x0;
    [~,df] = feval(obj,x);
    [c,dc] = feval(con,x);
    B = eye(n);
    m = size(c,1);
    z = ones(2*m+2*n,1);
    lid = 1:m;
    uid = (m+1):(2*m);
    clid = (2*m+1):(2*m+n);
    cuid = (2*m+n+1):(2*(n+m));
    
    % Log objects
    if log
        xHist = zeros(n,max_iter+1);
        pkHist = zeros(n,max_iter);
        timePerformence = zeros(1,max_iter);
        trustRegion = zeros(1,max_iter+1);
        
        xHist(:,1) = x0;
    end
    
    % Define options for quadprog
    options = optimset('Display', 'off');
    

        
    
    % Start for loop
    for i = 1:max_iter
        
        % Update lower and upper bounds for the quadratic approximation
        lk = -x+l;
        uk = -x+u;
        clk = -c+cl;
        cuk = -c+cu;
        
        j = 0;
        % Define infesibility program
        Hinf = [B zeros(n,2*m); zeros(2*m,2*m+n)];
        ginf = [df; penalty*ones(2*m,1)];
        iden = eye(m);
        Cinf = [dc -dc zeros(m,2*m) iden  -iden;iden zeros(m,m) iden zeros(m,m) zeros(m,2*m); zeros(m,m) iden zeros(m,m) iden zeros(m,2*m) ]';
        lkinf = [lk; zeros(2*m,1)];
        ukinf = [uk; inf(2*m,1)];
        while(true)
            j = j+1;
   
            dinf = [clk; -cuk; zeros(2*m,1); -dk*ones(2*m,1)];
            % Solves local QP program
            tic
            [pk,~,~,~,zhat] = quadprog(Hinf,ginf,-Cinf,-dinf,[],[],lkinf,ukinf,[],options);
            time = time + toc;
            if ~isempty(pk)
                if log
                   trustRegion(:,i) = dk; 
                end
                dk = dk0;
                break
            elseif j>30
                msg = ['The trust region is of size ',num2str(dk*2^(j-1),'%02d'),' and the step is still not feasible. Try the method SQP_ls'];
                error(msg)
            else
                dk = dk*2;
            end
        end
        

        zhat = [zhat.lower(1:n); zhat.upper(1:n); zhat.ineqlin(1:2*m)];
        pk = pk(1:n);
        pz = zhat-z;
        
        

        % Update the current point
        z = z + pz;
        x = x + pk;
        

        % For the quasi Newton update  
        dL = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

        % Update values for next iteration
        [~,df] = feval(obj,x);
        [c,dc] = feval(con,x);

        % Quasi newton update of the hessian
        dL2 = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

        p = pk;
        q = dL2-dL;
        theta = 1;
        Bp = (B*p);
        pBp = p'*Bp;
        
        if p'*q < 0.2*pBp
            theta = (0.8*pBp)/(pBp-p'*q);
        end
        r = theta*q+(1-theta)*(Bp);
        B = B + r*r'/(p'*r) - Bp*Bp'/pBp;
        
        if log
            pkHist(:,i) = pk;
            xHist(:,i+1) = x;
            timePerformence(1,i) = time;
            time = 0;
        end
        
        
        if norm(dL2, 'inf')<epsilon
            if log
                pkHist = pkHist(:,1:i);
                xHist = xHist(:,1:i+1);
                timePerformence = timePerformence(:,1:i);
                trustRegion = trustRegion(:,1:i);
                
                Hist = struct('xHist', xHist, 'pkHist', pkHist, 'timePerformence', timePerformence, 'Iterations' ,i, 'trustRegion', trustRegion);
            else
                Hist = struct();
            end
            
            return
        end
   end

end

