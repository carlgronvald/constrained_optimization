function [x,z,i] = intSQP_inf(B,df,dc,l,u,gl,gu,x0)
%{
This solver takes programs of the form 

            min_x f(x)
            s.t   gu>= cx >= gl
                   u>=  x >= l

The solver is based on Mehrotas predictor-corrector primal-dual algorithm. 

%}
    n = length(x0);
    m = length(gu);
    mc = n*2+2*m;
    epsilon = 0.0001;
    max_iter = 100;
    eta = 0.995;
    
    mu = 100;
    mu_vec = mu.*ones(2*m,1);
    
    tl = zeros(m,1);
    tu = zeros(m,1);
    x = [x0; tl; tu];
    
    d = [l;-u;gl-tl;-gu-tu; zeros(2*m,1)];
    
    mc = mc + m*2;
    z = ones(mc,1);
    s = ones(mc,1);
    
    e = ones(mc,1);
    
    while(any(s==0))
        x = x+0.000001;
        s = [x;-x; dc'*x; -dc'*x]-d;
    end
    
    tl = x(n+1:n+m);
    tu = x(n+m+1:n+2*m);
    xx =x(1:n);
    
    
    
    %Initilize constraint specific slacks and lagrange multipliers
    sl = s(1:n);
    su = s(n+1:n*2);
    scl = s(2*n+1:2*n+m);
    scu = s(m+2*n+1:n*2+2*m);
    stl = s(n*2+2*m+1:n*2+3*m);
    stu = s(n*2+3*m+1:n*2+4*m);
    
    zl = z(1:n);
    zu = z(n+1:n*2);
    zcl = z(2*n+1:2*n+m);
    zcu = z(m+2*n+1:n*2+2*m);
    ztl = z(n*2+2*m+1:n*2+3*m);
    ztu = z(n*2+3*m+1:n*2+4*m);
    
    %initial residuals
    rLx = B*xx+df-(zl-zu+dc*zcl-dc*zcu);
    rLt = mu_vec-[ztl;ztu];
    rL = [rLx;rLt];
    
    rC = s+d - [xx; -xx; dc'*xx; -dc'*xx; tl; tu];
    
    
    % Initial point heuristic
    zsl = diag(zl./sl);
    zsu = diag(zu./su);
    zslc = zcl./scl;
    zsuc = zcu./scu;
    zc = zslc + zsuc;
    Hbar = B + zsl + zsu + bsxfun(@times,zc',dc)*dc';
    
    zstl = ztl./stl;
    zstu = ztu./stu;
    zst = [zstl; zstu];
    
    zer = zeros(n,m*2);
    Hs = [Hbar zer; zer' diag(zst)];
    Hsp = sparse(Hs);
    
    [L,D, P] = ldl(Hsp);


    % Affine step
    rCs = (rC-s);
    
    rLx = [ zsl -zsu bsxfun(@times,zslc',dc) bsxfun(@times,zsuc',-dc)]*rCs(1:2*n+2*m);
    rLt = zst.*(rCs(1+2*n+2*m:end));
    rLbar = rL - [rLx;rLt];
    
    rhs = -rLbar;
    dxAff = P*(L' \ (D \ (L \ (P'*rhs) )));
    
    dxAffx = - [ zsl; -zsu; bsxfun(@times,zslc,dc'); bsxfun(@times,zsuc,-dc')]*dxAff(1:n);
    dxAfft = -zst.*dxAff(n+1:end);

    dzAff = [dxAffx;dxAfft] + (z./s).*rCs;
    dsAff = -s-(s./z).*dzAff;
  
    %Update of starting point
    z = max(1,abs(z+dzAff));
    s = max(1,abs(s+dsAff));
    
    %Update of initial residuals
    sl = s(1:n);
    su = s(n+1:n*2);
    scl = s(2*n+1:2*n+m);
    scu = s(m+2*n+1:n*2+2*m);
    stl = s(n*2+2*m+1:n*2+3*m);
    stu = s(n*2+3*m+1:n*2+4*m);
    
    zl = z(1:n);
    zu = z(n+1:n*2);
    zcl = z(2*n+1:2*n+m);
    zcu = z(m+2*n+1:n*2+2*m);
    ztl = z(n*2+2*m+1:n*2+3*m);
    ztu = z(n*2+3*m+1:n*2+4*m);

    
    rLx = B*xx+df-(zl-zu+dc*zcl-dc*zcu);
    rLt = mu_vec-[ztl;ztu];
    rL = [rLx;rLt];
    rC = s+d - [xx; -xx; dc'*xx; -dc'*xx; tl; tu];
    rSZ = s.*z;
    
    dualGap = (z'*s)/(mc);
    dualGap0 = dualGap;
    
    for i = 1:max_iter
        zsl = diag(zl./sl);
        zsu = diag(zu./su);
        zslc = zcl./scl;
        zsuc = zcu./scu;
        zc = zslc + zsuc;
        
        % Factorization
        Hbar = B + zsl + zsu + bsxfun(@times,zc',dc)*dc';
    
        zstl = ztl./stl;
        zstu = ztu./stu;
        zst = [zstl; zstu];

        zer = zeros(n,m*2);
        Hs = [Hbar zer; zer' diag(zst)];
        Hsp = sparse(Hs);

        [L,D, P] = ldl(Hsp);

        % Affine step
        rCs = (rC-s);

        rLx = [ zsl -zsu bsxfun(@times,zslc',dc) bsxfun(@times,zsuc',-dc)]*rCs(1:2*n+2*m);
        rLt = zst.*(rCs(1+2*n+2*m:end));
        rLbar = rL - [rLx;rLt];

        rhs = -rLbar;
        dxAff = P*(L' \ (D \ (L \ (P'*rhs) )));

        dxAffx = - [ zsl; -zsu; bsxfun(@times,zslc,dc'); bsxfun(@times,zsuc,-dc')]*dxAff(1:n);
        dxAfft = -zst.*dxAff(n+1:end);

        dzAff = [dxAffx;dxAfft] + (z./s).*rCs;
        dsAff = -s-(s./z).*dzAff;
        
        % Compute max alpha affine
        dZS = [dzAff; dsAff];
        alphas = (-[z;s]./dZS);
        alphaAff = min([1;alphas(dZS<0)]);
           
        dualGapAff = ((z+alphaAff*dzAff)'*(s+alphaAff*dsAff))/(mc);
        sigma = (dualGapAff/dualGap)^3;
        
        
        % Affine-Centering-Correction Direction
        rSZbar = rSZ + dsAff.*dzAff-sigma*dualGap*sigma*e;
        
        
        rCs = (rC-rSZbar./z);

        rLx = [ zsl -zsu bsxfun(@times,zslc',dc) bsxfun(@times,zsuc',-dc)]*rCs(1:2*n+2*m);
        rLt = zst.*(rCs(1+2*n+2*m:end));
        rLbar = rL - [rLx;rLt];

        rhs = -rLbar;
        dx = P*(L' \ (D \ (L \ (P'*rhs) )));
        
        
        dxx = - [ zsl; -zsu; bsxfun(@times,zslc,dc'); bsxfun(@times,zsuc,-dc')]*dx(1:n);
        dxt = -zst.*dx(n+1:end);

        dz = [dxx;dxt] + (z./s).*rCs;

        ds = -rSZbar./z-(s./z).*dz;
        
        %compute max alpha
        dZS = [dz; ds];
        alphas = (-[z;s]./dZS);
        alpha = min([1;alphas(dZS<0)]);
        
        alphaBar = eta*alpha;
        
        % Update of position 
        x = x + alphaBar * dx;
        z = z + alphaBar * dz;
        s = s + alphaBar * ds;
        
        tl = x(n+1:n+m);
        tu = x(n+m+1:n+2*m);
        xx =x(1:n);
        
        %d = [l;-u;gl-tl;-gu-tu; zeros(2*m,1)];
        
        % Update of residuals
        sl = s(1:n);
        su = s(n+1:n*2);
        scl = s(2*n+1:2*n+m);
        scu = s(m+2*n+1:n*2+2*m);
        stl = s(n*2+2*m+1:n*2+3*m);
        stu = s(n*2+3*m+1:n*2+4*m);

        zl = z(1:n);
        zu = z(n+1:n*2);
        zcl = z(2*n+1:2*n+m);
        zcu = z(m+2*n+1:n*2+2*m);
        ztl = z(n*2+2*m+1:n*2+3*m);
        ztu = z(n*2+3*m+1:n*2+4*m);


        rLx = B*xx+df-(zl-zu+dc*zcl-dc*zcu);
        rLt = mu_vec-[ztl;ztu];
        rL = [rLx;rLt];

        rC = s+d - [xx; -xx; dc'*xx; -dc'*xx; tl; tu];
        rSZ = s.*z;

        dualGap = (z'*s)/(mc);

        if(dualGap <= epsilon*0.01*dualGap0)
           return
        end
    end  
end

