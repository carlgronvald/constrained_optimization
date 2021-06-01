A = [2 3 5; 4 5 5; 5 6 5];
d = [2;3; 5];

% A*D
bsxfun(@times,d',A)
A*diag(d)

% D*A
bsxfun(@times,d,A')
diag(d)*A'

% A'*D*A
bsxfun(@times,d',A)*A'
A*diag(d)*A'




%%
M = randn(10000,10);
D = diag(randn(10000,1).^2);


tic
A = M'*D*M;
toc


tic
B = bsxfun(@times,M,sqrt(diag(D)));
B = B.'*B;
toc

d = diag(D);

tic
C = bsxfun(@times,d',M')*M;
toc

