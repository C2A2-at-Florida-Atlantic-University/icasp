function [Q, err, T, Z] = L1PCA_uniform_feature_preservation(X, K)
% Q is a K-dimentional suboptimal basis
% err = norm1(X - Proj_L1(Q, X)) % elementwise
% T = abs(X - Proj_L1(Q, X))
% Algorithm: Uniform Feature Preservation (UB)

szX = size(X);
[U, S] = svd(X);
s = diag(S);
% clear S
rnk = sum( s>(max(szX)*max(s)*eps) );
if rnk <= K
    Q = U(:, 1:K);
    err = 0;
    if nargout > 2
        Z = Q'*X;
        T = zeros(szX);
    end
    return
end

D = szX(1);
Is = combnk(1:D, D-K); % Caution: nchoosek( D, K ) could be huge.
dchk = size(Is, 1);
prf = zeros(dchk, 1);
Qs = zeros(D, K, dchk);
for ii = 1:dchk
    I = Is(ii, :);
    Ic = 1:D;
    Ic(I) = []; % compliment of `I`

    [Ct, ~] = solveCoreProb(X(Ic,:)', -X(I,:)');

    Qs(:,:,ii) = [eye(K); -Ct'];
    Qs([Ic I], :, ii) = Qs(:,:,ii);
    if true == false
        [Qs(:,:,ii), ~] = qr(Qs(:,:,ii), 0); % Optional
    end
    try 
    [~, prf(ii)] = solveCoreProb(Qs(:,:,ii), X);
    catch
        return;
    end
end

[err, I_prf] =  min(prf);
[Q, ~] = qr(Qs(:,:,I_prf), 0);
[Z, ~, T] = solveCoreProb(Q, X);


%% %%%% %%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%                   solveCoreProb                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%% %%%% %%

function [X, val, T, Nhopes] = solveCoreProb(A, B, orthogonalize) 

% f(X_) = norm1(A*X_-B) % elementwise L1 norm
% val =    min f
% X   = argmin f
% T   = abs(A*X-B)
% Algorithm: l1-magic (Emmanuel Cand`es and Justin Romberg, Caltech)

%%
N = size(B, 2);
[D, K]= size(A);

if ~exist('orthogonalize','var') || isempty(orthogonalize), orthogonalize = true; end % false works too
Nhopes = -1;

%%
if orthogonalize
    %% Orthogonalization and column dependencies avoidance (Begin)
    A_copy = A;
    B_copy = B;
    [Q, R, e] = qr(A, 0); % this is a RRQR on A (not on A')

    if K==1, dgR=R(1); else, dgR=diag(R); end
    p = find([abs(dgR); 0] < 1e-11, 1);
    prm   = e(1:p-1);
    prm_c = e(p:K);
    Z = A(:, prm)\A(:, prm_c); % i need this for later
    A = Q(:,1:p-1);
    R = R(1:p-1,:);

    B = B - A*(A'*B);
    non_zero_RHS_ind = find(sqrt(sum(abs(B).^2, 1)) >= 1e-11); % Non-zero columns
    B = B(:, non_zero_RHS_ind); % Drop the zero columns

    if isempty(B)
        opts.UT=true; X=linsolve(R,A'*B_copy,opts); % O(K^2) to retrieve X
        e(e) = 1:K; % invert permutation
        X = X(e,:);
        val = 0; % zero IS the minimum objective value
        T = zeros(D, N);
        return
    end
end

%%
szB2 = size(B, 2); % New `N` % These could have changed
szA2 = size(A, 2); % New `K` % during orthogonalization

if orthogonalize
    X = zeros(szA2, szB2);
else
    X = pinv(A)*B;
end

%%
Nhopes = zeros(szB2, 1);
for n = 1:szB2
    try
    [X(:,n), ~, Nhopes(n)] = l1decode_pd_matrix_only(A, B(:,n), X(:,n));
    catch
        return;
    end 
end

%%
if orthogonalize
    %% Orthogonalization and column dependencies avoidance (End)
    X_copy = zeros(szA2,N);
    X_copy(:,non_zero_RHS_ind) = X;
    X = X_copy;
    opts.UT=true; X=linsolve(R(:,1:p-1),X+A'*B_copy,opts); % plenty of solutions.. keep the ones with zeros down bellow, and chop these zeros.

    A = A_copy;
    X_ = randn(length(prm_c), N); % ones(length(prm_c),N); % zeros(length(prm_c),N); %
    X_copy = zeros(K, N);
    X_copy(prm, :) = X-Z*X_;
    X_copy(prm_c, :) = X_;
    X = X_copy;

    T = abs(A*X - B_copy);

    Nhopes_new = zeros(N, 1);
    Nhopes_new(non_zero_RHS_ind) = Nhopes;
    Nhopes = Nhopes_new; % If the RHS is in the span of A,
                         % then the number of hopes is zero.
else
    T = abs(A*X - B);
end

%%
val = sum(sum(T));



%% %%%% %%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%              l1decode_pd_matrix_only                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%% %%%%%%%% %%%% %%

% THIS FUNCTION IS AN EDITED VERSION OF THE FUNCTION
%         l1decode_pd.m
%
% Decoding via linear programming. Solve
% min_x  ||b-Ax||_1 .
%
% Recast as the linear program
% min_{x,u} sum(u)  s.t.  -Ax - u + y <= 0
%                          Ax - u - y <= 0
% and solve using primal-dual interior point method.
%
% A  - MxN matrix (M > N). 
% y  - Mx1 observed code.
% x0 - Nx1 vector, initial point.

function [xp, st, pditer] = l1decode_pd_matrix_only(A, y, x0)

st = 0;
[M, N] = size(A);

alpha=1e-2; beta=5e-1; mu=1e1;

gradf0 = [zeros(N,1); ones(M,1)];

if nargin<3, x=pinv(A)*y; else, x=x0; end

Ax = A*x;
u = .95*abs(y-Ax) + .1*max(abs(y-Ax));

fu1 =  Ax - y - u;
fu2 = -Ax + y - u;

lamu1 = -1./fu1;
lamu2 = -1./fu2;

Atv = A'*(lamu1-lamu2);

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*M/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [Atv; -lamu1-lamu2];
resnorm = norm([rdual; rcent]);

opts.POSDEF = true; 
opts.SYM = true;

pditer = 0;
done = (sdg < 1e-8) | (pditer >= 50);
while ~done    
    pditer = pditer + 1;

    w2 = -1 - 1/tau*(1./fu1 + 1./fu2);

    sig1 = -lamu1./fu1 - lamu2./fu2;
    sig2 =  lamu1./fu1 - lamu2./fu2;
    sigx = sig1 - sig2.^2./sig1;

    w1 = -1/tau*(A'*(-1./fu1 + 1./fu2));
    w1p = w1 - A'*((sig2./sig1).*w2);
    A_with_scaled_rows = repmat(sigx,1,N).*A;
%     H11p = A'*(sparse(diag(sigx))*A);
    H11p = A'*A_with_scaled_rows;
    try
    [dx, hcond] = linsolve(H11p, w1p,opts);
    catch
        return; 
    end
    if hcond < 1e-14
        disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
        st=-1;
        xp=x;
        return
    end
    Adx = A*dx;

    du = (w2 - sig2.*Adx)./sig1;

    dlamu1 = -(lamu1./fu1).*(Adx - du) - lamu1 - (1/tau)*1./fu1;
    dlamu2 =  (lamu2./fu2).*(Adx + du) - lamu2 - (1/tau)*1./fu2;
    Atdv = A'*(dlamu1-dlamu2);

    % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
    indl = find(dlamu1 < 0); 
    indu = find(dlamu2 < 0);
    s = min([1; -lamu1(indl)./dlamu1(indl); -lamu2(indu)./dlamu2(indu)]);
    indl = find(( Adx-du) > 0);  
    indu = find((-Adx-du) > 0);
    s = .99*min([s; -fu1(indl)./(Adx(indl)-du(indl)); -fu2(indu)./(-Adx(indu)-du(indu))]);

    % backtrack
    suffdec = 0;
    backiter = 0;
    while ~suffdec
        xp = x + s*dx;  
        up = u + s*du;
        Axp = Ax + s*Adx;  
        Atvp = Atv + s*Atdv;
        lamu1p = lamu1 + s*dlamu1; 
        lamu2p = lamu2 + s*dlamu2;
        fu1p = Axp - y - up; 
        fu2p = -Axp + y - up;
        rdp = gradf0 + [Atvp; -lamu1p-lamu2p];
        rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
        suffdec = (norm([rdp; rcp]) <= (1-alpha*s)*resnorm);
        s = beta*s;
        backiter = backiter + 1;
        if backiter > 32
            disp('Stuck backtracking, returning last iterate.  (See Section 4 of notes for more information.)')
            st=1;
            xp=x;
            return
        end
    end

    % next iteration
    x = xp; 
    u = up;
    Ax = Axp; 
    Atv = Atvp;
    lamu1 = lamu1p;  
    lamu2 = lamu2p;
    fu1 = fu1p; 
    fu2 = fu2p;

    % surrogate duality gap
    sdg = -(fu1'*lamu1 + fu2'*lamu2);
    tau = mu*2*M/sdg;
    rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
    rdual = rdp;
    resnorm = norm([rdual; rcent]);

    done = (sdg < 1e-8) | (pditer >= 50);

end