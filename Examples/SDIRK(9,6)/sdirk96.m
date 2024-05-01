% Initialize Variables
A=[sym('218127781944908/1000000000000000'),   0,0,0,0,0,0,0,0
  -sym('903514856119419/10000000000000000'),  sym('218127781944908/1000000000000000'),   0,0,0,0,0,0,0
   sym('172952039138937/1000000000000000'),  -sym('353655010362820/1000000000000000'),   sym('218127781944908/1000000000000000'), 0,0,0,0,0,0
   sym('511999875919193/1000000000000000'),   sym('289640332201925/10000000000000000'), -sym('144030945657094/10000000000000000'),sym('218127781944908/1000000000000000'),  0,0,0,0,0
   sym('465303495506782/100000000000000000'),-sym('075635818766597/1000000000000000'),   sym('217273030786712/1000000000000000'),-sym('206519428725472/10000000000000000'), sym('218127781944908/1000000000000000'),  0,0,0,0
   sym('896145501762472/1000000000000000'),   sym('139267327700498/1000000000000000'),  -sym('186920979752805/1000000000000000'), sym('672971012371723/10000000000000000'),-sym('350891963442176/1000000000000000'),  sym('218127781944908/1000000000000000'), 0,0,0
   sym('552959701885751/1000000000000000'),  -sym('439360579793662/1000000000000000'),   sym('333704002325091/1000000000000000'),-sym('339426520778416/10000000000000000'),-sym('151947445912595/1000000000000000'),  sym('213825661026943/10000000000000000'),sym('218127781944908/1000000000000000'), 0,0
   sym('631360374036476/1000000000000000'),   sym('724733619641466/1000000000000000'),  -sym('432170625425258/1000000000000000'), sym('598611382182477/1000000000000000'), -sym('709087197034345/1000000000000000'), -sym('483986685696934/1000000000000000'), sym('378391562905131/1000000000000000'), sym('218127781944908/1000000000000000'),0
   0,                                        -sym('155044525308690/1000000000000000'),   sym('194518478660789/1000000000000000'), sym('635156402792030/1000000000000000'),  sym('811722786641730/1000000000000000'),  sym('110736108691585/1000000000000000'),-sym('495304692414479/1000000000000000'),-sym('319912341007872/1000000000000000'),sym('218127781944908/1000000000000000')];
s = length(A);  p = 6;  % number of stages & order of method
e = ones(s,1);  I = eye(s);
b = A(s,:)';
M = e;
for i = 1:s-1
    M = [M A^i*e];
end

disp('Check tall tree conditions')
fact = sym(factorial((1:p)'));
fact_recip = 1./fact;
disp(logical(M(:,1:p)'*b == fact_recip)')

%% Strategy 1
syms z
syms y real
% Determine N and D functions
[N(z),D(z)] = numden(1 + z*b'*(I - z*A)^(-1)*e);

% Create E-poly
E(y) = collect(expand(D(1i*y)*D(-1i*y) - N(1i*y)*N(-1i*y)));

disp('Check sign of E-polynomial coefficients')
k = coeffs(E(y));
disp(logical(k >= 0))

% Save original coefficients
A0 = A;
b0 = b;

%% Strategy 2
% Rational Approximation of A
for i = 1:s
    for j = 1:i
        for k = 16:20
            a = str2sym(rat(A0(i,j),10^-k));
            if double(a) == double(A0(i,j))
                break
            end
        end
        A(i,j) = a;
    end
end

% Determine b by solving for tall tree conditions
b = sym('b',[s 1],'real');
b(1) = sym(0);
b(s) = A(s,s);
A(s,:) = b';
M = e;
for i = 1:s-1
    M = [M A^i*e];
end
sol = struct2cell(solve(M(:,1:p)'*b == fact_recip,b(3:s-1)));
bs = symfun(b,b(3:s-1));
b = bs(sol{:});

% Optimize for remaining free parameter
bs = symfun(b,b(2));
in = solve(diff(norm(bs - b0)));
b = bs(in);
B = diag(b);
A(s,:) = b';
M = e;
for i = 1:s-1
    M = [M A^i*e];
end
A_tilde = A
b_tilde = b

disp('Check tall tree conditions')
disp(logical(M(:,1:p)'*b == fact_recip)')

disp('Perturbation error bounds')
eps_A = double(max(max(abs(A - A0))))
eps_b = double(max(abs(b - b0)))

%% E-polynomial Conditions
syms z
syms y real
% Determine N and D functions
[N(z),D(z)] = numden(1 + z*b'*(I - z*A)^(-1)*e);

% Create E-poly
E(y) = collect(expand(D(1i*y)*D(-1i*y) - N(1i*y)*N(-1i*y)));

% Standardize E-poly
k = coeffs(E(y));
disp('Check sign of E_tilde-polynomial coefficients')
disp(logical(k >= 0))
P_tilde = diag(k)

%% CSTW Conditions
% Parameterize Ne=0
n = sym('n',[s s],'real');
n = triu(n);
N = sym(zeros(s,s));
for i = 1:s-1
    for j = i+1:s
        v = I(:,i) - I(:,j);
        N = N + n(i,j).*(v*v');
    end
end

% Define R and X wrt N
R = B + N;
X = R*A + A'*R - b*b';

% Parameterize nullvectors Xv=0 and update N,R,X
X_nullity = floor(p/2);
for i = 1:X_nullity
    in = struct2cell(solve(X*M(:,i), diag(n,i)));
    Ns = symfun(N, diag(n,i));
    N = Ns(in{:});
    R = B + N;
    X = R*A + A'*R - b*b';
end

% Degrees of freedom
dof = nchoosek(s-X_nullity,2);

% Create symbolic functions for R,X
inputs = sort(n(triu(n,X_nullity+1) ~= 0));
Rs = symfun(R,inputs);
Xs = symfun(X,inputs);

% Create double precision MATLAB functions for R,X
Rd = matlabFunction(R);
Xd = matlabFunction(X);

% SDP feasibility problem
cvx_begin sdp quiet
    variable eta(dof,1)
    minimize 1
    subject to
        in = num2cell(eta);
        Rd(in{:}) >= 0;
        Xd(in{:}) >= 0;
cvx_end

% Round eta to rational value that satisfies Xs,Rs>=0
etaCSTW = str2sym(rat(eta,1e-3))'

% Symbolic LDL factorization of X,R
in = sym2cell(etaCSTW);
X = Xs(in{:})
[LX,DX] = ldls(X)
R = Rs(in{:})
[LR,DR] = ldls(R)

%% Verify LDL factorization
% all(all(logical(X == LX*DX*LX')))
% all(all(logical(R == LR*DR*LR')))

% Verify A-stability Conditions
% all(all(logical(X == R*A + A'*R - b*b')))
% all(logical(R*e == b))
% for i = 1:X_nullity
%     all(logical(X*M(:,i) == 0))
% end
% all(logical(diag(DX) >= 0))
% all(logical(diag(DR) >= 0))

%% Clear Misc. Vars
clearvars -except A_tilde b_tilde P_tilde etaCSTW R LR DR X LX DX

%% Symbolic LDL factorization
function [L,D] = ldls(M)

    % Initialize L and D
    l = length(M);
    L = sym(eye(l));
    D = sym(zeros(l));
    
    % LDL^T decomposition
    for j = 1:l
        % Diagonal D entries
        sumD = 0;
        for k = 1:j-1
            sumD = sumD + L(j,k)^2*D(k,k);
        end
        D(j,j) = M(j,j) - sumD;
        
        % Off-diagonal L entries
        for i = j+1:l
            sumL = 0;
            for k = 1:j-1
                sumL = sumL + L(i,k)*L(j,k)*D(k,k);
            end
            if(D(j,j) ~= 0)
                L(i,j) = (M(i,j) - sumL)/D(j,j);
            end
        end
    end
end
