%% Initialize Variables
A = sym([ 1/4       0        0       0      0
          1/2       1/4      0       0      0
          17/50    -1/25     1/4     0      0
          371/1360 -137/2720 15/544  1/4    0
          25/24    -49/48    125/16 -85/12  1/4]);
s = length(A);  p = 4;       %Stage number and order
b = A(s,:)';    B = diag(b); 
e = ones(s,1);  I = eye(s);
M = e;
for i = 1:s-1
    M = [M A^i*e];
end

disp('Check tall tree conditions')
fact = sym(factorial((1:p)'));
fact_recip = 1./fact;
disp(logical(M(:,1:p)'*b==fact_recip)')

%% E-polynomial Conditions
syms z
syms y real
% Determine N and D functions
[N(z),D(z)] = numden(1+z*b'*(I-z*A)^(-1)*e);

% Create E-poly
E(y) = collect(expand(D(1i*y)*D(-1i*y)-N(1i*y)*N(-1i*y)))

% Standardize E-poly
k=coeffs(E(y));
P = diag(k)
Pd=double(P);
m = length(Pd);

%Create monomial vector
ys = 1;
for i=1:m-1
    ys = [ys;y^i];
end

% Define N such that y'Ny==0
n = triu(sym('n',[m m],'real'));
N=sym(zeros(m,m));
Im=eye(m);
for i=1:m-2
    for j=i+2:m
        % l=m*(i-1)-.5*i*(i+3)+j
        c=ceil((i+j)/2);
        f=floor((i+j)/2);
        N = N + n(i,j)*(Im(:,i)*Im(:,j)'+Im(:,j)*Im(:,i)'-Im(:,f)*Im(:,c)'-Im(:,c)*Im(:,f)');
    end
end

%Create symbolic function for N
inputs = sort(n(triu(n,2)~=0));
Ns=symfun(N,inputs);

%Create double precision function for N
Nd = matlabFunction(N);

%E-poly SDP
cvx_precision high 
cvx_begin sdp quiet
    variable eta
    minimize 1
    subject to
        Pd+Nd(eta)>=0;
cvx_end

%Symbolic LDL factorization of F
etaF=round(sym(eta),-1)
F=P+Ns(etaF)
[LF,DF]=ldls(F)

%% Verify LDL factoriztion
% all(all(logical(F==LF*DF*LF')))

%Verify A-stability conditions
% logical(expand(E(y)-y^6*ys'*F*ys)==0)
% all(logical(diag(DF)>=0))

%% CSTW Conditions
%Parameterize Ne=0
n = sym('n',[s s],'real');
n = triu(n);
N = sym(zeros(s,s));
for i = 1:s-1
    for j = i+1:s
        v=I(:,i)-I(:,j);
        N = N + n(i,j).*(v*v');
    end
end

%Define R and X wrt N
R = B+N;
X = R*A+A'*R-b*b';

%Parameterize nullvectors Xv=0 and update N,R,X
X_nullity = floor(p/2);
for i = 1:X_nullity
    in = struct2cell(solve(X*M(:,i), diag(n,i)));
    Ns = symfun(N, diag(n,i));
    N = Ns(in{:});
    R = B+N;
    X = R*A+A'*R-b*b';
end

%Degrees of freedom
dof = nchoosek(s-X_nullity,2);

%Create symbolic functions for R,X
inputs = sort(n(triu(n,X_nullity+1)~=0));
Rs=symfun(R,inputs);
Xs=symfun(X,inputs);

%Create double precision MATLAB functions for R,X
Rd = matlabFunction(R);
Xd = matlabFunction(X);

%SDP feasibility problem
cvx_begin sdp quiet
    variable eta(dof,1)
    minimize 1
    subject to
        in = num2cell(eta);
        Rd(in{:})>=0;
        Xd(in{:})>=0;
cvx_end

%Round eta to rational value, eta_, that satisfies Xs,Rs>=0
etaCSTW=sym(round(eta',2,'significant'))

%Symbolic LDL factorization of X,R
in=sym2cell(etaCSTW);
X=Xs(in{:})
[LX,DX] = ldls(X)
R=Rs(in{:})
[LR,DR] = ldls(R)

%% Verify LDL factorization
% all(all(logical(X==LX*DX*LX')))
% all(all(logical(R==LR*DR*LR')))

%Verify A-stability Conditions
% all(all(logical(X==R*A+A'*R-b*b')))
% all(logical(R*e==b))
% for i = 1:X_nullity
%     all(logical(X*M(:,i)==0))
% end
% all(logical(diag(DX)>=0))
% all(logical(diag(DR)>=0))

%% Clear Misc. Vars
clearvars -except etaF F LF DF etaCSTW R LR DR X LX DX

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
            if(D(j,j)~=0)
                L(i,j) = (M(i,j) - sumL)/D(j,j);
            end
        end
    end
end
