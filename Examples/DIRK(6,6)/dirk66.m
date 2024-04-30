%Initialize Variables
A=[sym('303487844706747/1000000000000000'),  0,0,0,0,0
   sym('-279756492709814/1000000000000000'), sym('500032236020747/1000000000000000'),  0,0,0,0
   sym('280583215743895/1000000000000000'),  sym('-438560061586751/1000000000000000'), sym('217250734515736/1000000000000000'),  0,0,0
   sym('-677678738539846/10000000000000000'),sym('984312781232293/1000000000000000'),  sym('-266720192540149/1000000000000000'), sym('2476680834526/10000000000000'),     0,0
   sym('125671616147993/1000000000000000'),  sym('-995401751002415/1000000000000000'), sym('761333109549059/1000000000000000'),  sym('-210281837202208/1000000000000000'),sym('866743712636936/1000000000000000'),0
   sym('-368056238801488/1000000000000000'), sym('-999928082701516/1000000000000000'), sym('534734253232519/1000000000000000'),  sym('-174856916279082/1000000000000000'),sym('615007160285509/1000000000000000'),sym('696549912132029/1000000000000000')];
b=[sym('257561510484877/1000000000000000')
   sym('234281287047716/1000000000000000')
   sym('126658904241469/1000000000000000')
   sym('252363215441784/1000000000000000')
   sym('396701083526306/1000000000000000')
   sym('-267566000742152/1000000000000000')];
s = length(A);  p = 6;  % number of stages & order of method
e = ones(s,1);  I = eye(s);
M = e;
for i = 1:s-1
    M = [M A^i*e];
end

disp('Check tall tree conditions')
fact = sym(factorial((1:p)'));
fact_recip = 1./fact;
disp(logical(M(:,1:p)'*b==fact_recip)')

%% Strategy 1
syms z
syms y real
% Determine N and D functions
[N(z),D(z)] = numden(1+z*b'*(I-z*A)^(-1)*e);

% Create E-poly
E(y) = collect(expand(D(1i*y)*D(-1i*y)-N(1i*y)*N(-1i*y)));

disp('Check sign of E-polynomial coefficients')
k=coeffs(E(y));
disp(logical(k>=0))

%Save original coefficients
A0=A;
b0=b;

%% Strategy 2
%Rational Approximation of A
for i = 1:s
    for j = 1:i
        for k = 16:30
            a = str2sym(rat(A0(i,j),10^-k));
            if double(a)==double(A0(i,j))
                break
            end
        end
        A(i,j)=a;
    end
end
A_tilde=A

% Create new M
M = e;
for i = 1:s-1
    M = [M A^i*e];
end

% Determine b by solving for tall tree conditions
b=sym('b',[s 1],'real');
b=cell2sym(struct2cell(solve(M(:,1:p)'*b==fact_recip)));
b_tilde = b
B=diag(b);

disp('Check tall tree conditions')
disp(logical(M(:,1:p)'*b==fact_recip)')

disp('Perturbation error bounds')
eps_A=double(max(max(abs(A-A0))))
eps_b=double(max(abs(b-b0)))

%% E-polynomial Conditions
syms z
syms y real
% Determine N and D functions
[N(z),D(z)] = numden(1+z*b'*(I-z*A)^(-1)*e);

% Create E-poly
E(y) = collect(expand(D(1i*y)*D(-1i*y)-N(1i*y)*N(-1i*y)));

% Standardize E-poly
k=coeffs(E(y));
disp('Check sign of E_tilde-polynomial coefficients')
disp(logical(k>=0))
P_tilde=diag(k)

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

%Round eta to rational value that satisfies Xs,Rs>=0
etaCSTW=sym(round(eta'))

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
            if(D(j,j)~=0)
                L(i,j) = (M(i,j) - sumL)/D(j,j);
            end
        end
    end
end
