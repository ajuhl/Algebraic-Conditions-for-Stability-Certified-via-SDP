%Initialize Variables
A=[sym('1290066345260422/10000000000000000')   0 0 0 0 0 0
   sym('3315354455306989/10000000000000000')   sym('1177478680001996/10000000000000000')    0 0 0 0 0
   sym('-8009819642882672/100000000000000000') sym('-2408450965101765/1000000000000000000') sym('9242630648045402/100000000000000000')  0 0 0 0
   sym('-1730636616639455/1000000000000000')   sym('1513225984674677/1000000000000000')     sym('1221258626309848/1000000000000000')    sym('2266279031096887/10000000000000000')   0 0 0  
   sym('1475353790517696/10000000000000000')   sym('3618481772236499/10000000000000000')    sym('-5603544220240282/10000000000000000')  sym('2455453653222619/1000000000000000')    sym('5742190161395324/10000000000000000')     0 0
   sym('2099717815888321/10000000000000000')   sym('7120237463672882/10000000000000000')    sym('-2012023940726332/100000000000000000') sym('-1913828539529156/100000000000000000') sym('-5556044541810300/1000000000000000000')  sym('3707277349712966/10000000000000000')  0
   sym('2387938238483883/10000000000000000')   sym('4762495400483653/10000000000000000')    sym('12339351512133/1000000000000000')      sym('6011995982693821/100000000000000000')  sym('6553618225489034/100000000000000000000') sym('-1270730910442124/10000000000000000') sym('3395048796261326/10000000000000000')];
s = length(A);  p = 4;  % number of stages & order of method
e = ones(s,1);  I = eye(s); 
b = A(s,:)';
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

% Standardize E-polynomial
kend = k(end);
k=k./kend;
P = diag(k);
Pd=double(P);
m = length(P);
dof=nchoosek(m-1,2);

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
    variable gam
    variable eta(dof,1)
    minimize 1
    subject to
        in = num2cell(eta);
        Pd+Nd(in{:})>=0;
cvx_end

% Symbolic LDL factorization of F
etaF=sym(round(eta,-1))'
in=sym2cell(etaF);
F=P+Ns(in{:})
[LF,DF]=ldls(F)

%% Verify LDL factoriztion
% all(all(logical(F==LF*DF*LF')))

%Verify A-stability conditions
% logical(expand(E(y)-kend*y^2*ys'*F*ys)==0)
% all(logical(diag(DF)>=0))

%% Clear Misc. Vars
clearvars -except etaF F LF DF

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