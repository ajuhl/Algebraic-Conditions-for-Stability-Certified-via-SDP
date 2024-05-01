%% Initialize variables
A = [0 0 0 0 0 0 0 0 0
     sym('1/6')                     sym('1/6') 0 0 0 0 0 0 0
     sym('11/96')                   sym('-1/32')                  sym('1/6') 0 0 0 0 0 0
     sym('1/12')                    sym('-1/4')                   sym('1/2')                      sym('1/6') 0 0 0 0 0
     sym('-2015/15072')             sym('-6987/5024')             sym('3271/1884')                sym('175/471')                  sym('1/6') 0 0 0 0
     sym('-326531/573678')          sym('-114988/31871')          sym('1208156/286839')           sym('132950/286839')            sym('68/203')             sym('1/6') 0 0 0
     sym('-331717945/2106545616')   sym('-480525599/416107776')   sym('2240951089/1404363744')    sym('394951619/2808727488')     sym('-5160553/26834976')  sym('35815/352512')     sym('1/6') 0 0
     sym("16264655341/73026914688") sym('9786099235/14425069568') sym('-34306812733/48684609792') sym('-15985588007/97369219584') sym('37652437/930279168') sym('-340747/12220416') sym('1/26') sym('1/6') 0
     sym('7/90')                    0                             0                               0                               sym('16/45')              sym('-4/45')            sym('2/15') sym('16/45') sym('1/6')];
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
syms alpha beta r y real
% Determine N and D functions
[N(z),DF(z)] = numden(1 + z*b'*(I - z*A)^(-1)*e);

% Create E-poly
E(y) = collect(expand(DF(1i*y)*DF(-1i*y) - N(1i*y)*N(-1i*y)));

disp('Check sign of E(pi/2)-polynomial coefficients')
k = coeffs(E(y));
disp(logical(k >= 0))

% Determine N and D functions
[N(z),DF(z)] = numden(1+z*b'*(I-z*A)^(-1)*e);

% Create E-poly
r=y^2;
E(y,alpha) = collect(simplify(expand(DF(-r*(cos(alpha) + 1i*sin(alpha)))*DF(-r*(cos(alpha) - 1i*sin(alpha))) - N(-r*(cos(alpha) + 1i*sin(alpha)))*N(-r*(cos(alpha) - 1i*sin(alpha))))));
E(y,beta) = collect(subs(E(y,alpha),{cos(alpha)},{beta}))

% Standardize E-poly
bta=sym('2218472195/100000000000')
k = coeffs(E(y,bta),y);

disp('Check sign of E(beta*)-polynomial coefficients')
disp(logical(k >= 0))

% Standardize E-poly
kend = k(end);
P = diag(k./kend);
Pd = double(P);
m = length(Pd);

%Create monomial vector
ys = 1;
for i = 1:m-1
    ys = [ys;y^i];
end

% Define N such that y'Ny==0
n = triu(sym('n',[m m],'real'));
N = sym(zeros(m,m));
Im = eye(m);
for i = 1:m-2
    for j = i+2:m
        % l = m*(i-1) - .5*i*(i+3) + j
        c = ceil((i+j)/2);
        f = floor((i+j)/2);
        if mod(i+j,2)   %Preset certain paramerters to zero
            n(i,j)=0;
        end
        N = N + n(i,j)*(Im(:,i)*Im(:,j)' + Im(:,j)*Im(:,i)' - Im(:,f)*Im(:,c)' - Im(:,c)*Im(:,f)');
    end
end

% Create symbolic function for N
inputs = sort(n(triu(n,2) ~= 0));
Ns = symfun(N,inputs);

% Create double precision function for N
Nd = matlabFunction(N);
dof = length(inputs);

% E-poly SDP
cvx_precision high 
cvx_begin sdp quiet
    variable eta(dof,1)
    minimize 1
    subject to
        in = num2cell(eta);
        Pd + Nd(in{:}) >=  0;
cvx_end

% Symbolic LDL factorization of F
etaF = str2sym(rat(eta,1e-6))'
in = sym2cell(etaF);
F = P + Ns(in{:});
[LF,DF] = ldls(F);

%% Verify LDL factoriztion
% all(all(logical(F == LF*DF*LF')))

% Verify A-stability conditions
% logical(expand(E(y,bta) - kend*y^2*ys'*F*ys) == 0)
% all(logical(diag(DF) >= 0))

%% Clear Misc. Vars
clearvars -except etaF F LF DF bta

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