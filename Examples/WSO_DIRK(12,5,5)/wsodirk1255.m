%% Initialize Variaibles
A=[sym('4113473525867655/100000000000000000')  0 0 0 0 0 0 0 0 0 0 0
   sym(' 1603459327727949/10000000000000000')  sym(' 6663913326722831/100000000000000000') 0 0 0 0 0 0 0 0 0 0
   sym('-3424389044264752/10000000000000000')  sym(' 8658006324816373/10000000000000000')  sym(' 9893519116923277/100000000000000000') 0 0 0 0 0 0 0 0 0
   sym(' 9437182028870806/1000000000000000')   sym('-1088783359642350/100000000000000')    sym(' 2644025436733866/1000000000000000')   sym(' 1846155800500574/10000000000000000')   0 0 0 0 0 0 0 0
   sym('-3425409029430815/10000000000000000')  sym(' 5172239272544332/10000000000000000')  sym(' 9163589909678043/10000000000000000')  sym(' 5225142808845742/100000000000000000')  sym(' 1165485436026433/10000000000000000')  0 0 0 0 0 0 0
   sym('-2094441177460360/1000000000000000')   sym(' 2577655753533404/1000000000000000')   sym(' 5704405293326313/10000000000000000')  sym(' 1213637180023516/10000000000000000')   sym('-4752289775376601/10000000000000000')  sym('5285605969257756/10000000000000000')  0 0 0 0 0 0
   sym(' 3391631788320480/10000000000000000')  sym('-2797427027028997/10000000000000000')  sym(' 1039483063369094/1000000000000000')   sym(' 5978770926212172/100000000000000000')  sym('-2132900327070380/10000000000000000')  sym('8344318363436753/100000000000000000') sym('2410106515779412/10000000000000000') 0 0 0 0 0
   sym(' 5904282488642163/1000000000000000')   sym(' 3171195765985073/1000000000000000')   sym('-1236822836316587/100000000000000')    sym('-4989519066913001/10000000000000000')   sym(' 2160529620826442/1000000000000000')   sym('1916104322021480/1000000000000000')   sym('1988059486291180/1000000000000000')  sym(' 2232092386922440/10000000000000000')   0 0 0 0
   sym(' 4616443509508975/10000000000000000')  sym('-1933433560549238/10000000000000000')  sym('-1212541486279519/10000000000000000')  sym(' 6662362039716674/100000000000000000')  sym(' 4254912950625259/10000000000000000')  sym('7856131647013712/10000000000000000')  sym('8369551389357689/10000000000000000') sym(' 1604780447895926/10000000000000000')   sym('3616125951766939/10000000000000000')  0 0 0
   sym('-7087669749878204/10000000000000000')  sym(' 6466527094491541/10000000000000000')  sym(' 4758821526542215/10000000000000000')  sym('-2570518451375722/10000000000000000')   sym(' 1123185062554392/1000000000000000')   sym('5546921612875290/10000000000000000')  sym('3192424333237050/10000000000000000') sym(' 3612077612576969/10000000000000000')   sym('5866779836068974/10000000000000000')  sym(' 2353799736246102/10000000000000000')   0 0
   sym(' 4264162484855930/10000000000000000')  sym(' 1322816663477840/1000000000000000')   sym(' 4245673729758231/10000000000000000')  sym('-2530402764527700/1000000000000000')    sym('-7822016897497742/100000000000000000') sym('1054463080605071/1000000000000000')   sym('4645590541391895/10000000000000000') sym(' 1145097379521439/1000000000000000')    sym('4301337846893282/10000000000000000')  sym(' 1499513057076809/1000000000000000')    sym('1447942640822165/100000000000000000')   0
   sym(' 1207394392845339/100000000000000000') sym(' 5187080074649261/10000000000000000')  sym(' 1121304244847239/10000000000000000')  sym('-4959806334780896/1000000000000000000') sym('-1345031364651444/1000000000000000')   sym('3398828703760807/10000000000000000')  sym('8159251531671077/10000000000000000') sym('-2640104266439604/1000000000000000000') sym('1439060901763520/100000000000000000') sym('-6556567796749947/1000000000000000000') sym('6548135446843367/10000000000000000000') sym('5454220210658036/10000000000000000')];

s = length(A);  p = 5;  % number of stages & order of method
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
b1 = A(s,:)';
b = sym('b',[s 1],'real');
A(s,:) = b';
M = e;
for i = 1:s-1
    M = [M A^i*e];
end
sol = struct2cell(solve(M(:,1:p)'*b == fact_recip,b(s-p+1:s)));
bs = symfun(b,b(s-p+1:s));
b = bs(sol{:});

bs = symfun(b,b(1:s-p));
in = sym2cell(b1(1:s-p));
b = bs(in{:});
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

% Standardize E-poly
kend = k(end);
P = diag(k./kend);
Pd = double(P);
m = length(P);
dof = nchoosek(m-1,2);

% Create monomial vector
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
        N = N + n(i,j)*(Im(:,i)*Im(:,j)' + Im(:,j)*Im(:,i)' - Im(:,f)*Im(:,c)' - Im(:,c)*Im(:,f)');
    end
end

% Create symbolic function for N
inputs = sort(n(triu(n,2) ~= 0));
Ns = symfun(N,inputs);

% Create double precision function for N
Nd = matlabFunction(N);

% E-poly SDP
cvx_precision high 
cvx_begin sdp quiet
    variable eta(dof,1)
    minimize 1
    subject to
        in = num2cell(eta);
        (Pd+Nd(in{:}))*1e-9 >= 0;
cvx_end

% Symbolic LDL factorization of F
etaF = sym(round(eta,-3))'
in = sym2cell(etaF);
F = P + Ns(in{:})
[LF,DF] = ldls(F)

%% Verify LDL factoriztion
% all(all(logical(F == LF*DF*LF')))

% Verify A-stability conditions
% logical(expand(E(y) - kend*y^6*ys'*F*ys) == 0)
% all(logical(diag(DF) >= 0))

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
            if(D(j,j) ~= 0)
                L(i,j) = (M(i,j) - sumL)/D(j,j);
            end
        end
    end
end

