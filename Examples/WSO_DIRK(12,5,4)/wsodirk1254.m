%Initialize Variables
A=[sym('2345371908646273/10000000000000000')   0 0 0 0 0 0 0 0 0 0 0
   sym('6874344413888787/10000000000000000')   sym('5515270980695153/100000000000000000')  0 0 0 0 0 0 0 0 0 0
   sym('-1183552669539587/10000000000000000')  sym('5463563002913454/1000000000000000000') sym('1458584459918280/10000000000000000')   0 0 0 0 0 0 0 0 0
   sym('-1832235204042292/10000000000000000')  sym('5269029412008775/100000000000000000')  sym('8203685085133529/10000000000000000')   sym('4812118949092085/100000000000000000')   0 0 0 0 0 0 0 0
   sym('9941572060659400/100000000000000000')  sym('4977904930055774/1000000000000000000') sym('5414758174284321/100000000000000000')  sym('-1666571741820749/1000000000000000000') sym('8078975617332473/100000000000000000')  0 0 0 0 0 0 0
   sym('-9896614721582678/10000000000000000')  sym('2860682690577833/1000000000000000')    sym('-1236119341063179/1000000000000000')   sym('2130219523351530/1000000000000000')     sym('-1260655031676537/1000000000000000')   sym('2457717913099987/10000000000000000')   0 0 0 0 0 0
   sym('-5656238413439102/100000000000000000') sym('1661985685769353/10000000000000000')   sym('6464600922362508/10000000000000000')   sym('6608854962269927/10000000000000000')    sym('3736054198873429/10000000000000000')   sym('6294456964407685/10000000000000000')   sym('5702752607818027/10000000000000000')  0 0 0 0 0
   sym('8048962104724392/10000000000000000')   sym('-6232034990249100/100000000000000000') sym('5737234603323347/10000000000000000')   sym('-9613723511489970/100000000000000000')  sym('5524106361737929/10000000000000000')   sym('5961002486833255/10000000000000000')   sym('1978411600659203/10000000000000000')  sym('3156238724024008/10000000000000000')    0 0 0 0 
   sym('-1606381759216300/10000000000000000')  sym('6833397073337708/10000000000000000')   sym('4734578665308685/10000000000000000')   sym('8037708984872738/10000000000000000')    sym('-1094498069459834/100000000000000000') sym('6151263362711297/10000000000000000')   sym('3908946848682723/10000000000000000')  sym('8966103265353116/100000000000000000')   sym('2973255537857041/100000000000000000')  0 0 0
   sym('7074283235644631/10000000000000000')   sym('4392037300952482/10000000000000000')   sym('-3623592480237268/100000000000000000') sym('7189990308645932/10000000000000000000') sym('5820968279166545/10000000000000000')   sym('3302003177175218/10000000000000000')   sym('-2394564021215881/10000000000000000') sym('-7540283547997615/1000000000000000000') sym('1702137469523672/10000000000000000')   sym('6268780138721711/10000000000000000')   0 0
   sym('1361197981133694/10000000000000000')   sym('-7486549901902831/10000000000000000')  sym('1893908350024949/1000000000000000')    sym('3940485196730028/10000000000000000')    sym('6240233526545023/100000000000000000')  sym('7511983862200027/10000000000000000')   sym('-5283465265730526/10000000000000000') sym('-1661625677872943/1000000000000000')    sym('9998723833190827/10000000000000000')   sym('1377776742457387/1000000000000000')    sym('8905676409277480/10000000000000000')    0
   sym('-7433675378768276/10000000000000000')  sym('1490594423766965/10000000000000000')   sym('-2042884056742363/100000000000000000') sym('8565329438087443/10000000000000000000') sym('1357261590983184/1000000000000000')    sym('2067512027776675/1000000000000000000') sym('9836884265759428/100000000000000000') sym('-1357936974507222/100000000000000000')  sym('-5428992174996300/100000000000000000') sym('-3803299038293005/100000000000000000') sym('-9150525836295019/1000000000000000000') sym('2712352651694511/10000000000000000')];

s = length(A);  p = 5;  % number of stages order of method
b = A(s,:)';    B=diag(b);
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

% Standardize E-poly
k=coeffs(E(y));
logical(k>=0);
F=diag(k./k(end))

%% Clear Misc. Vars
clearvars -except F