function [Pfail,inFid,r] = QSQnumericsSn
% QSQNUMERICSSN is a function that runs numerical analysis for 
% arXiv:2411.04215 for S_n model
% Output: 
%   Pfail -- average failing probability
%   inFid.val -- the largest infidelity among
%         the gates (average gate infidelity)
%         the state (infidelity)
%         the measurement (total var distance)
%   inFid.type -- specfies which of the infidelities is the largest:
%              = 1 if it's the state
%              = 2 if it's the gates
%              = 3 if it's the measurement

% -------------------------- main parameters ---------------------------
n = 4; % number of qubits
nn = 10^5; % number of samples
% levels of noise (relative to each other, should sum up to 1)
% the actual level of noise is adjusted such that maximal Pfail ~ Pfail_max
noise.U1 = 1; % for single-qubit unitary noise applied to gates
noise.Un = 0; % for n-qubit unitary noise applied to gates
noise.dep = 0; % depolarizing noise
noise.rho = 0; % for white noise applied to the state
noise.M = 0; % for the readout noise 
noise.alpha = .5; % scaling parameter (is adjusted automatically)
Pfail_max = 0.03; % maximal failling probability to be shown on the plot 
nmin = 20; % number of random points for the minimization
% set to 0 is no minimization is required

% -------------------------- parameters for determinig the ratio -------
eps = 10^(-5); % don't consider points below this threshold
ep = 10^(-3); % precision to which to determine the ratio
r = .1; % initial ratio

Pfail = zeros(1,nn);
inFid.val = zeros(1,nn);
inFid.type = zeros(1,nn);

% ------------------------- sequences, target model and outcomes -------
X = Sequences(n); % every cell contains a sequence
p = ones(1,numel(X))/numel(X); % distribution w.r.t. which sequences are drawn
Mod0 = TargetModel(n);
a_x = TargetOutcomes(X,Mod0,n); % every cell contains target outcomes

% adjust the noise level
noise = AdjustNoise(X,a_x,Mod0,noise,p,Pfail_max,n,100);

flag = 0;
disp(' Progress:   ')
for k=1:nn
    if round(4*k/nn)>flag
        flag = flag+1;
        noise.alpha = noise.alpha*0.5; % try to cover all ranges of noise levels
    end
    DisplayProgress(k,nn); % displays progress in percentage
    Mod = NoisyModel(noise,n,Mod0);
    % calculate the probability of the noisy model failing the protocol 
    Pfail(k) = 1-QuizzingProb(X,a_x,Mod)*p';
    % calculate the infidelity (partial unitary gauge optimization)
    if Pfail(k)>eps
        [inFid.val(k),r,inFid.type(k)] = Infidelity(Mod0,Mod,n,r,Pfail(k),ep,nmin);
        % here the ratio r is updated, if inFid.val is below the line
        % determined by r, no optimization is performed
    else
        % if the point is too close to the origin, simply igore it
        inFid.val(k) = r*eps; % (inFid_type(k) = 0)
    end
end

end

%% ---------- Quzzing probability ----------

function P = QuizzingProb(X,a_x,Mod)
% QUIZZINGPROB determines the quizzing probabilities
P = zeros(1,numel(X));
for j=1:numel(X)
    rho_j = Mod.rho;
    for l=1:numel(X{j})
        rho_j = ApplyChannel(rho_j,Mod.gates{X{j}(l)});
    end
    M_j = zeros(size(Mod.M{1}));
    for l=1:length(a_x{j})
        M_j = M_j+Mod.M{a_x{j}(l)};
    end
    P(j) = real(trace(M_j*rho_j));
end

end

%% ---------- Noise level adjustment ----------
function noise = AdjustNoise(X,a_x,Mod0,noise,p,Pfail_max,n,nn)
% ADJUSTNOISE automatically adjusts the noise level such that the maximal
% faling probability is (approx.) Pfail_max by sampling nn random models


alpha_d = noise.alpha/2;
while alpha_d>10^(-5)
% binary search
    Pfail_avg = 0;
    for k=1:nn
        Mod = NoisyModel(noise,n,Mod0);
        Pfail_avg = max(Pfail_avg,1-real(QuizzingProb(X,a_x,Mod))*p');
    end
    if Pfail_avg>Pfail_max
        noise.alpha = noise.alpha-alpha_d;
    else
        noise.alpha = noise.alpha+alpha_d;
    end
    alpha_d = alpha_d/2;
end
end

%% ---------- Target sequences ----------

function X = Sequences(n)
% SEQUENCES enumerates the sequences required for QSQ
% X is a cell array specifying the sequences
% e.g. X{1} = [1 2 1] corresponds to the sequence S_1 S_2 S_1
if n==1
    X = cell(0);
    X{1} = [];
    X{2} = [1 1];
    X{3} = [1 1 1 1];
else
    v = i2s(2*ones(1,n-1),1:2^(n-1))-1; % binary form of number 1:2^(n-1)
    X = cell(0);
    X{1} = [];
    for k=1:n
        % Squences X^loc_{N,k}
        v_k = [v(:,1:k-1),zeros(2^(n-1),1),v(:,k:end)];
        for j=1:2^(n-1) % go over all combinations for other qubits {0,2}
            seq = kron(find(v_k(j,:)),[1 1]);
            % now do all the combinations
            % X{end+1} = seq;
            X{end+1} = [seq, k];
            X{end+1} = [seq, k, k];
            X{end+1} = [seq, k, k, k];
            X{end+1} = [seq, k, k, k, k];
        end
        % sequence x_{N,k}
        X{end+1} = kron([1 1],[(k-1:-1:1),(n:-1:k)]);
        % sequence y_{N,k}
        X{end+1} = kron([1 1],[n:-1:1,k]);
    end
    % sequences X^G_{N,m}
    for m=1:n-1
        for l=m+1:n
            X{end+1} = [n:-1:m,l];
        end
    end
end
end

%% ---------- Target outcomes ----------

function a_x = TargetOutcomes(X,Mod0,n)
% TARGETOUTCOMES determine the target outcomes for the given sequences by
% calculating the probability distribution for the target model
a_x = cell(1,numel(X));
eps = .001/(2^n); % numerical precision (should be close to 0)

for j=1:numel(X)
    rho_j = Mod0.rho;
    for l=1:numel(X{j})
        rho_j = ApplyChannel(rho_j,Mod0.gates{X{j}(l)});
    end
    a_j = zeros(1,2^n);
    for l=1:2^n
        a_j(l) = real(trace(rho_j*Mod0.M{l}));
    end
    a_x{j} = find(a_j>eps);
end
end

%% ---------- Target model ----------

function Mod0 = TargetModel(n)
% TARGETGATES generates the target S_n model for n qubits
% Here as in the proof of S_n model, the basis is the computational one and
% the "S gate" is the S_y gate, i.e., sqrt(Y)
% gates:
Mod0.gates = cell(1,n);
S = [1 1i; -1i 1]/2+1i*[1 -1i; 1i 1]/2; % S_y gate
for k=1:n
    U_k = 1;
    for l=1:n
        if l==k
            U_k = kron(U_k,S);
        else
            U_k = kron(U_k,eye(2));
        end
    end
    Mod0.gates{k} = U_k;
end
% measurement:
v = i2s(2*ones(1,n),1:2^n);
m(:,:,1) = [1 0; 0 0];
m(:,:,2) = [0 0; 0 1];
Mod0.M = cell(1,2^n);
for k=1:2^n
    M_k = 1;
    for l=1:n
        M_k = kron(M_k,m(:,:,v(k,l)));
    end
    Mod0.M{k} = M_k;
end
% state:
Mod0.rho = Mod0.M{1};
end

%% ---------- Infidelity ----------

function [inFid,r,s] = Infidelity(Mod0,Mod,n,r,Pfail,ep,nn)
% INFIDELITY calculates the infidelity between the channels, states and the
% measurement
%   First, the infidelity is calculated. If it is below the one
%   corresponding to the worst ratio (up to ep), no optimization happens.
%
%   nn -- number of random point for the minimization
%   s -- type of the indifelity that is the largest

V = eye(2^n); % the initial unitary gauge (here the identity)

[~,inFid,s] = funMinPhases(zeros(1,2^n-1),Mod,Mod0,V,n);
if inFid>r*Pfail && nn>0 % the point is above the line given by the ratio r
    fun = @(x)funMinPhases(x,Mod,Mod0,V,n);
    options = optimoptions('fminunc','Display','off');
    fmin = 1;
    for k=1:nn
        [xk,fk] = fminunc(fun,rand(1,2^n-1),options);
        if fk<fmin
            x = xk;
            fmin = fk;
            [~,inFid,s] = funMinPhases(x,Mod,Mod0,V,n);
            if inFid<=r*Pfail*(1+ep)
                break
            end
        end
    end
    if inFid>r*Pfail % if it's still above
        r = inFid/Pfail; % update the coefficient
    end
end
end

%% ---------- Objective function Infidelity ----------

function [fmin,inFid,s] = funMinPhases(x,Mod,Mod0,V,n)
% FUNMINPHASES is the objective function for the minimization in Infidelity
d = 2^n;
% construct the gauge
U = V*diag([1,exp(1i*2*pi*x)]);

% calculate the infidelities
inFid_rho = real(1-trace(Mod.rho*U*Mod0.rho*U'));
inFid_gates = zeros(1,n);
for k=1:n
    % construct the Choi of the ideal gates
    inFid_gates(k) = (1-real(trace(CJ(Mod.gates{k},d)*CJ(U*Mod0.gates{k}*U',d))))*d/(d+1);
end
dist_M = 0; % ! this only works for the readout noise and V = eye(d) !
for k=1:d
    dist_M = dist_M+abs(trace(Mod0.rho*(Mod.M{k}-U*Mod0.M{k}*U')))/2;
end

fmin = (inFid_rho+sum(inFid_gates)/n+dist_M)/3; % function for the minimization
[inFid,s] = max([inFid_rho,max(inFid_gates),dist_M]);
% maybe add state and measurement
end

%% ---------- Noisy model ----------

function Mod = NoisyModel(noise,n,Mod0)
% NOISYMODEL generates a random noise added to state, channels and measurement
%   noise -- parameters controlling the ammount of noise (see description
%   above)
%
d = 2^n;
% state (white noise)
e = noise.rho*noise.alpha*rand(1); % a random level of noise between 0 and noise.rho
Mod.rho = (1-e)*Mod0.rho+e*eye(d)/d;
% measurement (stochastical noise)
v = i2s(2*ones(1,n),1:d);
e = noise.M*noise.alpha*rand(1,d); % a random level of noise for each qubit
m = cell(1,n);
for k=1:n
    m{k}(:,:,1) = [1-e(k) 0; 0 e(k)];
    m{k}(:,:,2) = [e(k) 0; 0 1-e(k)];
end
Mod.M = cell(1,d);
for k=1:d
    M_k = 1;
    for l=1:n
        M_k = kron(M_k,m{l}(:,:,v(k,l)));
    end
    Mod.M{k} = M_k;
end
% gates
Mod.gates = cell(1,n);
e_u1 = noise.U1*noise.alpha*rand(2,n,n); % level of single-qubit unitary noise
e_un = noise.Un*noise.alpha*rand(2,n); % level of n-qubit unitary noise
% the first index is = 1 for left and = 2 for right noise

for k=1:n
    % first the unitary noise
    Ul = 1; %left 
    Ur = 1; % right
    for l=1:n
        Ul = kron(Ul,RandUnitary(e_u1(1,l,k),2));
        Ur = kron(Ur,RandUnitary(e_u1(2,l,k),2));
    end
    Ul = Ul*RandUnitary(e_un(1,k),d);
    Ur = Ur*RandUnitary(e_un(2,k),d);
    Mod.gates{k} = Ur*Mod0.gates{k}*Ul;
end
if noise.dep~=0
    e_dep = noise.dep*noise.alpha*rand(1,n);
    E = eye(d);
    for k=1:n
        Mod.gates{k}(:,:,1) = sqrt(1-e_dep(k))*Mod.gates{k}(:,:,1);
        for i=1:d
            for j=1:d
                Mod.gates{k}(:,:,(i-1)*d+j+1) = sqrt(e_dep(k)/d)*E(:,i)*E(j,:);
            end
        end
    end
end
end

%% ---------- Random unitary ----------

function U = RandUnitary(alpha,d)
% RANDUNITARY generates a random unitary close to the identity
%   alpha -- parameter that determines how close U is to the identity
Z = randn(d)+1i*randn(d);
U = expm(alpha*1i*(Z+Z'));
end


%% ---------- Apply channel ----------

function rho = ApplyChannel(rho0,Lambda)
% APPLYCHANNEL applies channel Lambda to the state rho
% Lambda is given as a list of Kraus operators
rho = zeros(size(rho0));
for k=1:size(Lambda,3)
    rho = rho+Lambda(:,:,k)*rho0*Lambda(:,:,k)';
end
end


%% ---------- Choi-Jamiolkowski state ----------

function J = CJ(Lambda,d)
% CJ calculates the Choi-Jamiolkowski state of the channel Lambda
J = zeros(d^2);
phi = reshape(eye(d),[],1);

for k=1:size(Lambda,3)
    phi_k = kron(Lambda(:,:,k),eye(d))*phi;
    J = J+phi_k*phi_k'/d;
end
end

%% ---------- Index-to-string ----------

function s = i2s(D,in)
% I2S is essentially the same as ind2sub, but faster. 
%    D -- vector of dimensions.
%    For D = [2 2 .. 2] use I2V, which is much faster.

if any(in>prod(D))
    error('Index cannot exceed prod(D).');
end
if any(in<=0)
    error('Index must be a positive integer.');
end
D = fliplr(D);
s = zeros(numel(in),numel(D));
for k=1:numel(in)
    for n=numel(D)-1:-1:1
        r = rem(in(k)-1,prod(D(1:n)))+1;
        s(k,n+1) = (in(k)-r)/prod(D(1:n))+1;
        in(k) = in(k) - (s(k,n+1)-1)*prod(D(1:n));
    end
    s(k,1) = rem(in(k)-1,D(1))+1;
    s(k,:) = fliplr(s(k,:));
end
end

%% ---------- Display progress ----------

function DisplayProgress(k,N)
% DISPLAYPROGRESS is used to display progress in the precentage while k
% runs from 1 to N.
intervals_k = floor(N*(.01:.01:1));
if isempty(find(intervals_k==k,1))==0 % if k is one of the intervals
    percent_to_diplay = floor(100*k/N);
    if percent_to_diplay==floor(100/N)
        fprintf(1,' %i%%',percent_to_diplay);
    elseif percent_to_diplay<=10
        fprintf(1,'\b\b%i%%',percent_to_diplay);
    elseif percent_to_diplay==100
        fprintf(1,'\b\b\b\b\n');
    else
        fprintf(1,'\b\b\b%i%%',percent_to_diplay);
    end
end
end
