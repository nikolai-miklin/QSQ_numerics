function [Pfail,inFid,r] = QSQnumericsCl2(n_sample)
% QSQNUMERICSCL2 is a function that runs numerical analysis for 
% arXiv:2411.04215 for Cl_2 model
% 
% This program samples random noisy models for plots of the dependency
% between model infidelity vs. failing probability.
% Partial unitary gauge optimization is performed. In the paper, this
% corresponds to the final step of the proof for S_n model where the global
% phases are fixed.
% Each noise model is parameterized by the noise strength, which is a single
% general parameter, alpha, which is takes n_alpha values from 0 to alpha_max
% Input:
%   n_sample -- number of randomly sampled noisy models
% Output: 
%   Pfail -- average failing probability
%   inFid.max -- the largest infidelity among
%         the gates (average gate infidelity)
%   inFid.avg -- the average gate set infidelity
%   inFid.rho -- infidelity for the initial state
%   inFid.meas -- the total variation distance for the measurement
%   inFid.type -- specifies which of the infidelities is the largest:
%              = 1 if it's the gates (averaged)
%              = 2 if it's the state
%              = 3 if it's the measurement
%   r.avg -- the worst ratio between the average gate set infidelity and
%   the failing probability
%   r.max -- the same for the worst infidelity among the gates

n = 2; % this program is for 2-qubit Clifford model only
% fix the general parameters for noise
noise_param.U = 5; % parameter for misalignment
noise_param.eig = .1; % parameter for over-rotation noise
noise_param.list = [1:4]; % list of gates to which the noise is applied
noise_param.dep = 0; % depolarizing noise
noise_param.rho = 0; % white noise applied to the state
noise_param.M = 0; % readout noise 

% pre-allocate the memory
Pfail = zeros(1,n_sample);
inFid.max = zeros(1,n_sample);
inFid.type = zeros(1,n_sample);
inFid.avg = zeros(1,n_sample);
inFid.rho = zeros(1,n_sample);
inFid.meas = zeros(1,n_sample);

% ------------------------- sequences, target model and outcomes -------
X = Sequences(n); % every cell contains a sequence
p = ones(1,numel(X))/numel(X); % distribution w.r.t. which sequences are drawn
Mod0 = TargetModel(n);
a_x = TargetOutcomes(X,Mod0,n); % every cell contains target outcomes

% general noise scaling parameter alpha
alpha_max = .1;
n_alpha = 10;
d_alpha = alpha_max/n_alpha;
n_noise = round(n_sample/n_alpha); % number of different noise models
% gauge optimization parameters
r0 = 2; % points below this are not optimized
n_opt = 5; % number of random initial points for partial gauge optimization

disp(' Progress:   ')
for k=1:n_noise
    DisplayProgress(k,n_noise); % displays progress in percentage
    noise = Noise(n,numel(Mod0.gates)); % generate the noise model
    % now take n_alpha points between d_alpha and alpha_max
    for l=1:n_alpha
        ind = n_alpha*(k-1)+l;
        Mod = NoisyModel(noise,noise_param,n,Mod0,d_alpha*(l-rand(1)));
        Pfail(ind) = 1-QuizzingProb(X,a_x,Mod)*p';
        % calculate the infidelity
        [inFid.max(ind),inFid.avg(ind),inFid.rho(ind),inFid.meas(ind),inFid.type(ind)] = Infidelity(Mod0,Mod,n,r0*Pfail(ind),n_opt);
    end
end

% calculate the ratios
eps = 10^(-6); % tolerance region around 0
r.max = max(max([inFid.max(Pfail>eps);inFid.rho(Pfail>eps);inFid.meas(Pfail>eps)],[],1)./Pfail(Pfail>eps));
r.min = min(max([inFid.max(Pfail>eps);inFid.rho(Pfail>eps);inFid.meas(Pfail>eps)],[],1)./Pfail(Pfail>eps));

end

%% ---------- Quizzing probability ----------

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


%% ---------- Target sequences ----------

function X = Sequences(n)
% SEQUENCES enumerates the sequences required for QSQ
% X is a cell array specifying the sequences
% e.g. X{1} = [1 2 1] corresponds to the sequence s_1 s_2 s_1
if n~=2
    error(['This program needs to be modified to work for ',num2str(n),' qubits'])
end

% enumeration of the gates:
% 1 - S_1
% 2 - S_2
% 3 - CX
% 4 - H_1

% first, sequences for S_2 model
X{1} = [];
X{2} = 1;
X{3} = [1 1];
X{4} = [1 1 1];
X{5} = [1 1 1 1];
X{6} = [2 2 1];
X{7} = [2 2 1 1];
X{8} = [2 2 1 1 1];
X{9} = [2 2 1 1 1 1];
X{10} = 2;
X{11} = [2 2];
X{12} = [2 2 2];
X{13} = [2 2 2 2];
X{14} = [1 1 2];
X{15} = [1 1 2 2];
X{16} = [1 1 2 2 2];
X{17} = [1 1 2 2 2 2];
X{18} = [1 2 1 2];
X{19} = [2 1 2 1];
% now CX
X{20} = 3;
X{21} = [1 1 3];
X{22} = [2 2 3];
X{23} = [1 1 2 2 3];
X{24} = [2 3 2];
X{25} = [1 3 1];
X{26} = [1 2 2 3 1];
% now H_1
X{27} = 4;
X{28} = [1 1 4];
X{29} = [4 4];
X{30} = [1 4 1];
X{31} = [1 1 1 4 1];
X{32} = [4 1 4];
X{33} = [2 2 4];
X{34} = [2 2 1 1 4];
X{35} = [2 2 4 4];
X{36} = [2 2 1 4 1];
X{37} = [2 2 1 1 1 4 1];
X{38} = [2 2 4 1 4];
X{39} = [1 2 4 1 2];
X{40} = [1 2 4 4 1 2];
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
% TARGETMODEL generates the target S_n model for n qubits
% Here as in the proof of S_n model in the Appendix the model is given 
% in the computational basis, i.e., the "S" gate is actually Rx(pi/2)
% gates:
Mod0.gates = cell(1,4);
H = [1, 1; 1, -1]/sqrt(2);
S = H*[1, 0; 0, 1i]*H'; % S gate
Mod0.gates{1} = kron(S,eye(2));
Mod0.gates{2} = kron(eye(2),S);
Mod0.gates{3} = diag([1 1 1 -1]); % CX
Mod0.gates{4} = kron(H,eye(2)); % H
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

function [inFid_max,inFid_avg,inFid_rho,inFid_meas,inFid_type] = Infidelity(Mod0,Mod,n,inFid0,n_opt)
% INFIDELITY calculates the infidelity between the channels, states and the
% measurement
% n_opt -- number of points for optimization

flag = 0;
[fmin,inFid_max,inFid_avg,inFid_rho,inFid_meas,inFid_type] = funMinPhases(zeros(1,2^n-1),Mod,Mod0,n);
if n_opt>0 && inFid_avg>inFid0
    fun = @(x)funMinPhases(x,Mod,Mod0,n);
    options = optimoptions('fminunc','Display','off');
    for k=1:n_opt
        [xk,fk] = fminunc(fun,rand(1,2^n-1),options);
        if fk<fmin
            x = xk;
            fmin = fk;
            flag = 1;
        end
        if fmin<inFid0
            break
        end
    end
    if flag==1
        [~,inFid_max,inFid_avg,inFid_rho,inFid_meas,inFid_type] = funMinPhases(x,Mod,Mod0,n);
    end
end

end

%% ---------- Objective function Infidelity ----------

function [fmin,inFid_max,inFid_avg,inFid_rho,inFid_meas,inFid_type] = funMinPhases(x,Mod,Mod0,n)
% FUNMINPHASES is the objective function for the minimization in Infidelity
d = 2^n;
% construct the gauge
V = eye(d); % possibly add a general gauge
U = V*diag([1,exp(1i*2*pi*x)]); % diagonal in the computational basis

% calculate the infidelities
inFid_rho = real(1-trace(Mod.rho*U*Mod0.rho*U'));
inFid_gates = zeros(1,n);
inFid_gates_c = zeros(1,n);
for k=1:numel(Mod.gates)
    % construct the Choi of the ideal gates
    inFid_gates(k) = (1-real(trace(CJ(Mod.gates{k},d)*CJ(U*Mod0.gates{k}*U',d))))*d/(d+1);
    inFid_gates_c(k) = (1-real(trace(CJ(Mod.gates{k},d)*CJ(U*conj(Mod0.gates{k})*U',d))))*d/(d+1);
end
inFid_meas = 0; % ! this only works for the readout noise and V = eye(d) !
for k=1:d
    inFid_meas = inFid_meas+abs(trace(Mod0.rho*(Mod.M{k}-U*Mod0.M{k}*U')))/2; % total var distance
end
inFid_max = min(max(inFid_gates),max(inFid_gates_c));
inFid_avg = min(sum(inFid_gates),sum(inFid_gates_c))/numel(Mod.gates);
[~,inFid_type] = max([inFid_avg,inFid_rho,inFid_meas]);

fmin = inFid_avg; % change here what the objective function of the optimization is
end



%% ----------- Noise model ---------
function noise = Noise(n,n_gates)
% NOISE generates the noise model
d = 2^n;
% global noise for gates
noise.U = cell(1,n_gates);
noise.D = cell(1,n_gates);
noise.eig = cell(1,n_gates);
for k=1:n_gates
    noise.U{k} = RandomUnitary(d);
    noise.D{k} = rand(1,d-1);
    noise.eig{k} = rand(1,d-1);
end
noise.distr = RandDistr(n_gates);
% depolarizing noise
noise.dep = rand(1,n_gates);
% statistical noise for the measurement
noise.M = rand(1,n);
end


%% ---------- Noisy model ----------

function Mod = NoisyModel(noise,noise_param,n,Mod0,alpha)
% NOISYMODEL generates a random noise added to state, channels and measurement
%   noise -- Noise model
%   noise_param -- parameters controlling the noise strength between
%   different elements
%   alpha -- general noise strength parameter
%
d = 2^n;
% state (white noise)
Mod.rho = (1-noise_param.rho*alpha)*Mod0.rho+noise_param.rho*alpha*eye(d)/d;
% measurement (readout noise)
v = i2s(2*ones(1,n),1:d);
if noise_param.M~=0  
    e = noise_param.M*alpha*noise.M; % level of noise for each qubit
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
else
    Mod.M = Mod0.M;
end
% gates
h = [1, 1; 1, -1]/sqrt(2);
H = h;
for k=2:n
    H = kron(H,h);
end
Mod.gates = cell(1,numel(Mod0.gates));
% coherent noise
for k=1:numel(Mod0.gates)
    if ismember(k,noise_param.list)
        D = diag(exp(1i*2*pi*[0,noise_param.U*alpha*noise.D{k}*noise.distr(k)]));
        U = noise.U{k}*D*noise.U{k}'; % misalignment unitary
        % now over-rotation, given as diagonal in the respective eigenbasis
        if ismember(k,[1,2]) % if it's one of the S gates
            V = H*diag(exp(1i*2*pi*[0,noise_param.eig*alpha*noise.eig{k}*noise.distr(k)]))*H';
        elseif k==3 % if it's CX
            V = diag(exp(1i*2*pi*[0,noise_param.eig*alpha*noise.eig{k}*noise.distr(k)]));
        else % if it's H
            K = kron(sqrt(ones(2)/2+[1 -1; -1 1]/sqrt(8)).*[1 1; 1 -1],eye(2));
            V = K*diag(exp(1i*2*pi*[0,noise_param.eig*alpha*noise.eig{k}*noise.distr(k)]))*K';
        end

        Mod.gates{k} = U*V*Mod0.gates{k}*U';
    else
        Mod.gates{k} = Mod0.gates{k};
    end
end

if noise_param.dep~=0
    e_dep = noise.dep*alpha*noise_param.dep.*noise.distr;
    E = eye(d);
    for k=noise_param.list
        Mod.gates{k}(:,:,1) = sqrt(1-e_dep(k))*Mod.gates{k}(:,:,1);
        for i=1:d
            for j=1:d
                Mod.gates{k}(:,:,(i-1)*d+j+1) = sqrt(e_dep(k)/d)*E(:,i)*E(j,:);
            end
        end
    end
end
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
% DISPLAYPROGRESS is used to display progress in the percentage while k
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

%% ---------- Random Unitary -------------
function U = RandomUnitary(d)
% RANDOMUNITARY constructs a random unitary from U(d)
% Borrowed from http://www.qetlab.com/RandomUnitary

% construct the Ginibre ensemble
gin = randn(d)+1i*rand(d);

% QR decomposition of the Ginibre ensemble
[Q,R] = qr(gin);

% compute U from the QR decomposition
R = sign(diag(R));
R(R==0) = 1; % protect against potentially zero diagonal entries
U = bsxfun(@times,Q,R.'); % much faster than the naive U = Q*diag(R)
end

%% ----------- Random distribution -------
function p = RandDistr(d)
% RANDDISTR generates a random distribution of a d-valued variable
p = [sort(rand(1,d-1)),1];
p(2:end) = p(2:end)-p(1:end-1);
p = p(randperm(d));
end
