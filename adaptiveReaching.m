function simout = adaptiveReaching(loads,gamma,perturbation,tvia,buildup,simin)

% SIMOUT = ADAPTIVEREACHING(LOADS,GAMMA,PERTURBATION,TVIA,BUILDUP,SIMIN)
% 
%   LOADS: 1x2 vector of Force field parameters, LOADS(1): true force field
%          parameter, LOADS(2): estimated force field parameter if SIMIN is
%          empty.
%   GAMMA: Online learning rate
%   PERTURBATION: Amplitude of externally applied step load (not used)
%   TVIA: Time for reaching the via point (see previous study, eNeuro, 2020)
%   BUILDUP: Degree of polynomial buildup for the cost matrices. 
%   SIMIN: Input structure, empty: defaut values, non empty: takes the
%          values of a previous simulation run as input parameters
%
%   SIMOUT: Output structure with state variables, control vector,
%           estimated model matrices.


% script_adaptiveControl
Gx = .1;
Gy = Gx;
m = 1;
tau = 0.1;
L = loads(1); 
lambda = 0;
k = 0;


A = [0 0 1 0 0 0;0 0 0 1 0 0;-k/m 0 -Gx/m L/m m^-1 0;...
    0 -k/m -0/m -Gy/m 0 m^-1;0 0 0 0 -tau^-1 lambda/tau;0 0 0 0 lambda/tau -tau^-1];

ANull = [0 0 1 0 0 0;0 0 0 1 0 0;-k/m 0 -Gx/m 0/m m^-1 0;...
    0 -k/m -0/m -Gy/m 0 m^-1;0 0 0 0 -tau^-1 lambda/tau;0 0 0 0 lambda/tau -tau^-1];

AClamp = [0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0/m 0 0;...
    0 -k/m -0/m -Gy/m 0 m^-1;0 0 0 0 0 0;0 0 0 0 lambda/tau -tau^-1];

B = [0 0;0 0;0 0;0 0;tau^-1 lambda/tau; lambda/tau tau^-1];
%A

simdata.A = A;
simdata.ANull = ANull;
simdata.AClamp = AClamp;
simdata.B = B;

% AParam = -.02*eye(size(A));
xp = mvnrnd(zeros(size(A)),eye(6));
AParam = -eye(6);

% Initialize A matrix
Aest = A;

if isempty(simin)
    Aest(3,4) = loads(2)/m;
    simdata.Clamp = 0;
else
    Aest(3:4,1:4) = simin.AestCont(3:4,1:4,end);
    simdata.Clamp = simin.Clamp;
end


% Model Parameters
simdata.AestCont = Aest;
simdata.BestCont = B;
simdata.AParamCont = AParam;
simdata.time = .4; 
simdata.tvia = tvia;
simdata.stab = .01;
simdata.delta= .01;
simdata.m = m;

% Cost and Learning Rate
simdata.alpha = [200 200 10 10 0 0];
if ~isempty(buildup)
    simdata.buildup = buildup; % Use n-order polynomial builup during movement- 3 good
else
    simdata.buildup = 3;
end
    
simdata.r = 10^-5;
simdata.p = 1;
simdata.gamma = gamma;

% Perturbation Flow
simdata.pert = perturbation;

xinit = [0 0 0 0 0 0]';
xfinal = [0 .15 0 0 0 0]'; 
xvia = [0 .12 0 0 0 0]'; % Previous simulation 
simout = adaptiveLQG(xinit,xfinal,xvia,simdata);

