%% script_AC
close all

% The script calls adaptiveReaching.m, and adapts the online learning rate
% on a trial-by-trial basis. Two simulations are compared: one without
% online learning, one with online learning. The function calls nlinfit and
% nlparci.m (Mathworks) to estimate the exponential fits performed on the
% path length and initial angles of each simulated trajectory. The script
% cass the custom functions expfit.m and expfitdual.m for fitting the
% exponential models with one or two time constants. 

nsimu = 75;
simout = adaptiveReaching([0 0],[0 0],0,0,25,[]);
simout_TBT = simout;
nStep = size(simout.state,2);

allx_TBT = zeros(nsimu,nStep);
ally_TBT = zeros(nsimu,nStep);
allx_Online = zeros(nsimu,nStep);
ally_Online = zeros(nsimu,nStep);

alldx_TBT = zeros(nsimu,nStep);
alldy_TBT = zeros(nsimu,nStep);
alldx_Online = zeros(nsimu,nStep);
alldy_Online = zeros(nsimu,nStep);

pathLength = zeros(nsimu,2);
angle = zeros(nsimu,2);
 
% Good simulations
gamma1 = .001;
gamma2 = 0;
gamma1Max = 0;
gamma1Rate  = .02;
gamma2Max = .2;
gamma2Rate = 0.4;

for i = 1:nsimu

    simout = adaptiveReaching([13 0],[gamma1 0],0,0,25,simout);
    
    simout_TBT = adaptiveReaching([13 0],[0 0],0,0,25,simout);
    simout_Online = adaptiveReaching([13 0],[gamma2 0],0,0,25,simout);
    
    allx_TBT(i,:) = simout_TBT.state(1,:);
    ally_TBT(i,:) = simout_TBT.state(2,:);
    
    alldx_TBT(i,:) = simout_TBT.state(3,:);
    alldy_TBT(i,:) = simout_TBT.state(4,:);
    
    allx_Online(i,:) = simout_Online.state(1,:);
    ally_Online(i,:) = simout_Online.state(2,:);
    
    alldx_Online(i,:) = simout_Online.state(3,:);
    alldy_Online(i,:) = simout_Online.state(4,:);
    
    %calculate pathLength and angle: Fix theta
    ind = find(simout_TBT.state(2,:)>.05,1,'first');
    angle(i,1) = -(atan2(simout_TBT.state(2,ind),simout_TBT.state(1,ind))-pi/2)*180/pi;
    pathLength(i,1) = .01*sum((simout_TBT.state(3,:).^2 + simout_TBT.state(4,:).^2).^0.5);
    
    % calculate pathLength and angle: variable theta
    ind = find(simout.state(2,:)>.05,1,'first');
    angle(i,2) = -(atan2(simout_Online.state(2,ind),simout_Online.state(1,ind))-pi/2)*180/pi;
    pathLength(i,2) = .01*sum((simout_Online.state(3,:).^2 + simout_Online.state(4,:).^2).^0.5);
    
    % update online learning rate
    gamma2 = gamma2 + gamma2Rate*(gamma2Max-gamma2);
    gamma1 = gamma1 + gamma1Rate*(gamma1Max-gamma1);

end


figure
subplot(131)
ht1 = plot(allx_TBT(end,:),ally_TBT(end,:),'k:','linewidth',2); hold on
ht2 = plot(allx_Online(end,:),ally_Online(end,:),'k','linewidth',2);
axis square; axis([-.1 .1 -.02 .18]);
title('Trajectories'); xlabel('x [m]'); ylabel('y [m]'); 

subplot(132)
h1 = plot(angle(:,1),'ko'); hold on;
h2 = plot(angle(:,2),'k.','MarkerSize',20); hold on;
xplot = [1:75];
axis square;
xlabel('Trial Number'); ylabel('Initial Angle [deg]');

subplot(133), 
plot(pathLength(:,1),'ko'); hold on;
plot(pathLength(:,2),'k.','MarkerSize',20); hold on;

% fit 2 time constants to path length - online learning on
[beta2PL,r2PL,J2PL,sigma2PL] = nlinfit([1:nsimu]',pathLength(:,2),@expfitdual,[1 1 -1 1 -1]/10);
ci2PL = nlparci(beta2PL,r2PL,'covar',sigma2PL);

% fit one time constant to the initial angle, online learning on
[beta1Angle,r1Angle,~,s1Angle] = nlinfit([1:75]',angle(:,2),@expfit,[1 20 -.1/10]);
ci1Angle = nlparci(beta1Angle,r1Angle,'covar',s1Angle);

% fit one time constant to the initial angle, online learning off
[beta1AngleOFF,r1AngleOFF,~,s1AngleOFF] = nlinfit([1:75]',angle(:,1),@expfit,[1 20 -.1/10]);
ci1AngleOFF = nlparci(beta1AngleOFF,r1AngleOFF,'covar',s1AngleOFF);

% fit one time constant to the path length, online learning off
[beta1PLOFF,r1PLOFF,~,s1PLOFF] = nlinfit([1:75]',pathLength(:,1),@expfit,[1 .5 -1/10]);
ci1PLOFF = nlparci(beta1PLOFF,r1PLOFF,'covar',s1PLOFF);


subplot(132)
h3 = plot(xplot,expfit(beta1AngleOFF,xplot'),'b','Linewidth',2);
plot(xplot,expfit(beta1Angle,xplot'),'b','Linewidth',2);
axis square;

subplot(133)
hs = plot(xplot,expfit(beta1PLOFF,xplot'),'b','Linewidth',2);
hd = plot(xplot,expfitdual(beta2PL,xplot'),'r','Linewidth',2);
axis square;
xlabel('Trial Number'); ylabel('Path Length [m]');

legend([ht1,ht2,h1,h2,h3,hd],'Offline','Offlin and Online','Offline',...
    'Offline and Online','Single Rate','Dual Rate','Location','northOutside');
