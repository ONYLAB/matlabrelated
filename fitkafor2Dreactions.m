function fitkafor2Dreactions()
% Concentrations are in units of molecules per nm^{2}
% Length scale s is in units of nm
% Data file should be tab delimited in the format
% t({\mu}s)  A(t)(molecules per nm^{2})
% 0          A(0)
% ..         ..
% ..         ..
% ..         ..
% k and D are in units of nm^2/({\mu}s)
% 7/2015 Yogurtcu&Johnson

% Clean out everything
killEmAll;
close all hidden;
close all force;
clc;

global prev;
prev.x = [];
prev.fval = [];

% % initialization
tsample = 100; %used for data file refinement
maxfuneval = 300; %used for optimization
tolfun = 1e-12; %for optimization
klist = logspace(-3,3,7); %k_a list for multiple initial starting points
timelowerbndmultip = 1e2; %used for asymptotic fit timelowerbndmultip*s^2/D is the time lowerbound for curve fitting long t scale
timeupperbndmultip = 1e-3; %used for asymptotic fit timeupperbndmultip*s^2/D is the time upperbound for curve fitting short t scale
lb = [1e-3; 5e-3;  0.9;]; %lower bounds k D s
ub = [1e3;  20;   20;]; %upper bounds k D s
lbkb = [1e-3; 5e-3;  0.9; 0]; %lower bounds k D s kb
ubkb = [1e3;  20;   20; 1e2*1e-6]; %upper bounds k D s kb
optveci = [1 1 1 1e-6]; % optveci = [k D s kb]; initial optim vec
nummultivaroptim = 20; %number of multivariable optimization iterations
fixestimates = 0;
A0 = 0;
B0 = 0;
k  = 0;
kb = 0;
D  = 0;
s  = 0;
warning on;

% Construct questdlg with two options
fitchoice = questdlg('Do you have initial estimates for D and sigma?', ...
    'Any estimates?', ...
    'Yes','No','No');
switch fitchoice
    case 'Yes'
        estimates = 1;
        
        % Construct questdlg with two options
        fitchoice2 = questdlg('Would you like to keep these estimates fixed?', ...
            'Fix estimates?', ...
            'Yes','No','No');
        [choice,B0,D,s,rxntype] = constructquestdlgALL();
        
        switch fitchoice2
            case 'Yes'
                fixestimates = 1;
                optveci = [1 1e-6]; % optveci = [k kb]; new initial optim vec created
            case 'No'
                fixestimates = 0;
                optveci = [1 D s 1e-6]; % optveci = [k D s kb]; initial optim vec changed
                lbkb = [1e-3; max(5e-3,D/3);  max(0.9,s/3); 0]; %lower bounds k D s kb
                ubkb = [1e3;  min(20,D*3);   min(20,s*3); 1e2*1e-6]; %upper bounds k D s kb
        end
        
    case 'No'
        estimates = 0;
        [choice,B0,rxntype] = constructquestdlg();
end


% Read time data file
[A0,timedat] = readtimedata(tsample);
disp([choice ' optimization starting....'])
disp('________________');

% Now let's do multivariable optimization
h = waitbar(0,'Please wait...','Name','Variable fitting...');
warning off
disp('Optimization Fit Result!');
optionsbase = optimoptions('fmincon','Display','off','TolFun', tolfun,'Algorithm','sqp'); %choice of sqp is not arbitrary, it is for to stay away from possible NaNs/Infs on the costfunction

if fixestimates == 1
    
    if rxntype<3
        % Irreversible Simulations
        
        options = optimoptions(optionsbase,'MaxFunEvals',ceil(maxfuneval/3),'OutputFcn',@outfun1D); % run interior-point algorithm
        [fits{1},fval(1),~,output] = fmincon(@(optvec) costfunctionvectorizedSMALL(optvec,timedat,rxntype,A0,B0,D,s,timelowerbndmultip,timeupperbndmultip),optveci(1),[],[],[],[],lb(1),ub(1),[],options);
        frstordopt(1) = output.firstorderopt;
        % multi-start optim test
        for i = 2:nummultivaroptim
            optveci = abs(fits{i-1}.*(1.4*rand(1,length(optveci(1)))+0.3));
            [fits{i},fval(i),~,output] = fmincon(@(optvec) costfunctionvectorizedSMALL(optvec,timedat,rxntype,A0,B0,D,s,timelowerbndmultip,timeupperbndmultip),optveci,[],[],[],[],lb(1),ub(1),[],options);
            waitbar(i / length(klist))
            frstordopt(i) = output.firstorderopt;
        end
        [~,ind] = min(fval);
        fitall = fits{ind};
        fitall(4) = 0;
        fitfirstorderoptimality = frstordopt(ind);
    else
        % Reversible Simulations
        
        options = optimoptions(optionsbase,'MaxFunEvals',ceil(maxfuneval/3*2),'OutputFcn',@outfun2D); % run interior-point algorithm
        [fits{1},fval(1),~,output] = fmincon(@(optvec) costfunctionvectorizedSMALL(optvec,timedat,rxntype,A0,B0,D,s,timelowerbndmultip,timeupperbndmultip),optveci,[],[],[],[],[lbkb(1) lbkb(4)],[ubkb(1) ubkb(4)],[],options);
        frstordopt(1) = output.firstorderopt;
        
        % multi-start optim test
        for i = 2:ceil(nummultivaroptim/3*2)
            optveci = abs(fits{i-1}.*(1.4*rand(1,length(optveci))+0.3));
            [fits{i},fval(i),~,output] = fmincon(@(optvec) costfunctionvectorizedSMALL(optvec,timedat,rxntype,A0,B0,D,s,timelowerbndmultip,timeupperbndmultip),optveci,[],[],[],[],[lbkb(1) lbkb(4)],[ubkb(1) ubkb(4)],[],options);
            waitbar(i / length(klist))
            frstordopt(i) = output.firstorderopt;
        end
        [~,ind] = min(fval);
        fitall = fits{ind};
        disp(['Fitted k_b ' num2str(fitall(2)*1e6) ' s^{-1}'])
        fitall(4) = fitall(2);
        fitfirstorderoptimality = frstordopt(ind);
        
        erplot2D(prev);
        
    end
    disp(['Fitted k_a ' num2str(fitall(1)) ' nm^2/{\mu}s'])
    disp('________________');
    close(h)
    
    fitall(2) = D;
    fitall(3) = s;
    
else
    
    if rxntype<3
        % Irreversible Simulations
        
        options = optimoptions(optionsbase,'MaxFunEvals',maxfuneval); % run interior-point algorithm
        [fits{1},fval(1),~,output] = fmincon(@(optvec) costfunctionvectorized(optvec,timedat,rxntype,A0,B0,timelowerbndmultip,timeupperbndmultip),optveci(1:3),[],[],[],[],lb,ub,[],options);
        frstordopt(1) = output.firstorderopt;
        
        % multi-start optim test
        for i = 2:nummultivaroptim
            optveci = abs(fits{i-1}.*(1.4*rand(1,length(optveci(1:3)))+0.3));
            [fits{i},fval(i),~,output] = fmincon(@(optvec) costfunctionvectorized(optvec,timedat,rxntype,A0,B0,timelowerbndmultip,timeupperbndmultip),optveci,[],[],[],[],lb,ub,[],options);
            waitbar(i / length(klist))
            frstordopt(i) = output.firstorderopt;
        end
        [~,ind] = min(fval);
        fitall = fits{ind};
        fitall(4) = 0;
        fitfirstorderoptimality = frstordopt(ind);
                
    else
        % Reversible Simulations
        
        options = optimoptions(optionsbase,'MaxFunEvals',ceil(maxfuneval/3*4)); % run interior-point algorithm
        [fits{1},fval(1),~,output] = fmincon(@(optvec) costfunctionvectorized(optvec,timedat,rxntype,A0,B0,timelowerbndmultip,timeupperbndmultip),optveci,[],[],[],[],lbkb,ubkb,[],options);
        frstordopt(1) = output.firstorderopt;
        
        % multi-start optim test
        for i = 2:ceil(nummultivaroptim/3*4)
            optveci = abs(fits{i-1}.*(1.4*rand(1,length(optveci))+0.3));
            [fits{i},fval(i),~,output] = fmincon(@(optvec) costfunctionvectorized(optvec,timedat,rxntype,A0,B0,timelowerbndmultip,timeupperbndmultip),optveci,[],[],[],[],lbkb,ubkb,[],options);
            waitbar(i / length(klist))
            frstordopt(i) = output.firstorderopt;
        end
        [~,ind] = min(fval);
        fitall = fits{ind};
        disp(['Fitted k_b ' num2str(fitall(4)*1e6) ' s^{-1}'])
        fitfirstorderoptimality = frstordopt(ind);
                
    end
    disp(['Fitted k_a ' num2str(fitall(1)) ' nm^2/{\mu}s'])
    disp(['Fitted D ' num2str(fitall(2)) ' nm^2/{\mu}s'])
    disp(['Fitted s ' num2str(fitall(3)) ' nm'])
    disp('________________');
    close(h)
    
end

% plot the curves and k(t)
[R2,corrcoeff,MSE,MAE,MAXE] = plotall(fitall,A0,B0,timedat,rxntype,timelowerbndmultip,timeupperbndmultip);

%show relative errors
writestats(fitfirstorderoptimality,R2,corrcoeff,MSE,MAE,MAXE);
% reportrelativeerror(fitall)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % END  MAIN % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function stop = outfun1D(x, optimValues, state)
stop = false;
hold on;
plot(x,optimValues.fval,'.');
ylabel('Cumulative Residual Error');
xlabel('k_a [nm^2/{\mu}s]')
box on
drawnow

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function stop = outfun2D(x, optimValues, state)
global prev;
stop = false;
prev.fval = [prev.fval; optimValues.fval];
prev.x = [prev.x; x];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function erplot2D(prev)
xs = prev.x;
scatter3(xs(:,1),xs(:,2)*1e6,prev.fval,'*')
title('Cumulative Residual Error');
xlabel('k_a [nm^2/{\mu}s]');
ylabel('k_b [s^{-1}]');
box on
axis square

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function kt = givekt(t,ktvec,tvec)
% find kt interpolated

kt = interp1(tvec,ktvec,t);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime)
%creates sampled k(t)

numsample = 20;

ttemp = logspace(log10(time(2)),log10(time(end)),numsample);
ttemp(1) = time(2);
ttemp(end) = time(end);
tvec = [0;ttemp(:)];
ktvec = createkt(tvec,D,k,s,scalarlongtime,scalarshorttime);

function kt = createkt(ts,D,k,s,scalarlongtime,scalarshorttime)
% gives k(t), uses approximate k(t) and shot-time and long-time asymptotics

kt = zeros(length(ts),1);
cutofflongtime = scalarlongtime*s*s/D;
cutoffshorttime = scalarshorttime*s*s/D;

for i = 1:length(ts)
    t = ts(i);
    
    if t<=cutoffshorttime
        
        % This is the short time asymptotic
        T = D.*t./s./s;
        
        if t == 0
            kt(i) = k;
        else
            if isinf(k)
                kt(i) = 2*pi*D*(1./(sqrt(pi*T)) + 1/2 - 1/4*sqrt(T/pi) + 1/8*T - 100/384*T.*sqrt(T/pi));
            else
                kt(i) = 2*pi*D*(1/2 - k/(2*pi*D)*sqrt(T/pi) - k/(4*pi*D)*T.*sqrt(T/pi) + k/(2*pi*D)*erfcx(k/(2*pi*D)*sqrt(T)).*(1 - pi*D/k + k/(2*pi*D)*T + 1/2*(k/(2*pi*D)*T).^2)  );
            end
        end
        
    elseif t>=cutofflongtime
        
        % This is the long time asymptotic
        
        B = 4*pi*D/k;
        g = -psi(1);
        
        C = D*t/s/s;
        % Eq longtassymptotic integral of Yogurtcu and Johnson 2015
        kt(i) = 4*pi*D*((1./(log(4*C)-2*g+B))-g*(1./(log(4*C)-2*g+B).^2)-1.331*(2./(log(4*C)-2*g+B).^3)+0.25*(2./(log(4*C)-2*g+B).^4));
        
    else
        kt(i) = ktnum(t,s,D,k);
        
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [A0,timedat] = readtimedata(tsample)
% Read time data file

[filename, pathname] = ...
    uigetfile({'*.*'},'Please select your tab delimited data file for A(t)');
timedat = dlmread([pathname filename]);
disp(['Finished reading ' [pathname filename] '.....']);
disp('Analyzing and refining data file..');

newt = logspace(log10(timedat(2,1)),log10(timedat(end,1)),tsample);
newt = [timedat(1,1); newt(:)];
newt(end)=timedat(end,1);

newA = interp1(timedat(:,1),timedat(:,2),newt);
timedat = [newt newA(:)];

% initial concentration
A0 = timedat(1,2);

if sum(find(timedat(:,2)==0))
    disp('Zero(s) detected in data file for A(t). Curve fitting might not be perfect :(')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [choice,B0,rxntype] = constructquestdlg()
% Construct questdlg with two options

B0 = 0;
choice = questdlg('What is the type of your reaction?', ...
    'Reaction Type', ...
    'Reversible','Irreversible','Reversible');
switch choice
    case 'Irreversible'
        % Construct another questdlg with two options
        choice = questdlg('Which irreversible reaction is it?', ...
            'Reaction Type', ...
            'A+A => 0','A+B => 0','A+A => 0');
        
        switch choice
            case 'A+A => 0'
                rxntype = 1;
                
            case 'A+B => 0'
                rxntype = 2;
                options.Interpreter='tex';
                parms = inputdlg({'Initial Concentration of B (molecules/nm^{2})'},...
                    'Params', [1 50;],{'0'},options);
                B0 = str2num(parms{1});
        end
        
    case 'Reversible'
        % Construct another questdlg with two options
        choice = questdlg('Which reversible reaction is it?', ...
            'Reaction Type', ...
            'A+A <=> C','A+B <=> C','A+A <=> C');
        
        % Handle response
        switch choice
            case 'A+A <=> C'
                rxntype = 3;
                
            case 'A+B <=> C'
                rxntype = 4;
                options.Interpreter='tex';
                parms = inputdlg({'Initial Concentration of B (molecules/nm^{2})'},...
                    'Params', [1 50;],{'0'},options);
                B0 = str2num(parms{1});
        end
        
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [choice,B0,D,s,rxntype] = constructquestdlgALL()
% Construct questdlg with two options

B0 = 0;
D = 0;
s = 0;

choice = questdlg('What is the type of your reaction?', ...
    'Reaction Type', ...
    'Reversible','Irreversible','Reversible');
switch choice
    case 'Irreversible'
        % Construct another questdlg with two options
        choice = questdlg('Which irreversible reaction is it?', ...
            'Reaction Type', ...
            'A+A => 0','A+B => 0','A+A => 0');
        switch choice
            case 'A+A => 0'
                rxntype = 1;
            case 'A+B => 0'
                rxntype = 2;
        end
        
    case 'Reversible'
        % Construct another questdlg with two options
        choice = questdlg('Which reversible reaction is it?', ...
            'Reaction Type', ...
            'A+A <=> C','A+B <=> C','A+A <=> C');
        % Handle response
        switch choice
            case 'A+A <=> C'
                rxntype = 3;
            case 'A+B <=> C'
                rxntype = 4;
        end
        
end

options.Interpreter='tex';

if rxntype == 2 || rxntype == 4
    
    parms = inputdlg({'Estimated Total Diffusion Coefficient (nm^2/{\mu}s)','Estimated Reaction Length Scale, {\sigma}, (nm)','Initial Concentration of B (molecules/nm^{2})'},...
        'Params', [1 50; 1 50; 1 50],{'0','0','0'},options);
    D = str2num(parms{1});
    s = str2num(parms{2});
    B0 = str2num(parms{3});
    
else
    
    parms = inputdlg({'Estimated Total Diffusion Coefficient (nm^2/{\mu}s)','Estimated Reaction Length Scale, {\sigma}, (nm)'},...
        'Params', [1 50; 1 50],{'0','0'},options);
    D = str2num(parms{1});
    s = str2num(parms{2});
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function At = rxntype1(A0,k,s,D,time,scalarshorttime,scalarlongtime)
%%% FOR RXN TYPE: A+A->0

At = zeros(length(time),1);
[ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime);

for i = 1:length(time)
    
    ts = time(i);
    if ts == 0
        ktint = 0;
    else
        ktint = integral(@(t) givekt(t,ktvec,tvec),0,ts);
    end
    
    At(i) = A0/(1+A0*ktint);
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function At = rxntype2(A0,B0,k,s,D,time,scalarshorttime,scalarlongtime)
%%% FOR RXN TYPE: A+B->0

At = zeros(length(time),1);
[ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime);

if A0~=B0
    
    for i = 1:length(time)
        
        ts = time(i);
        if ts == 0
            ktint = 0;
        else
            ktint = integral(@(t) givekt(t,ktvec,tvec),0,ts);
        end
        At(i) = A0*(B0-A0)/(B0*exp((B0-A0)*ktint)-A0);
    end
    
else
    
    [~,res] = ode23s(@(t,y) ODEforrxntype2(t,y,k,ktvec,tvec),time,[A0 B0]);
    At = res(:,1);
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function At = rxntype3(A0,k,kb,s,D,time,scalarshorttime,scalarlongtime)
%%% FOR RXN TYPE: A+A<->C

[ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime);
[~,res] = ode23s(@(t,y) ODEforrxntype3(t,y,k,kb,ktvec,tvec),time,[A0 0]);
At = res(:,1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function At = rxntype4(A0,B0,k,kb,s,D,time,scalarshorttime,scalarlongtime)
%%% FOR RXN TYPE: A+B<->0

[ktvec,tvec] = ktassister(time,D,k,s,scalarlongtime,scalarshorttime);
[~,res] = ode23s(@(t,y) ODEforrxntype4(t,y,k,kb,ktvec,tvec),time,[A0 B0 0]);
At = res(:,1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function cost = costfunctionvectorized(optvec,timedat,rxntype,A0,B0,scalarlongtime,scalarshorttime)

cost = 0;
small = 1e-10; %to prevent possible divby0
k=optvec(1);
D=optvec(2);
s=optvec(3);

if rxntype == 1
    cost = sum(abs(timedat(:,2)-rxntype1(A0,k,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
elseif rxntype == 2
    cost = sum(abs(timedat(:,2)-rxntype2(A0,B0,k,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
elseif rxntype == 3
    kb=optvec(4);
    cost = sum(abs(timedat(:,2)-rxntype3(A0,k,kb,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
elseif rxntype == 4
    kb=optvec(4);
    cost = sum(abs(timedat(:,2)-rxntype4(A0,B0,k,kb,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function cost = costfunctionvectorizedSMALL(optvec,timedat,rxntype,A0,B0,D,s,scalarlongtime,scalarshorttime)

cost = 0;
small = 1e-10; %to prevent possible divby0
k=optvec(1);

if rxntype == 1
    cost = sum(abs(timedat(:,2)-rxntype1(A0,k,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
elseif rxntype == 2
    cost = sum(abs(timedat(:,2)-rxntype2(A0,B0,k,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
elseif rxntype == 3
    kb=optvec(2);
    cost = sum(abs(timedat(:,2)-rxntype3(A0,k,kb,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
elseif rxntype == 4
    kb=optvec(2);
    cost = sum(abs(timedat(:,2)-rxntype4(A0,B0,k,kb,s,D,timedat(:,1),scalarshorttime,scalarlongtime))./(small+timedat(:,2)));
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function dy = ODEforrxntype2(t,y,ktvec,tvec)

dy = zeros(2,1);    % a column vector
thekt = givekt(t,ktvec,tvec);
dy(1) = -thekt * y(1) * y(2);
dy(2) = -thekt * y(1) * y(2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function dy = ODEforrxntype3(t,y,k,kb,ktvec,tvec)

dy = zeros(2,1);    % a column vector
thekt = givekt(t,ktvec,tvec);
dy(1) = -thekt * y(1) * y(1) + 2*kb/k * thekt * y(2);
dy(2) =  1/2*thekt * y(1) * y(1) - kb/k * thekt * y(2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function dy = ODEforrxntype4(t,y,k,kb,ktvec,tvec)

dy = zeros(3,1);    % a column vector
thekt = givekt(t,ktvec,tvec);
dy(1) = -thekt * y(1) * y(2) + kb/k * thekt * y(3);
dy(2) = -thekt * y(1) * y(2) + kb/k * thekt * y(3);
dy(3) =  thekt * y(1) * y(2) - kb/k * thekt * y(3);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [R2,corrcoeff,MSE,MAE,MAXE] = plotall(fitall,A0,B0,timedat,rxntype,timelowerbndmultip,timeupperbndmultip)
% Plot all curves... A(t), fitted A(t) and the predicted k(t)

kfit = fitall(1);
Dfit = fitall(2);
sfit = fitall(3);
kbfit = fitall(4);

whatsmyregime(kfit,Dfit);

figure
% plot fit/optim results multivariable fit
semilogx(timedat(:,1), timedat(:,2), 'r','LineWidth',2)
hold on

[AT,~]=plotfunction(kfit,kbfit,A0,B0,sfit,Dfit,timedat,rxntype,timelowerbndmultip,timeupperbndmultip);
h = legend('Exp. Data', 'Fit','Location','NorthEast');
set(h,'FontSize',12)
title('Fit result')

times = logspace(-6,5,100);
kt = createkt(times,Dfit,kfit,sfit,timelowerbndmultip,timeupperbndmultip);

figure
semilogx(times, kt,'LineWidth',2)
xlabel('Time ({\mu}s)','FontSize',16);
set(gca,'FontSize',16);
ylabel('k(t) (nm^{2}/{\mu}s)','FontSize',16);
set(gcf,'color','w');

[R2,corrcoeff,MSE,MAE,MAXE] = calcstats(timedat(:,2),AT);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [R2,corrcoeff,MSE,MAE,MAXE] = calcstats(obs,fit)
%calculates fit stats

numdatpts = length(fit);
ybar = mean(obs);
fbar = mean(fit);
difv = obs-ybar;
SStot = dot(difv,difv);
difo = fit-ybar;
SSreg = dot(difo,difo);
difs = fit-obs;
SSres = dot(difs,difs);
R2 = 1-SSres/SStot;

corrcoeff = dot(obs-ybar,fit-fbar) / sqrt(dot(obs-ybar,obs-ybar)) / sqrt(dot(fit-fbar,fit-fbar));

MSE = dot(obs-fit,obs-fit)/numdatpts;
MAE = sum(abs(obs-fit))/numdatpts;
MAXE = max(obs-fit);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function writestats(fitfirstorderoptimality,R2,corrcoeff,MSE,MAE,MAXE)

disp('________________');
disp('Fit Statistics:');
disp(['First-Order Optimality Measure: ' num2str(fitfirstorderoptimality)]);
disp(['R^2 Goodness of Fit: ' num2str(R2)]);
disp(['Correlation Coefficient: ' num2str(corrcoeff)]);
disp(['Maximum Error: ' num2str(MAXE)]);
disp(['Mean Squared Error: ' num2str(MSE)]);
disp(['Mean Absolute Error: ' num2str(MAE)]);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [AT,T] = plotfunction(k,kb,A0,B0,s,D,timedat,rxntype,scalarlongtime,scalarshorttime)
% plot A(t)

AT = zeros(length(timedat),1);
AT(1:length(timedat)) = nan;
T = timedat(:,1);

if rxntype == 1
    
    AT = rxntype1(A0,k,s,D,timedat(:,1),scalarshorttime,scalarlongtime);
    
elseif rxntype == 2
    
    AT = rxntype2(A0,B0,k,s,D,timedat(:,1),scalarshorttime,scalarlongtime);
    
elseif rxntype == 3
    
    AT = rxntype3(A0,k,kb,s,D,timedat(:,1),scalarshorttime,scalarlongtime);
    
elseif rxntype == 4
    
    AT = rxntype4(A0,B0,k,kb,s,D,timedat(:,1),scalarshorttime,scalarlongtime);
    
end

semilogx(T, AT,'--','LineWidth',2);
xlabel('Time ({\mu}s)','FontSize',16);
set(gca,'FontSize',16);
ylabel('A(t) (molecules/nm^{2})','FontSize',16);
set(gcf,'color','w');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function killEmAll
% Close all open dialog boxes

set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
set(0,'ShowHiddenHandles','off');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function whatsmyregime(kfit,Dfit)
temp = kfit/Dfit;
if temp<0.05
    disp('This reaction is in the Rate Limited Regime (ka/D < 0.05).');
    disp('For k(t), less than 10% deviation from k(t=0) expected.');
elseif temp<0.5
    disp('This reaction is in the Weakly Diffusion Influenced Regime (0.05 < ka/D < 0.5).');
    disp('For k(t), less than 50% (more than 10%) deviation from k(t=0) expected.');
elseif temp<20
    disp('This reaction is in the Diffusion Influenced Regime (0.5 < ka/D < 20).');
    disp('For k(t), less than 95% (more than 50%) deviation from k(t=0) expected.');
else
    disp('This reaction is in the Diffusion Controlled Regime (ka/D > 20).');
    disp('For k(t), more than 95% deviation from k(t=0) expected.');
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function f = ktnumtobeintegrated(x,t,s,D,k)

h = 2.0 * pi * s * D;
alp = h * x .* j1(x * s) + k * j0(x * s);
bet = h * x .* y1(x * s) + k * y0(x * s);
tet = alp .* alp + bet .* bet;
P = j0(x * s) .* y1(x * s) - j1(x * s) .* y0(x* s);
T = (j0(x * s) .* bet - y0(x * s) .* alp) ./ tet;

f = T .* P .* (1 - exp(-D * t * x .* x));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function kt = ktnum(ts,s,D,k)

for i = 1:length(ts)
    t = ts(i);
    integra = integral(@(x)ktnumtobeintegrated(x,t,s,D,k),0,Inf,'RelTol',1e-10,'AbsTol',1e-10);
    kt(i) = k*(1-integra*s*k);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = j0(inp)

r = besselj(0,inp);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = j1(inp)

r = besselj(1,inp);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = y0(inp)

r = bessely(0,inp);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function r = y1(inp)

r = bessely(1,inp);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function ktintegrated = ktlongtassym(tmax,D,k,s)
% This is the long time asymptotic integrated
% Not used

B = 4*pi*D/k;
g = -psi(1);

C = D*tmax/s/s;
% Eq longtassymptotic integral of Yogurtcu and Johnson 2015
ktmax = 4*pi*s*s*(C*(1/(log(4*C)-2*g+B)+1/(log(4*C)-2*g+B)^2 +2/(log(4*C)-2*g+B)^3)-g*C*(1/(log(4*C)-2*g+B)^2 +2/(log(4*C)-2*g+B)^3)-1.331/2*C*(2/(log(4*C)-2*g+B)^3));

C = 0;%D*tmin/s/s;
% Eq longtassymptotic integral of Yogurtcu and Johnson 2015
ktmin = 4*pi*s*s*(C*(1/(log(4*C)-2*g+B)+1/(log(4*C)-2*g+B)^2 +2/(log(4*C)-2*g+B)^3)-g*C*(1/(log(4*C)-2*g+B)^2 +2/(log(4*C)-2*g+B)^3)-1.331/2*C*(2/(log(4*C)-2*g+B)^3));

ktintegrated = ktmax-ktmin;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function ktintegrated = ktshorttassym(tm,d,k,s)
% This is the short time asymptotic integrated
% Not used

if tm>0
    ktintegrated = (1/30).*k.^(-3).*pi.^(-3/2).*s.^(-2).*tm.^(-1).*(tm.*(480.*d.^3.*( ...
        (-1)+exp(1).^((1/4).*d.^(-1).*k.^2.*pi.^(-2).*s.^(-2).*tm)).*pi.^( ...
        11/2).*s.^4+3.*k.^3.*(d.*s.^(-2).*tm).^(1/2).*(40.*pi.^2.*s.^4+ ...
        k.^2.*tm.^2)+(-120).*d.^2.*k.*pi.^(7/2).*s.^2.*((-3).*pi.*s.^2+( ...
        -4).*pi.^(1/2).*s.^2.*(d.*s.^(-2).*tm).^(1/2)+exp(1).^((1/4).*d.^( ...
        -1).*k.^2.*pi.^(-2).*s.^(-2).*tm).*(3.*pi.*s.^2+k.*tm))+d.*k.^2.* ...
        pi.*(15.*exp(1).^((1/4).*d.^(-1).*k.^2.*pi.^(-2).*s.^(-2).*tm).* ...
        pi.^(1/2).*(8.*pi.^2.*s.^4+4.*k.*pi.*s.^2.*tm+k.^2.*tm.^2)+(-2).*( ...
        3.*k.^2.*tm.^2.*(d.*s.^(-2).*tm).^(1/2)+60.*pi.^2.*s.^4.*(pi.^( ...
        1/2)+3.*(d.*s.^(-2).*tm).^(1/2))+5.*k.*pi.*s.^2.*tm.*((-3).*pi.^( ...
        1/2)+4.*(d.*s.^(-2).*tm).^(1/2)))))+(-15).*d.^(1/2).*exp(1).^(( ...
        1/4).*d.^(-1).*k.^2.*pi.^(-2).*s.^(-2).*tm).*pi.^(3/2).*s.*tm.^( ...
        1/2).*(d.*s.^(-2).*tm).^(1/2).*((-24).*d.*k.*pi.^3.*s.^4+32.* ...
        d.^2.*pi.^4.*s.^4+4.*k.^3.*pi.*s.^2.*tm+k.^4.*tm.^2+8.*k.^2.* ...
        pi.^2.*s.^2.*(s.^2+(-1).*d.*tm)).*erf((1/2).*d.^(-1/2).*k.*pi.^( ...
        -1).*s.^(-1).*tm.^(1/2)));
else
    ktintegrated = 0;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function reportrelativeerror(fitall)
% Reports relative error NOT USED ANYMORE
k = fitall(1);
D = fitall(2);
s = fitall(3);

load('relerrdata.mat');
datapoints = relerr(:,1:3);
mypoint = [s D k];

diffpts = datapoints - repmat(mypoint,length(datapoints),1);
aa = diffpts * diffpts';
[~,ind]=min(sqrt(diag(aa)));
disp('________________');
disp(['The estimated relative error in k(t) at t=10^-3 sigma^2/D is ' num2str(relerr(ind,4))]);
disp(['The estimated relative error in k(t) at t=10 sigma^2/D is ' num2str(relerr(ind,5))]);
