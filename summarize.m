function summarize()

close all;
[~,~]=unix('rm -rf Summary/');
[~,~]=unix('rm *~');
[~,~]=unix('chmod 777 *');
timefreq = 1000;

% % READ IN PARAMETER FILE
fileinf = dir('*parms*');
fileID = fopen(fileinf.name);
C = textscan(fileID,'%f %*[^\n]');
C = C{:};
fclose(fileID);

fileinf = dir('Correlation*');
cfname = fileinf.name;

Nprotypes = C(1);
nint = C(2);
maxtime = C(3);
Nspecies = C(4);
Nrxn = C(5);
vol = C(6);
Vsphere = C(7);
checkpointiter = C(8);
stepwrite = C(9);
niter = C(10);
delt = C(11);

choice = questdlg('Would you like to run a simulation?', ...
    'Simulation?', ...
    'Yes','No, thank you!','No, thank you!');
% Handle response
switch choice
    case 'Yes'
        fileinf = dir('*parms*');
        f1 = fileinf.name;
        fileinf = dir('concfile*');
        f2 = fileinf.name;
        fileinf = dir('inet*');
        f3 = fileinf.name;
        fileinf = dir('progill*');
        f4 = fileinf.name;
        fileinf = dir('rxn*');
        f5 = fileinf.name;
        fileinf = dir('status*');
        f6 = fileinf.name;
                
        fileID = fopen(f1,'w');
        texts = {'#Npro types',...
            '#n protein interfaces',...
            '#maxtime',...
            '#Nspecies',...
            '#Nrxn',...
            '#volume um^3',...
            '#Vsphere 1=sphere geom, else box',...
            '#checkpoint iterations',...
            '#stepwrite',...
            '#niterations (if repeat)',...
            '#delt for histogramming (sec)'};
        
        dlg_title = 'Simulation Parameters';
        num_lines = 1;
        def = {num2str(C(1)),num2str(C(2)),num2str(C(3)),num2str(C(4)),num2str(C(5)),num2str(C(6)),num2str(C(7)),num2str(C(8)),num2str(C(9)),num2str(C(10)),num2str(C(11))};
        answer = inputdlg(texts,dlg_title,num_lines,def);
        
        for i = 1:11
            
            if isempty(find(i==[3 6]))
                fprintf(fileID,'%d %s\n',str2num(answer{i}),texts{i});
            else
                fprintf(fileID,'%f %s\n',str2num(answer{i}),texts{i});
            end
            
        end
        
        fclose(fileID);       
        
        hwarn = warndlg('Running a new simulation!','!! Warning !!');
        [~,~] = unix(['./ssa_complex_double_repeat.exe ' f1 ' ' f2 ' ' f3 ' ' f4 ' ' f5 ' ' f6 ' >outfile']);
        delete(hwarn);
    case 'No, thank you!'
        disp('No new simulation.')
end

if ~exist('Summary', 'dir')
    mkdir('Summary');
end
% % READ IN PARAMETER FILE
fileinf = dir('*parms*');
fileID = fopen(fileinf.name);
C = textscan(fileID,'%f %*[^\n]');
C = C{:};
fclose(fileID);

fileinf = dir('Correlation*');
cfname = fileinf.name;

Nprotypes = C(1);
nint = C(2);
maxtime = C(3);
Nspecies = C(4);
Nrxn = C(5);
vol = C(6);
Vsphere = C(7);
checkpointiter = C(8);
stepwrite = C(9);
niter = C(10);
delt = C(11);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

fileinf = dir('concfile*');
fileID = fopen(fileinf.name);
C = textscan(fileID,'%s %f %s %*[^\n]');
fclose(fileID);

concs = C{1,2};
um2dm = 1e-5;
volinLitre = vol*um2dm^3;
NA = 6.022e23;

if Vsphere == 1
    R = (3/4/pi*vol)^(1/3);
    A = 4*pi*R^2;
else
    R = vol^(1/3);
    A = R^2;
end

for i = 1:length(C{1,1})
    
    if C{1,1}{i,1} == 'M'
        nummol(i) = round(concs(i)*A);
    else
        nummol(i) = round(concs(i)*volinLitre*NA);
    end
    
end

fileinf = dir('ProPart*');
fid = fopen(fileinf.name);

% make sure the file is not empty
finfo = dir(fileinf.name);
fsize = finfo.bytes;

tempd = ' ';
tempf = ' ';

for i = 1:Nprotypes
    
    tempd = [tempd ' ' '%d'];
    tempf = [tempf ' ' '%f'];
    
end

A = [];
B = [];

if fsize > 0
    
    % read the file
    it = 0;
    while ~feof(fid)
        
        it = it + 1;
        
        if mod(it,Nprotypes)==1
            sk = fscanf(fid, 'total complexes: %d iter: %d time: %f\n', [3 1]);
            %             total complexes: 25810 iter: 10000 time: 16.1967
            A = [A; sk(:)'];
        else
            sk = fscanf(fid, ['Pro: %d Nbound:' tempd ' AvgNumBound:' tempf ' \n'], [9 1]);
            %             Pro: 1 Nbound: 	0	0	0	0	 AvgNumBound: 	-nan	-nan	-nan	-nan
            B = [B; sk(:)'];
        end
        
    end
    
end

% close the file
fclose(fid);

dirs = diff(A(:,end));
list = find(dirs<0);

list = [list(:);length(A(:,1))];
Btemp = 0;

for i = 1:length(list)
    
    Btemp = Btemp + B(list(i)*(Nprotypes-1)-2:list(i)*(Nprotypes-1),:);
    
end

Btemp = Btemp/length(list);

for i = 2:Nprotypes
    
    Btemp(i-1,2:2+Nprotypes-1) = Btemp(i-1,2:2+Nprotypes-1)/nummol(i);
    
end

headertext = '#Protein ';
for i = 1:Nprotypes
    headertext = [headertext '%partnerswithpro' num2str(i-1) ' '];
end
for i = 1:Nprotypes
    headertext = [headertext 'AvgNumBoundwithpro' num2str(i-1) ' '];
end
fid = fopen('Summary/newPropartnersfileFINALtime.dat','wt');
fprintf(fid, '%s\n', headertext);

dlmwrite('Summary/newPropartnersfileFINALtime.dat',Btemp,'delimiter','\t','precision','%.4f', '-append');
fclose(fid);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 
t = A(1:list(1),end);
tmin = t(1);
tmax = t(end);

t=A(1:list(1),end);

if(t(1)>tmin)
    tmin = t(1);
end
if(t(end)<tmax)
    tmax = t(end);
end

for i = 1:length(list)-1
    
    t=A(list(i)+1:list(i+1),end);
    
    if(t(1)>tmin)
        tmin = t(1);
    end
    if(t(end)<tmax)
        tmax = t(end);
    end
    
end

gridtime = logspace(log10(tmin), log10(tmax), timefreq);

dumpmat = zeros(timefreq,1+(Nprotypes*(Nprotypes-1)));
dumpmat(:,1) = gridtime;

headertext = 'Time ';

for i = 1:Nprotypes-1
    for j = 1:Nprotypes
        
        headertext = [headertext 'pro' num2str(i) '&' 'pro' num2str(j-1) ' '];
        
        mytr{i,j} = B(((1:length(A(:,end)))-1)*(Nprotypes-1)+i,1+j);
        numcomp = mytr{i,j};
        
        t = A(1:list(1),end);
        if gridtime(1)<t(1)
            t(1)=gridtime(1);
        end
        if gridtime(end)>t(end)
            t(end)=gridtime(end);
        end
        
%         disp([i j (i-1)*Nprotypes+j])
        
        if(sum(numcomp))==0
            dumpmat(:,1+(i-1)*Nprotypes+j) = zeros(timefreq,1);
        else
            dumpmat(:,1+(i-1)*Nprotypes+j) = interp1(t,numcomp(1:list(1)),gridtime)';
        end
        
        for ii = 1:length(list)-1
            
            t = A(1:list(1),end);
            if gridtime(1)<t(1)
                t(1)=gridtime(1);
            end
            if gridtime(end)>t(end)
                t(end)=gridtime(end);
            end
            if(sum(numcomp))==0
                dumpmat(:,1+(i-1)*Nprotypes+j) = dumpmat(:,1+(i-1)*Nprotypes+j) + zeros(timefreq,1);
            else
                dumpmat(:,1+(i-1)*Nprotypes+j) = dumpmat(:,1+(i-1)*Nprotypes+j) + interp1(t,numcomp(list(ii)+1:list(ii+1)),gridtime)';
            end
            
        end
    end
end

dumpmat(:,2:end) = dumpmat(:,2:end)/niter;

fid = fopen('Summary/newProPartnerstime.dat','wt');
fprintf(fid, '%s\n', headertext);
% have a matrix M(N,3) ready to go
dlmwrite('Summary/newProPartnerstime.dat',dumpmat,'delimiter','\t','precision','%.4f', '-append');
fclose(fid);

% % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % %
complex = dlmread(cfname);
dirs = diff(complex(:,1));
list = find(dirs<0);

t = complex(1:list(1),2);
tmin = t(1);
tmax = t(end);

for i = 1:length(list)-1
    
    t=complex(list(i)+1:list(i+1),2);
    
    if(t(1)>tmin)
        tmin = t(1);
    end
    if(t(end)<tmax)
        tmax = t(end);
    end
    
end

gridtime = logspace(log10(tmin), log10(tmax), timefreq);
% % % % % % % % % % % % % % % % % % % % % % % %

dumpmat = zeros(timefreq,Nspecies+1);
dumpmat(:,1) = gridtime;

for j = 1:Nspecies
    
    t = complex(1:list(1),2);
    if gridtime(1)<t(1)
        t(1)=gridtime(1);
    end
    if gridtime(end)>t(end)
        t(end)=gridtime(end);
    end
    numcomp = complex(1:list(1),3+j-1);
    if(sum(numcomp))==0
        dumpmat(:,1+j) = zeros(timefreq,1);
    else
        dumpmat(:,1+j) = interp1(t,numcomp,gridtime)';
    end
end

for i = 1:length(list)-1
    
    for j = 1:Nspecies
        
        t = complex(list(i)+1:list(i+1),2);
        numcomp = complex(list(i)+1:list(i+1),3+j-1);
        if(sum(numcomp))==0
            temp = zeros(timefreq,1);
        else
            temp = interp1(t,numcomp,gridtime)';
        end
        
        dumpmat(:,1+j) = dumpmat(:,1+j)+temp;
        
    end
    
end

dumpmat(:,2:end) = dumpmat(:,2:end)/niter;
dlmwrite('Summary/newCorrelationfile.dat',dumpmat,'delimiter','\t','precision','%.4f');
% % % % % % % % % % % % % % % %
loglog(dumpmat(:,1),dumpmat(:,2:end));
xlabel('Time (s.)');
ylabel('# Molecules Avg.');
saveas(gcf, 'Summary/AvgMolecs.fig', 'fig')


% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
fileinf = dir('ComplexMeanSize*');
cmsfname = fileinf.name;

complex = dlmread(cmsfname);
dirs = diff(complex(:,1));
list = find(dirs<0);

timefreq = 1000;
t = complex(1:list(1),2);
tmin = t(1);
tmax = t(end);

for i = 1:length(list)-1
    
    t=complex(list(i)+1:list(i+1),2);
    
    if(t(1)>tmin)
        tmin = t(1);
    end
    if(t(end)<tmax)
        tmax = t(end);
    end
    
end

gridtime = logspace(log10(tmin), log10(tmax), timefreq);

% % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %

dumpmat = zeros(timefreq,3);
dumpmat(:,1) = gridtime;

for j = 1:2
    
    t = complex(1:list(1),2);
    numcomp = complex(1:list(1),3+j-1);
    if(sum(numcomp))==0
        dumpmat(:,1+j) = zeros(timefreq,1);
    else
        dumpmat(:,1+j) = interp1(t,numcomp,gridtime)';
    end
end


for i = 1:length(list)-1
    
    for j = 1:2
        
        t = complex(list(i)+1:list(i+1),2);
        numcomp = complex(list(i)+1:list(i+1),3+j-1);
        if(sum(numcomp))==0
            temp = zeros(timefreq,1);
        else
            temp = interp1(t,numcomp,gridtime)';
        end
        
        dumpmat(:,1+j) = dumpmat(:,1+j)+temp;
        
    end
    
end

dumpmat(:,2:end) = dumpmat(:,2:end)/niter;
headertext='Time (s) Mean#ProteinsinaComplex Tot#Complex';
fid = fopen('Summary/newComplexMeanSize.dat','wt');
fprintf(fid, '%s\n', headertext);
dlmwrite('Summary/newComplexMeanSize.dat',dumpmat,'delimiter','\t','precision','%.4f', '-append');
fclose(fid);
% % % % % % % % % % % % % % % %
figure
semilogx(dumpmat(:,1),dumpmat(:,2));
xlabel('Time (s.)');
ylabel('Mean Complex Size');
saveas(gcf, 'Summary/MeanCompSize.fig', 'fig')


figure
semilogx(dumpmat(:,1),dumpmat(:,3));
xlabel('Time (s.)');
ylabel('Total # Complexes');
saveas(gcf, 'Summary/TotalnumComp.fig', 'fig')

% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
fileinf = dir('ComplexHist*');
fileID = fopen(fileinf.name);
C = textscan(fileID,'%s %s %s %s %s %s %s %*[^\n]');
fclose(fileID);

locs = C{1,7};
list = [];
times = [];
dat = [];

for i=1:length(locs)
    
    if ~isempty(locs{i})
        list = [list i];
        times = [times str2num(locs{i})];
        if(length(list)-1)>0
            datcont{length(list)-1} = dat;
        end
        dat = [];
    else
        dat = [dat; str2num(C{1,1}{i,1}) ...
            str2num(C{1,2}{i,1}) ...
            str2num(C{1,3}{i,1}) ...
            str2num(C{1,4}{i,1})];
    end
end

datcont{length(list)} = dat;

dirs = diff(times);
list = find(dirs<0);

for i = 1:length(list)
    
    finalcont{i} = datcont{list(i)};
    
end

finalcont{length(list)+1} = datcont{end};

bigmat = zeros(3,4);

for i=1:length(finalcont)
    
    mat = finalcont{i};
    sizeofmat = length(mat(:,1));
    finalmat(1:2,:) = mat(1:2,:);
    
    finalmat(3,:) = zeros(1,length(mat(1,:)));
    
    for j = 3:sizeofmat
        finalmat(3,:) = finalmat(3,:) + mat(j,:);
    end
    
    bigmat = bigmat + finalmat/length(finalcont);
    
end

bigmat(3,1) = 3;

headertext='ComplexSize %ComplexThisSize TotComplexThisSizeComplex %TotalProteinsThisSizeComplex';
fid = fopen('Summary/newComplexHist.dat','wt');
fprintf(fid, '%s\n', headertext);
dlmwrite('Summary/newComplexHist.dat',bigmat,'delimiter','\t','precision','%.4f', '-append');
fclose(fid);

[~,~] = unix('mv outfile Summary/');
