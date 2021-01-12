% This is the main function for simulating austenite-ferrite phase transformations
% 
% Be able to coupled with magnetic and ND models
% First established in 2017
% updated on September 9th, 2020
% This model can be seen as a further development of a previous model using
% the effective interface mobity.
% Features: Voronoi structure, PBC, spherical ferrite grains, classical nucleation theory,
% mixed-mode,cyclic thermal cycle,
% Gibbs energy balance between chemical driving force and energy
% dissipations due to interface friction and solute drag
%
% Suitable for continuous cooling, isothermal heating as well as cycling in
% the gamma-alpha two phase region.
% Copyright by Haixing Fang

% At present only valid for Fe-xC-yMn ternary alloys
% Thermodynamic data has to be available and certain changes have to be made
% before this model can be applied to other alloys
% for running on PC

% Last update on Jan 4, 2021

clear all;
close all;

% load Multi-Parametric Toolbox 3
% https://www.mpt3.org/
% addpath('C:\Users\hfang\Documents\MATLAB\tbxmanager');
if exist('.\tbxmanager','dir')
    addpath('.\tbxmanager');
    tbxmanager restorepath;
    mpt_init;
else
    fprintf(['Error: the path of mpt3 toolbox is not found. \n', ...
        'Please install mpt3 toolbox first! \n', ...
        'Download from https://www.mpt3.org/ \n']);
end
% cite using MPT3:
% M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, pages 502?10, Zurich, Switzerland, July 17?9 2013.

% % % characteristics of steels with predefined composition
% Read the simulation conditions: sample information
[Comp T_a3 T_a1 Tc a_eq(5,:) b_eq(5,:) kafang_p T_wC Temp_eq(5,:) Feq_max RT hold_t cycN cycT CyclicFlag]=SimulCond();

% modelling parameters
T0=T_a3; % starting temperature [K]
rou_bcc=1e15; % threshold for ferrite nuclei density [m-3]
SoluteDragFlag=1; % 0-not consider of solute drag; 1-include solte drag effect
R=8.314; % gas constant [J/(mol.K)]
deltaS=3.5e5;  % Entropy [J/(K.m3)]
eps=1e-6; % minimum error

% starting microstrucuture
Lb=70; % box size [um]
Npot=fix(rou_bcc*(Lb*1e-6)^3); % maximum number of ferrite nuclei
f_N=1;       % scaling factor
N0=Npot*f_N; % scaled number of potential nucleation sites [-]
dfcc=20; % average austenite grain diameter [um]
dmin=0.25*dfcc;     % The minimum distance required between the centroids of austenite grains
N_fcc=fix(Lb^3/(pi*(dfcc^3)/6)); % number of austenite grains
rou_fcc=N_fcc/Lb^3; % Number density of austenite grain [um-3]
shieldD=1/4*(Lb^3/(N0/5))^(1/3); % 1/4 of the average ferrite spacing [um]
sigma_FA=0.62; % interfacial energy between frrite and austenite [J/m^2], H.I. Aaronson Metall.Mater.Trans.A 1988
%%%%% in the article it equals to 0.5 J/m^2 M. Segawa 2017 Comp. Mater

enlarge_Lb=Lb+2*(0.5*1/rou_fcc^(1/3));% enlarged by 2*average neibouring distance
enlarge_N=fix(rou_fcc*enlarge_Lb^3); % Number of austenite grains in enlargement box (enlarged by 2*average neibouring distance) that have the same number density
minus_edge=0.5*1/rou_fcc^(1/3);

% Define the anistropy axis for ferrite nucleation sites
MA_theta0=2*pi*rand(1,Npot);
MA_phi0=pi*rand(1,Npot);

% Assign the initial value for magnetization of each ferrite grain
alpha0=2*pi*rand(Npot,1);
beta0=pi*rand(Npot,1);

%General constantants
M_Fe=56;
M_C=12;
M_Mn=55;
M_Si=28;
Mole=Comp(1)/M_C+Comp(2)/M_Mn+Comp(3)/M_Si+(100-Comp(1)-Comp(2)-Comp(3))/M_Fe;%total mole number of elements
Comp_m(1)=Comp(1)/M_C/Mole;Comp_m(2)=Comp(2)/M_Mn/Mole;Comp_m(3)=Comp(3)/M_Si/Mole;%mole fraction of each element
Comp_m=100.*Comp_m; % [mol%]
Ux=Comp_m(2)/(Comp_m(2)+(100-Comp_m(1)-Comp_m(2)-Comp_m(3))); % U-fraction X/(Fe+X)

% Meff=3.5e-7; % effective interface mobility [m3.m/(J.s)] in Pina and C.Bos Meff=2.5~5e-7, E.Gamsjager think Meff=6.9~10e-7
M0=2.4e-7; % [m3.m/(J.s)]
Meff=M0/7.1e-6; % [mol.m/(J.s)]
QM=140e3; % [J/mol]

%Generate positions potential nucleation sites and plot the 3d voronoi diagram
generate_new_seed=0; % generate new seeding points: 1-yes; 0-no; by defaul to load the old seeding points
if generate_new_seed==1
    A1=enlarge_Lb*rand([1,3])-minus_edge;
    i=1;
    A=A1;
    while i<=enlarge_N-1
       A2=enlarge_Lb*rand([1,3])-minus_edge-Lb/2;
       comb=[A;A2];
       distance=pdist(comb,'euclidean');
       if min(distance)>dmin
          A=[A;A2];
          i=i+1
       end
    end
    % export the predefined austenite centroid coordinates
    dlmwrite(strcat(num2str(length(A)),'AusteniteCentroidsSymmetric.txt'),[A(:,1) A(:,2) A(:,3)],'delimiter',' ');
else
    % load the predefined austenite centroid coordinates
	if dfcc==20
        FileName_prefix='153AusteniteCentroidsSymmetric'; % randomly generate 153 symmetric centers
	else
	    FileName_prefix='151AusteniteCentroidsSymmetric_Lb175_d50'; % primarily for dfcc = 50 um
	end
    FileName_pattern='.txt';
    baseFileName=[FileName_prefix FileName_pattern];
    fullFileName=fullfile(pwd, baseFileName);
    fileID=fopen(fullFileName,'r');
    A=[];
    while(~feof(fileID))
        textdata=str2num(fgetl(fileID));
        if isempty(textdata)
            continue;
        else
            A=[A; textdata];
        end
    end
    fclose(fileID);
end
dis=pdist(A,'euclidean'); % Pairwise distance of austenite grain center
B=Polyhedron([-minus_edge-Lb/2 -minus_edge-Lb/2 -minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 -minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2; ...
    -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2; ...
    -minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2]);% Boundary vertices
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
vertices = V.Set.forEach(@(e) e.V, 'UniformOutput', false);% find the vertices for each cell
new_vertices = cat(1, vertices{:}); % Combine the vertices into one matrix
potential1=unique(round(new_vertices*1e6), 'rows')/1e6; % Get rid of the numerical noise

% common faces for each grain
for i=1:length(V.Set)
    CommonFaces{i}=[];
    for j=1:length(V.Set)
        if i~=j
            [ts iP iQ]=isAdjacent(V.Set(i),V.Set(j));
            if ts==1
                CommonFaces{i}=[CommonFaces{i};i j ts iP iQ]; % grain i j, isAdjacent, faceNr in grain i and j
            end
        end
    end
    i
end

% grain edge (on the line) and boundary (on the face)
edges=[];
boundaries=[];
for i=1:length(V.Set)
    faces{i}=V.Set(i).getFacet; % get all faces of each polyhedron
    for j=1:length(faces{i})
        edges_points=faces{i}(j).facetInteriorPoints;
        ShareGrNr=find(CommonFaces{i}(:,4)==j);
        for k=1:length(edges_points(:,1))
            if ~isempty(ShareGrNr)
                edges=[edges;edges_points(k,:) CommonFaces{i}(ShareGrNr,[1:2 4:5])]; % grain line [xyz, grain number in share, face nr in both grains]
            else
                edges=[edges;edges_points(k,:) i 0 j 0];
            end
        end
        bou=V.Set(i).chebyCenter(j);
        if ~isempty(ShareGrNr)
            boundaries=[boundaries;bou.x' CommonFaces{i}(ShareGrNr,[1:2 4:5])]; % grain face [xyz, grain number in share, face nr]
        else
            boundaries=[boundaries;bou.x' i 0 j 0];
        end
    end
    i
end

potential2=unique(edges(:,1:3),'rows');
potential3=unique(boundaries(:,1:3),'rows');
potential0=[potential1;potential2;potential3];
potential0(:,4)=[ones(length(potential1),1);2*ones(length(potential2),1);3*ones(length(potential3),1)];

to_remove = [-Lb/2 -Lb/2 -Lb/2;Lb/2 -Lb/2 -Lb/2;Lb/2 Lb/2 -Lb/2;-Lb/2 Lb/2 -Lb/2; ...
    -Lb/2 -Lb/2 Lb/2;Lb/2 -Lb/2 Lb/2;Lb/2 Lb/2 Lb/2;-Lb/2 Lb/2 Lb/2];% Remove the corners of the Lb*Lb*Lb cubic
[potential ia]=setdiff(potential0(:,1:3), to_remove, 'rows');
potential(:,4)=potential0(ia,4);

for i=1:length(potential(:,1))
    for j=1:3
    if potential(i,j)>Lb/2
       potential(i,j)=Lb/2+1;
    else if potential(i,j)<-Lb/2
        potential(i,j)=-Lb/2-1;
        end
    end
    end
end
site_all=[];
% Remove points outside of the Lb*Lb*Lb cubic box
for i=1:length(potential(:,1))
  if all(potential(i,:)-Lb/2-1) && all(potential(i,:)+Lb/2+1)
%      eval(['site',num2str(i),'=','potential(i,:)']); % Convert the number to string
%      eval(['site=[site;site',num2str(i),'];']);  % Combine site1, site 2,...
     site_all=[site_all;potential(i,:)];
  end
end
site_all=site_all(randperm(length(site_all(:,1))),:);

% define the fraction of types of potential nucleation sites:
% 1-corners, 2-edges, 3-boundaries
% fsite=[length(find(site_all(:,4)==1)) length(find(site_all(:,4)==2)) length(find(site_all(:,4)==3))]/length(site_all(:,1));
fsite=[1 0 0]; % manually defined, by default, only corners are potential nucleation sites
% fsite=[0.5 0.25 0.25]; % manually defined
site_potential=length(find(site_all(:,4)==1)); % total potential nucleation sites
fsite(2:3)=fix(fsite(2:3)*site_potential);
fsite(1)=site_potential-sum(fsite(2:3));
site=[];
for i=1:3
    site_index=find(site_all(:,4)==i);
    site=[site;site_all(site_index(1:fsite(i)),:)];
end

p1=randperm(length(site));
if N0>length(site)
    N0=length(site);
end
position=site(p1(1:N0),1:3);% Randomly select N0 potential nucleation site for ferrite
site=site(p1,:); % randomly sort the possible nucleation sites
site_type=site(:,4);
site=site(:,1:3);

figure('Name','Austenite geometry');
subplot(1,2,1);
V.plot('alpha', 0.13); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'b+','LineWidth',2);% Plot the centroids
axis([-enlarge_Lb/2 enlarge_Lb/2 -enlarge_Lb/2 enlarge_Lb/2 -enlarge_Lb/2 enlarge_Lb/2]);
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
% legend('+ centroids of austenite grains','. corners of austenite grains');
box on;
grid off;
hold off;

subplot(1,2,2);
V.plot('alpha', 0.13); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'b+','LineWidth',2);% Plot the centroids
% plot3(site_all(:,1),site_all(:,2),site_all(:,3),'r.','MarkerSize',20);
% plot3(site(:,1),site(:,2),site(:,3),'b*','LineWidth',2);% Plot all the potential nucleation sites
plot3(site(:,1),site(:,2),site(:,3),'r.','MarkerSize',20);
% siteA=site(find(site_type==1),:);
% siteB=site(find(site_type==2),:);
% siteC=site(find(site_type==3),:);
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
% axesLabelsAlign3D;
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
box on;
hold off;

% Define initial features for each austenite
for i=1:length(A)
   A_P(i,1)=i; % order number
   A_P(i,2)=V.Set(i).Data.voronoi.seed(1,1); % x coordinate
   A_P(i,3)=V.Set(i).Data.voronoi.seed(2,1); % y coordinate
   A_P(i,4)=V.Set(i).Data.voronoi.seed(3,1); % z coordinate
   A_P(i,5)=V.Set(i).volume; % volume [um3]
   A_P(i,6)=0; % number of nucleated ferrite in its vertex
   A_P(i,7)=length(V.Set(i).V); % total number of vertex
   A_P(i,8)=0;
   for j=1:length(V.Set(i).V)
       if abs(V.Set(i).V(j,1))>Lb/2 || abs(V.Set(i).V(j,2))>Lb/2 || abs(V.Set(i).V(j,3))>Lb/2
          A_P(i,8)=A_P(i,8)+1; % total number of vertex outside Lb
       end
   end
   A_P(i,9)=A_P(i,7)-A_P(i,8); % effective number of vertex(<Lb)
   A_P(i,10)=Comp(1); % C content in matrix [wt.%]
   A_P(i,11)=Comp(2); % Mn content in matrix [wt.%]
   A_P(i,12)=0; % ferrite volume in this specific austenite [um3]
   A_P(i,13)=A_P(i,12)/A_P(i,5); % ferrite fraction of each austenite
   A_P(i,14)=0; % soft impingement_1 or not_0
   A_P(i,15)=Comp_m(1); % C content in matrix [mol%]
   A_P(i,16)=Comp_m(2); % C content in matrix [mol%]
end
[position_MXY position_MXZ position_MYZ position_MO]=MirrorPoints(site,Lb); % Get the coordinates of the mirror points

% Generate a N0*34 matrix:siteID,x,y,z,active or not,nucleation time,radius,effective radius,nucleated flag
% actual volume,ratio of actual volume to the original volume, actual fraction/equilibrium fraction
% impinge 4 flag, average distance between vertex and center ...
for i=1:length(site)
    N_p(i,1)=i;%site number
    N_p(i,2)=site(i,1);%x coordinate
    N_p(i,3)=site(i,2);%y coordinate
    N_p(i,4)=site(i,3);%z coordinate
    N_p(i,5)=1; % active=1, not active=0
    N_p(i,6)=-1;% nucleation time, not nucleated=-1
    N_p(i,7)=0; % original radius
    N_p(i,8)=0; % effective radius
    N_p(i,9)=0; % not nucleated=0, nucleated=1
    N_p(i,10)=0;% actual volume (after substract the overlay volume)
    N_p(i,11)=-1;% actual volume/original volume,=-1 when original volume=0
    N_p(i,12)=0;% actual volume fraction/equilibrium volume fraction predicted by phase diagram
    N_p(i,13)=0;% The flag of impingement of 4 spheres (0 denotes no impingement while 1 denotes impingement)
    N_PR{i}=[];
    if site_type(i)==1
        for j=1:length(A)
              for k=1:length(V.Set(j).V)
               if abs(site(i,1)-V.Set(j).V(k,1))<=eps && abs(site(i,2)-V.Set(j).V(k,2))<=eps && abs(site(i,3)-V.Set(j).V(k,3))<=eps
                  N_PR{i}=[N_PR{i};j k]; % neighboring grain sharing the grain vertice and vertice number
               end
              end
        end
    elseif site_type(i)==2
        edges_index=find(edges(:,1)==site(i,1) & edges(:,2)==site(i,2) & edges(:,3)==site(i,3));
%         N_PR{i}=[N_PR{i};edges(edges_index,5:6)]; % neighboring grain sharing the grain line and face number
        if length(edges_index)==1
            N_PR{i}=[N_PR{i};edges(edges_index,[4 6]);edges(edges_index,[5 7])]; % can be shared up to 3 grains
        else
            edges_potential=[edges(edges_index,[4 6]);edges(edges_index,[5 7])];
            N_PR{i}=[N_PR{i};edges_potential([1 3:end],:)];
        end
    else
        boundaries_index=find(boundaries(:,1)==site(i,1) & boundaries(:,2)==site(i,2) & boundaries(:,3)==site(i,3));
        N_PR{i}=[N_PR{i};boundaries(boundaries_index,[4 6]);boundaries(boundaries_index,[5 7])]; % neighboring grain sharing the grain face and face number
    end           
    for m=1:length(N_PR{i}(:,1))
        N_PR{i}(m,3)=sqrt((site(i,1)-A_P(N_PR{i}(m,1),2))^2+(site(i,2)-A_P(N_PR{i}(m,1),3))^2+(site(i,3)-A_P(N_PR{i}(m,1),4))^2);% identical distance to neighbouring centers
    end
    N_p(i,14)=mean(N_PR{i}(:,3)); % average distance between vertex and austenite centers [um]
    N_p(i,15)=mean(A_P(N_PR{i}(:,1),10)); % average remote C content in austenite [wt.%]
    N_p(i,16)=mean(A_P(N_PR{i}(:,1),11)); % average remote Mn content in austenite [wt.%]
    N_p(i,17)=0.004; % average remote C content in ferrite at T_a3 [wt.%]
    N_p(i,18)=mean(A_P(N_PR{i}(:,1),11)); % average remote Mn content in ferrite [wt.%]
    N_p(i,19)=T0-T_a3; % local undercooling [K]
    N_p(i,20)=0; % G_fcc-G_bcc
    N_p(i,21)=1; % only becomes 1 means it can be potential nucleation sites
    N_p(i,22)=Comp(1); % C content at the site [wt.%]
    N_p(i,23)=Comp(2); % Mn content at the site [wt.%]
    N_p(i,24)=0; % interface velocity [um/s]
    N_p(i,25)=0; % flag of soft impingement either 0 or 1
    N_p(i,26)=0; % hard impingement either 0 or 1
    N_p(i,27)=N_p(i,14); % initial soft impingement distance [um]
    if i==1
       N_p(i,28)=min([mean(pdist([site(i,:);site(i+1:end,:)])) mean(pdist([position_MXY(i,:);site(i+1:end,:)])) ...
           mean(pdist([position_MXZ(i,:);site(i+1:end,:)])) mean(pdist([position_MYZ(i,:);site(i+1:end,:)])) ...
           mean(pdist([position_MO(i,:);site(i+1:end,:)]))]); % mean distance between vertices [um]
    else if i==length(site)
       N_p(i,28)=min([mean(pdist([site(i,:);site(1:i-1,:)])) mean(pdist([position_MXY(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MXZ(i,:);site(1:i-1,:)])) mean(pdist([position_MYZ(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MO(i,:);site(1:i-1,:)]))]); % mean distance between vertices [um]
        else
       N_p(i,28)=min([mean(pdist([site(i,:);site(1:i-1,:)])) mean(pdist([position_MXY(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MXZ(i,:);site(1:i-1,:)])) mean(pdist([position_MYZ(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MO(i,:);site(1:i-1,:)]))]); % mean distance between vertices [um]
       N_p(i,28)=mean([N_p(i,28) min(pdist([site(i,:);site(i+1:end,:)])) min(pdist([position_MXY(i,:);site(i+1:end,:)])) ...
           min(pdist([position_MXZ(i,:);site(i+1:end,:)])) min(pdist([position_MYZ(i,:);site(i+1:end,:)])) ...
           min(pdist([position_MO(i,:);site(i+1:end,:)]))]); % minimum distance between vertices [um]
        end
    end
    N_p(i,29)=0; 
    for m=1:length(N_PR{i}(:,1))
        N_p(i,29)=N_p(i,29)+A_P(N_PR{i}(m,1),5); % volume that it can maximum grow [um^3]
    end
    N_p(i,30)=-1; % restore the moment when extended volume is corrected [s]
    N_p(i,31)=0;  % restore the site number of pair growing ferrite
    N_p(i,32)=(3.*sum(A_P(N_PR{i}(:,1),5))/(4*pi)).^(1/3); % remaining untransformed surrounding austenite volume [um]
    N_p(i,33)=0;
    N_p(i,34)=Comp_m(1);
    N_pD(i,7)=0; % define an initial value for grain radius in diffusion-controlled mode [um]
    i;
end
NucPool=[];
for i=1:length(site)
    if i==1
        NucPool=[NucPool;N_p(i,1:4)];
        N_p(i,5)=1;
    else
        mindis=sqrt((N_p(i,2)-NucPool(:,2)).^2+(N_p(i,3)-NucPool(:,3)).^2+(N_p(i,4)-NucPool(:,4)).^2);
        mindis=min(mindis(mindis>0));
        if mindis<=shieldD
            N_p(i,5)=0;
        else
            NucPool=[NucPool;N_p(i,1:4)];
        end
    end
end

for j=1:length(A)
    A_PR{j}=[];
    for k=1:length(V.Set(j).V)
        for i=1:length(N_p(:,1))
           if abs(N_p(i,2)-V.Set(j).V(k,1))<=eps && abs(N_p(i,3)-V.Set(j).V(k,2))<=eps && abs(N_p(i,4)-V.Set(j).V(k,3))<=eps
%               A_PR{j}=[A_PR{j};N_p(i,1)]; % restore the site order number
              A_PR{j}(k)=N_p(i,1); % restore the site order number
           end
        end
    end
end

% parameters of thermal cycling
cyct=(cycT(2)-cycT(1))/RT(2)+hold_t(2)+(cycT(2)-cycT(1))/RT(3)+hold_t(3); % time of one cycle [s]
% critical time node [s]
tcr(1)=0; % start moment
if CyclicFlag>=0
    tcr(2)=(T0-cycT(1))/RT(1); % start of holding
    tcr(3)=hold_t(1)+tcr(2); %
else
    tcr(2)=(T0-Temp_eq(5,2))/RT(1);
    tcr(3)=tcr(2);
end
if CyclicFlag==1
    for i=1:cycN
       tcr(3+(i-1)*4+1)=(cycT(2)-cycT(1))/RT(2)+tcr(3+(i-1)*4);
       tcr(3+(i-1)*4+2)=hold_t(2)+tcr(3+(i-1)*4+1);
       tcr(3+(i-1)*4+3)=(cycT(2)-cycT(1))/RT(3)+tcr(3+(i-1)*4+2);
       tcr(3+(i-1)*4+4)=hold_t(3)+tcr(3+(i-1)*4+3);
    end
    tcr(3+(i-1)*4+5)=(cycT(1)-Temp_eq(5,2))/RT(4)+tcr(3+(i-1)*4+4); % time to reach lowest temperature [s]
end

%%%%%%%%%%%%%%%%%%%%%%
Timer(1)=0; % total crystallization time
i=0;j=1;k=1; % initial definition
start=0; % initiation of the first nuclei time
dd=1; % Count number for the calculation of magnetic_m
stop_cycle=0;
F(1)=0; % initial ferrite fraction
sumN(1)=0; % total ferrite nucleis
t=0;    % define initial t
if CyclicFlag==1
%     dt(1)=((cycT(2)-cycT(1))/RT(2))/50; % make sure there are 20-500 data points during the heating or cooling segment of cycling
    dt(1)=tcr(3)/1000;
    if dt(1)<0.5
        dt(1)=0.5;
    end
elseif CyclicFlag==0
    dt(1)=tcr(end)/1000; % for isothermal, 3000 iterations take about 24 h, keep iterations between 1000 and 3000
else
    dt(1)=tcr(end)/100; % for continuous cooling, about 3.5 h
end
dx=[1 0.005]; % max and min step size in length [um]
DiffInfo=cell([1 length(N_p(:,1))]);

while stop_cycle~=1 % main loop to calculate nucleation,growth and call impingement function
    t
    i=i+1;
    n=1;
    Nucleated{i}=[]; % restore all the nucleated ferrite
    N_pA{i}=[];
    Index_minD{i}=[];
    Timer(i)=t; % timer [s]
    
    %%%% temperature profile [K]
    for j=1:length(tcr)-1
        if t-tcr(j)>=0 && t-tcr(j+1)<=0
            t_index(1)=j;
            t_index(2)=j+1;
        end
    end
    if t_index(1)==1
        T(i)=T0-RT(1)*t;
    else if t_index(1)==2
            T(i)=cycT(1);
            else if t_index(1)>=3 && t_index(1)<length(tcr)-1
                    t_cycle(1)=fix((t_index(1)-3)/4)+1; % which cycle
                    t_cycle(2)=mod((t_index(1)-3),4); % which stage of a particular cycle
                    if t_cycle(2)==0
                        T(i)=cycT(1)+RT(2)*(t-tcr(3+4*(t_cycle(1)-1)));
                    else if t_cycle(2)==1
                            T(i)=cycT(2);
                        else if t_cycle(2)==2
                                T(i)=cycT(2)-RT(3)*(t-tcr(3+4*(t_cycle(1)-1)+2));
                            else if t_cycle(2)==3
                                    T(i)=cycT(1);
                                end
                            end
                        end
                    end
                else if t_index(1)==length(tcr)-1
                        T(i)=cycT(1)-RT(4)*(t-tcr(end-1));
                    end
                end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for m=1:length(N_p(:,1))
        N_p(m,19)=T(i)-(T_wC(1)*N_p(m,15)^2+T_wC(2)*N_p(m,15)+T_wC(3)); % updated on Sep 10, 2018
        N_p(m,20)=-deltaS*N_p(m,19);% driving force [J/m3] positive means enabling austenite to ferrite
        if N_p(m,20)<0
            N_p(m,21)=0; % only becomes 1 means it can be potential nucleation sites
        else
            N_p(m,21)=1;
        end
        if N_p(m,9)==0 && N_p(m,5)==1
           N_pA{i}(n,:)=N_p(m,:); % restore the potential available site
           n=n+1;
        end
    end
    if ~isempty(N_pA{i})
        deltaGV(i)=mean(N_pA{i}(:,20)); % driving force averaging the whole structure [J/m3]
    else
        deltaGV(i)=0;
    end
    
    %%%%% classical nucleation theory
    if Timer(i)==0
        Nt_cnt(i)=0;
        EnergyB(i)=0;
        Freq(i)=0;
        Scailing(i)=0;
    else
       [Nt_cnt(i) EnergyB(i) Freq(i) Scailing(i)]=CNT_cyclic_Gs(T(i),F(i-1),length(site), ...
	         deltaGV(i),Timer(i),dt(i),xC_F_eq(i-1),Comp_m); % [s-1]
    end
    if CyclicFlag==1 && Timer(i)>=tcr(3)
        %%%% no accumulation
        N_cnt(i)=Nt_cnt(i)*dt(i);
        if N_cnt(i)<1 || i==1
            N(i)=0;
        else
            N(i)=fix(N_cnt(i));
        end
        N_eff_cnt(i)=N_eff_cnt(i-1);
    else
        %%% accumulation
        if i==1
            Nreal(i)=Nt_cnt(i)*dt(i);
            Neff(i)=fix(Nt_cnt(i)*dt(i));
            Nrest(i)=Nreal(i)-Neff(i);
            N_eff_cnt(i)=Neff(i);
        else
            Nreal(i)=Nt_cnt(i)*dt(i)+Nrest(i-1);
            Neff(i)=fix(Nreal(i));
            Nrest(i)=Nreal(i)-Neff(i);
            N_eff_cnt(i)=N_eff_cnt(i-1)+Neff(i);
        end
    end
    if i==1
        N(i)=0;
    else
        N(i)=N_eff_cnt(i)-N_eff_cnt(i-1);
    end

    if i>=2
        sumN(i)=sumN(i-1)+N(i);
    %%%%% when using CNT or LinearNucleation, remember to comment this part
        if sumN(i)>=Npot
            sumN(i)=Npot;
            N(i)=sumN(i)-sumN(i-1);
        end
    end

    N_eff(i)=0;% N_eff is the real nucleation number
    if i>1 && sumN(i-1)==0 && sumN(i)>0
        start=Timer(i);% start moment of the 1st nuclei
    end
    numb=1;
    if ~isempty(N_pA{i})
    [sortN_p{i},sortIndex{i}] = sortrows(N_pA{i}, 22);% min wC should nucleat first
    if t>=start
        j=1;
       while (isempty(Nucleated{i})&&N(i)>0) || (~isempty(Nucleated{i})&&length(Nucleated{i}(:,1))<=N(i)-1)
           Index=N_pA{i}(sortIndex{i}(j),1); % site oder number
           % Update the distance from the nuclei
           if ~isempty(Nucleated{i-1}) 
              N_p(Index,28)=min([N_p(Index,28) min(pdist([N_p(Index,2:4);Nucleated{i-1}(1:end,2:4)])) ...
                   min(Index_minD{i})]);
               if ~isempty(Nucleated{i})
                   for k=1:length(Nucleated{i}(:,1))
                       kk=Nucleated{i}(k,1);
                       Index_minD{i}(k)=min([sqrt((N_p(Index,2)-position_MXY(kk,1))^2+(N_p(Index,3)-position_MXY(kk,2))^2+(N_p(Index,4)-position_MXY(kk,3))^2) ...
                       sqrt((N_p(Index,2)-position_MXZ(kk,1))^2+(N_p(Index,3)-position_MXZ(kk,2))^2+(N_p(Index,4)-position_MXZ(kk,3))^2) ...
                       sqrt((N_p(Index,2)-position_MYZ(kk,1))^2+(N_p(Index,3)-position_MYZ(kk,2))^2+(N_p(Index,4)-position_MYZ(kk,3))^2) ...
                       sqrt((N_p(Index,2)-position_MO(kk,1))^2+(N_p(Index,3)-position_MO(kk,2))^2+(N_p(Index,4)-position_MO(kk,3))^2)]);
                   end
                   N_p(Index,28)=min([N_p(Index,28) min(Index_minD{i})]);
               end
           end
           if N_p(Index,5)==1 && N_p(Index,21)==1 && N_p(Index,9)==0 && N_p(Index,28)>shieldD% active state and driving force > 0 and not nucleated
             N_p(Index,6)=Timer(i);%the specific nucleation time for particle j
             N_p(Index,7)=2*sigma_FA/N_p(Index,20)*1e6; % nucleus radius [um]
             if N_p(Index,7)>0.05
                 N_p(Index,7)=0.05;
             end
             N_pD(Index,7)=5e-3; % nucleus radius [um] for diffusion-controlled mode
             N_p(Index,9)=1;% particle j nucleated
             N_p(Index,34)=100*(N_p(Index,15)/M_C)/(N_p(Index,15)/M_C+N_p(Index,16)/M_Mn+(100-N_p(Index,15)-N_p(Index,16))/M_Fe); % the remote matrix C content when it starts to nucleate [at.%]
             Nucleated{i}=[Nucleated{i};N_p(Index,:)]; % restore the nucleated ferrite
             N_eff(i)=N_eff(i)+1;%N_eff is the real nucleation number
             j=j+1;
           else
              j=j+1;
           end
           if j>length(sortIndex{i})
               break;
           end
       end
       if i>1
          N_eff(i)=N_eff(i)+N_eff(i-1);
          Nucleated{i}=[Nucleated{i-1};Nucleated{i}]; % combine all the nucleated ferrite
       end
    end
    else
        Nucleated{i}=[Nucleated{i-1};Nucleated{i}];
        N_eff(i)=N_eff(i-1);
    end
    
    if ~isempty(Nucleated{i})
       numb=length(Nucleated{i}(:,1));% let numb=j,reveals the number of nucleated particle
       if numb>N0
          numb=N0;
       end
    else
        numb=1;
    end
      if T(i)<T_a3 && T(i)>=T_a1(2)  %When T is above A3-temperature, there is no C redistribution
         wC_A_eq(i)=a_eq(5,1)*T(i)^2+a_eq(5,2)*T(i)+a_eq(5,3);% This equation of A3-line is fitted from TC data
         wC_F_eq(i)=b_eq(5,1)*T(i)+b_eq(5,2);% This equation of C solubility in ferrite is fitted from TC data
      else if T(i)<T_a1(2)
              wC_A_eq(i)=a_eq(5,1)*T(i)^2+a_eq(5,2)*T(i)+a_eq(5,3);% This equation of A3-line is fitted from TC data
              wC_F_eq(i)=b_eq(5,1)*T(i)+b_eq(5,2);% T
          else
             wC_A_eq(i)=Comp(1);
             wC_F_eq(i)=0;
          end
      end
    xC_A_eq(i)=100*(wC_A_eq(i)/M_C)/(wC_A_eq(i)/M_C+(100-wC_A_eq(i))*Ux/M_Mn+(100-wC_A_eq(i))*(1-Ux)/M_Fe); % [mol%]
    xC_F_eq(i)=100*(wC_F_eq(i)/M_C)/(wC_F_eq(i)/M_C+(100-wC_F_eq(i))*Ux/M_Mn+(100-wC_F_eq(i))*(1-Ux)/M_Fe); % [mol%]
          
    F_eq(i)=(wC_A_eq(i)-Comp(1))/(wC_A_eq(i)-wC_F_eq(i));%Equilibrium ferrite fraction predicated by the phase diagram at time t
    if F_eq(i)>Feq_max
        F_eq(i)=Feq_max;
    end
    if F_eq(i)<0
        F_eq(i)=1e-12;
    end

    % below is to calculate the ferrite growth
    if i==1
       x_C=Comp_m(1)/100; % x_C in at.
    else
       x_C=(Comp_m(1)-F(i-1)*xC_F_eq(i-1))/(1-F(i-1))/100; % x_C in at.
    end
    y_C=x_C/(1-x_C);
    
    D_C(i)=4.53e-7*((1+y_C*(1-y_C)*8339.9/T(i))*exp(-(1/T(i)-2.221e-4)*(17767-26436*y_C)))*1e12;% Volume diffusion of Carbon in um2/s, J Agren 1986
    D_C1(i)=2.343e-5*exp(-148e3/(R*T(i)))*1e12; % Volume diffusion of Carbon in austenite um2/s, J Agren 1986
    D_C2(i)=1.5e-5*exp(-142.1e3/(R*T(i)))*1e12; % Volume diffusion of Carbon in austenite um2/s, R.C.Weast, 1989; C.Bos and J. Sietsma, 2007
    Mobility(i)=Meff*exp(-QM/(R*T(i)));% effective interface mobility [mol.m/(J.s)]
    Kafang(i)=kafang_p(1)*T(i)^2+kafang_p(2)*T(i)+kafang_p(3); % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
%     plot(T,D_C,'ro-',T,D_C1,'b.-',T,D_C2,'m*-')
%     legend('site fraction-J Agren 1986','Arrhenius-J Agren 1986','Arrhenius-R C Weast 1989')
    if i>1
        if T(i)-T(i-1)<=0
            TransDir=1; % transformation direction, positive: ferrite formation
        else
            TransDir=-1;
        end
    else
        TransDir=1;
    end
    ll=1;
    if i>1 && ~isempty(Nucleated{i})
    while ll<=length(Nucleated{i}(:,1))
      l=Nucleated{i}(ll,1); % site oder number
         if t>=N_p(l,6) && N_p(l,6)>0  % && N_p(l,13)==0
             if SoluteDragFlag==0 % no solute drag effect
                 for r=1:length(N_PR{l}(:,1))
                     if N_p(l,26)==0  % no hard impingement
                         [DiffInfo{l}(i,r,:)]=Mixmode_diffusion_profile_vpasolve(xC_F_eq(i),xC_A_eq(i),N_p(l,34),D_C(i)*1e-12,Mobility(i),Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32),N_p(l,26));% mix-mode update in April 2018
                         HardFlag{l}(i,r)=0; % hard flag
                         if DiffInfo{l}(i,r,2)<0 %|| DiffInfo{l}(i,r,2)>xC_A_eq(i) % fail to identify the hard impingement
                             DiffInfo{l}(i,r,1)=(Comp_m(1)-N_p(l,12)*xC_F_eq(i))/(1-N_p(l,12));
                             DiffInfo{l}(i,r,2)=DiffInfo{l}(i,r,1);
                             DiffInfo{l}(i,r,3)=0;
                             DiffInfo{l}(i,r,7)=1;
                             N_p(l,25)=1;
                             N_p(l,26)=1;
                         end
                         if DiffInfo{l}(i,r,3)<0
                             DiffInfo{l}(i,r,3)=N_p(l,27)-DiffInfo{l}(i,r,4);
                         end
                         % DiffInfo=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag];
                         % mol%, mol%, um, um, mol%, mol%, 1or0;
                     else
                         DiffInfo{l}(i,r,1)=(Comp_m(1)-N_p(l,12)*F_eq(i)*xC_F_eq(i))/(1-N_p(l,12)*F_eq(i));
                         DiffInfo{l}(i,r,2)=DiffInfo{l}(i,r,1); % Xip [mol%] 
                         [DiffInfo{l}(i,r,:)]=Mixmode_diffusion_profile_vpasolve(xC_F_eq(i),xC_A_eq(i),N_p(l,34),D_C(i)*1e-12,Mobility(i),Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32),N_p(l,26));% mix-mode update in April 2018
                         DiffInfo{l}(i,r,7)=1; % soft flag
                         HardFlag{l}(i,r)=1; % hard flag
                     end
                 end
                % mixed-mode model
                v_t(i,l)=Mobility(i)*Kafang(i)*(xC_A_eq(i)-mean(DiffInfo{l}(i,:,2)))*1e6; % velocity [um/s] mix-mode
                N_p(l,24)=v_t(i,l);
                N_p(l,7)=N_p(l,7)+v_t(i-1,l)*dt(i); %% grain radius [um]
             else % consider solute drag effect
                     % this is only for debugging
%                      r=1;
%                      i=74;
%                      l=1;
%                      ll=find(Nucleated{i}(:,1)==l);
%                      N_p(l,34)=Nucleated{i}(ll,34);
%                      N_p(l,7)=Nucleated{i-1}(ll,7);
%                      N_p(l,8)=Nucleated{i-1}(ll,8);
%                      N_p(l,32)=Nucleated{i-1}(ll,32);
%                      N_p(l,12)=Nucleated{i-1}(ll,12);
%                      N_p(l,26)=Nucleated{i-1}(ll,26);
%                      N_p(l,15)=Nucleated{i-1}(ll,15);
%                      N_p(l,17)=Nucleated{i-1}(ll,17);
                          DragIter=1;
                          StopDrag=0;
                          if ~exist('v_t') || length(v_t(i-1,:))<l
                                if Timer(i)>tcr(3) && CyclicFlag==1
                                    V0=[5e-8 -5e-8]; % [m/s]
                                else
                                    V0=[5e-8 1e-8]; % [m/s]
                                end
                          else
                                if Timer(i)>tcr(3) && CyclicFlag==1
                                    V0=[v_t(i-1,l)/1e6 -5e-8]; % [m/s]
                                else
                                    V0=[v_t(i-1,l)/1e6 5e-8]; % [m/s]
                                end
                          end
                          if V0(1)==0 || abs(V0(1))>1e-6
                              V0(1)=5e-8; % [m/s]
                          end
                          DragTemp=[];
                          DiffTemp=[];
                          while StopDrag~=1
                              if isempty(DiffInfo{l})
                                 [DragInfo{l}(i,1,:) DiffInfo{l}(i,1,:) GibbsE]=loopV_solute_drag(V0(DragIter),T(i),Comp_m(2),N_p(l,15),N_p(l,17),Ux, ...
                                   xC_F_eq(i),xC_A_eq(i),N_p(l,34),D_C(i)*1e-12,Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32),N_p(l,26), ...
                                   F(i-1),CyclicFlag,tcr(3),Timer(i),TransDir);
                              else
                                 [DragInfo{l}(i,1,:) DiffInfo{l}(i,1,:) Gdiff0{l}(i)]=solute_drag_fsolver(V0(DragIter),T(i),Comp_m(2),N_p(l,15),N_p(l,17),Ux, ...
                                   xC_F_eq(i),xC_A_eq(i),N_p(l,34),D_C(i)*1e-12,Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32),N_p(l,26),F(i-1),DiffInfo{l}(i-1,1,2), ...
                                   CyclicFlag,tcr(3),Timer(i),TransDir);
                              end
                              if DragInfo{l}(i,1,6)==1
                                  StopDrag=1;
                              end
                              if DragInfo{l}(i,1,6)>0 && DragInfo{l}(i,1,6)~=1
                                  DragTemp=[DragTemp;DragInfo{l}(i,1,1) DragInfo{l}(i,1,2) DragInfo{l}(i,1,3) DragInfo{l}(i,1,4) ...
                                      DragInfo{l}(i,1,5) DragInfo{l}(i,1,6)];
                                  DiffTemp=[DiffTemp;DiffInfo{l}(i,1,1) DiffInfo{l}(i,1,2) DiffInfo{l}(i,1,3) DiffInfo{l}(i,1,4) ...
                                      DiffInfo{l}(i,1,5) DiffInfo{l}(i,1,6) DiffInfo{l}(i,1,7)];
                              end
                              if DragIter==2
                                  StopDrag=1;
                                  if ~isempty(DragTemp)
                                      [SolDrag,SolIndex]=min(abs(DragTemp(:,2)-DragInfo{l}(i-1,1,2)));
                                      DragInfo{l}(i,1,:)=DragTemp(SolIndex,:);
                                      DiffInfo{l}(i,1,:)=DiffTemp(SolIndex,:);
                                  end
                              end
                              DragIter=DragIter+1;
                          end

                  if DragInfo{l}(i,1,6)<=0 || (DragInfo{l}(i,1,1)<0 && CyclicFlag<=0) ...
                          || ((abs(DragInfo{l}(i,1,1))>2e-6 || (Timer(i)<=tcr(3) && DragInfo{l}(i,1,1)<0)) && DragInfo{l}(i,1,6)>0)
                   %%%%%% this function use vpasolve: slower
                   %%%%%% tests show the following output the same result,
                   %%%%%% but, it can stay for further validation
                   if abs(DragInfo{l}(i,1,1))>=0.8e-6 || (DragInfo{l}(i,1,1)<=0 && CyclicFlag<=0) || (DragInfo{l}(i,1,6)>0 && CyclicFlag==1) || ...
                           (DragInfo{l}(i,1,6)<=0 && DragInfo{l}(i,1,3)~=0 && CyclicFlag==1)
                       [DragInfo{l}(i,1,:) DiffInfo{l}(i,1,:)]=solute_drag_dissipation_InterX_drag_analytical(T(i),Comp_m(2),N_p(l,15),N_p(l,17),Ux, ...
                        xC_F_eq(i),xC_A_eq(i),N_p(l,34),D_C(i)*1e-12,Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32),N_p(l,26),F(i-1),tcr(3),Timer(i),TransDir);
                       
                        [DragInfo{l}(i,1,:) DiffInfo{l}(i,1,:) GibbsE]=loopV_solute_drag(V0(1),T(i),Comp_m(2),N_p(l,15),N_p(l,17),Ux, ...
                            xC_F_eq(i),xC_A_eq(i),N_p(l,34),D_C(i)*1e-12,Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32),N_p(l,26), ...
                            F(i-1),CyclicFlag,tcr(3),Timer(i),TransDir); % option to use conventional loopV method
                   end
                   if DragInfo{l}(i,1,6)<=0
                        DragInfo{l}(i,1,1)=0;
                   end
                  end
                  for r=2:length(N_PR{l}(:,1))
                      DragInfo{l}(i,r,:)=DragInfo{l}(i,1,:);
                      DiffInfo{l}(i,r,:)=DiffInfo{l}(i,1,:);
                  end
                  % DragInfo=[Velocity Gchem Gfriction Gdiff Gdiss exitflag];
                  % Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag]';
                  Velocity=mean(DragInfo{l}(i,:,1)); % [m/s]
                  Gchem(i,l)=mean(DragInfo{l}(i,:,2)); % chemical driving force [J/mol]
                  Gfriction(i,l)=mean(DragInfo{l}(i,:,3)); % friction dissipation [J/mol]
                  Gdiff(i,l)=mean(DragInfo{l}(i,:,4)); % trans-diffusion dissipation [J/mol]
                  % solute drag model
                  v_t(i,l)=Velocity*1e6; % velocity [um/s]
                  N_p(l,24)=v_t(i,l); % [um/s]
                  N_p(l,7)=N_p(l,7)+v_t(i-1,l)*dt(i); %% grain radius [um]
                  if N_p(l,7)>N_p(l,32)*Feq_max
                      N_p(l,7)=N_p(l,32)*Feq_max;
                  end
                  if N_p(l,7)<0
                      N_p(l,7)=0;
                  end
             end
         end
             if ll==length(Nucleated{i}(:,1))
                 for aa=1:length(A)
                      SoftPair{aa}=[];
                      HardPair{aa}=[];
                      matrixC=[];
                     for bb=1:length(A_PR{aa})
                         if A_PR{aa}(bb)~=0 && N_p(A_PR{aa}(bb),9)==1
                             bb_dis=0;
                             dis_ferrite=[];
                             GrowDis=[];
                             DiffDis=[];
                         for cc=1:length(A_PR{aa})
                             if cc~=bb && A_PR{aa}(cc)~=0 && N_p(A_PR{aa}(cc),9)==1
                                bb_dis=bb_dis+1;
                                dis_ferrite(bb_dis)=sqrt((N_p(A_PR{aa}(cc),2)-N_p(A_PR{aa}(bb),2))^2+(N_p(A_PR{aa}(cc),3)-N_p(A_PR{aa}(bb),3))^2+ ...
                                    (N_p(A_PR{aa}(cc),4)-N_p(A_PR{aa}(bb),4))^2); % distance [um]
                                if (bb_dis>1 && dis_ferrite(bb_dis)<dis_ferrite(bb_dis-1)) || bb_dis==1
                                if (mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(cc)}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))<dis_ferrite(bb_dis) && N_p(A_PR{aa}(bb),25)==0
                                    N_p(A_PR{aa}(bb),25)=0; % no soft impingement
                                else if ((mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))<dis_ferrite(bb_dis) && (mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(cc)}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))>=dis_ferrite(bb_dis)) || N_p(A_PR{aa}(bb),25)==1
                                    GrowDis(aa,A_PR{aa}(bb),A_PR{aa}(cc))=mean(DiffInfo{A_PR{aa}(bb)}(i,:,3))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4));
                                    DiffDis(aa,A_PR{aa}(bb),A_PR{aa}(cc))=dis_ferrite(bb_dis)-mean(DiffInfo{A_PR{aa}(cc)}(i,:,3))-mean(DiffInfo{A_PR{aa}(cc)}(i,:,4));
                                    N_p(A_PR{aa}(bb),25)=1; % soft impingement
                                    SoftPair{aa}=[SoftPair{aa};A_PR{aa}(bb) A_PR{aa}(cc) dis_ferrite(bb_dis)]; % pair of soft impingement
                                    if ((mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))>=dis_ferrite(bb_dis)) || (N_p(A_PR{aa}(bb),25)==1 && N_p(A_PR{aa}(bb),26)==1)
                                        N_p(A_PR{aa}(bb),25)=1; % soft flag
                                        N_p(A_PR{aa}(bb),26)=1; % hard impingement flag
                                        HardPair{aa}=[HardPair{aa};A_PR{aa}(bb) A_PR{aa}(cc) dis_ferrite(bb_dis)]; % pair of hard impingement
                                    end
                                    end
                                end
                                end
                             end
                         end
                         
                         if isempty(dis_ferrite) && N_p(A_PR{aa}(bb),31)~=0
                             cc=N_p(A_PR{aa}(bb),31);
                             bb_dis=bb_dis+1;
                                dis_ferrite(bb_dis)=min([d(A_PR{aa}(bb),cc) d_MO(A_PR{aa}(bb),cc) d_MXY(A_PR{aa}(bb),cc) d_MXZ(A_PR{aa}(bb),cc) d_MYZ(A_PR{aa}(bb),cc)]);
                                if (bb_dis>1 && dis_ferrite(bb_dis)<dis_ferrite(bb_dis-1)) || bb_dis==1
                                if (mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{cc}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))<dis_ferrite(bb_dis) && N_p(A_PR{aa}(bb),25)==0
                                    N_p(A_PR{aa}(bb),25)=0; % no soft impingement
                                else if ((mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))<dis_ferrite(bb_dis) && (mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{cc}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))>=dis_ferrite(bb_dis)) || N_p(A_PR{aa}(bb),25)==1
                                    GrowDis(aa,A_PR{aa}(bb),cc)=mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3));
                                    DiffDis(aa,A_PR{aa}(bb),cc)=mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3));
                                    N_p(A_PR{aa}(bb),25)=1; % soft impingement
                                    SoftPair{aa}=[SoftPair{aa};A_PR{aa}(bb) cc dis_ferrite(bb_dis)]; % pair of soft impingement
                                    if ((mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))>=dis_ferrite(bb_dis)) || (N_p(A_PR{aa}(bb),25)==1 && N_p(A_PR{aa}(bb),26)==1)
                                        N_p(A_PR{aa}(bb),25)=1; % soft flag
                                        N_p(A_PR{aa}(bb),26)=1; % hard impingement flag
                                        HardPair{aa}=[HardPair{aa};A_PR{aa}(bb) cc dis_ferrite(bb_dis)]; % pair of hard impingement
                                    end
                                    end
                                end
                                end
                         end
                           matrixC(aa,A_PR{aa}(bb))=mean(DiffInfo{A_PR{aa}(bb)}(i,:,1)); % matrix C [mol%]
                             if ~isempty(GrowDis) && ~isempty(DiffDis)     
                                N_p(A_PR{aa}(bb),27)=min([N_p(A_PR{aa}(bb),27) min(nonzeros(GrowDis(aa,A_PR{aa}(bb),:)))]);
                                DiffInfo{A_PR{aa}(bb)}(i,:,3)=N_p(A_PR{aa}(bb),27)-mean(DiffInfo{A_PR{aa}(bb)}(i,:,4));
                                if DiffInfo{A_PR{aa}(bb)}(i,:,3)<0
                                    DiffInfo{A_PR{aa}(bb)}(i,:,3)=0;
                                end
                             end
                         end
                     end
                        if ~isempty(matrixC)
                         A_P(aa,15)=mean(nonzeros(matrixC(aa,:))); % matrix C [mol%]
                         A_P(aa,16)=100*(100-A_P(aa,15))*Ux/(A_P(aa,15)+(100-A_P(aa,15))*Ux+(100-A_P(aa,15))*(1-Ux)); % Mn concentent [mol%]
                         A_P(aa,10)=100*A_P(aa,15)*M_C/(A_P(aa,15)*M_C+(100-A_P(aa,15))*Ux*M_Mn+(100-A_P(aa,15))*(1-Ux)*M_Fe); % maxtrix C [wt.%]
                         A_P(aa,11)=100*(100-A_P(aa,15))*Ux*M_Mn/(A_P(aa,15)*M_C+(100-A_P(aa,15))*Ux*M_Mn+(100-A_P(aa,15))*(1-Ux)*M_Fe); % maxtrix Mn [wt.%]
                        end
                 end
             end                    
        ll=ll+1;
    end
    Ndt=[300 30]; % number of calculation points for one heating/cooling segment
    % too small could lead to oscillation results. Recommend values are [300 30]
         if max(abs(v_t(i,:)))~=0
            dt1=dx(1)/max(abs(v_t(i,:))); % update time step [s]
            dt2=dx(2)/max(abs(v_t(i,:))); % update time step [s]
            if dt1>=dt(1) && dt2<=dt(1)
                dt(i+1)=dt(1);
            elseif dt1<dt(1)
                    dt(i+1)=dt1;
            elseif dt2>dt(1)
                    dt(i+1)=dt(1);
            end
         elseif all(v_t(i,:)==0)==1 && Timer(i)>tcr(2)+(tcr(3)-tcr(2))*1/4 && Timer(i)<tcr(3)-1-dt(i)
             dt(i+1)=tcr(3)-1-Timer(i)-dt(i);
            if Timer(i)+dt(i)+dt(i+1)>tcr(3)
                 if RT(2)<5
                    dt(i+1)=min([dt(i+1) (cycT(2)-cycT(1))./RT(2)/Ndt(1)]);
                 else
                    dt(i+1)=min([dt(i+1) (cycT(2)-cycT(1))./RT(2)/Ndt(2)]);
                 end
            end
         else
             dt(i+1)=dt(1);
             if Timer(i)+dt(i)+dt(i+1)>tcr(3)
                 if RT(2)<5
                    dt(i+1)=min([dt(i+1) (cycT(2)-cycT(1))./RT(2)/Ndt(1)]);
                 else
                    dt(i+1)=min([dt(i+1) (cycT(2)-cycT(1))./RT(2)/Ndt(2)]);
                 end
             end
         end
         if CyclicFlag==1 && Timer(i)>=tcr(3)-dt(1) && Timer(i)<=tcr(end-1)+dt(1)
             if RT(2)<5
                dt(i+1)=min([dt(i+1) (cycT(2)-cycT(1))./RT(2)/Ndt(1)]);
             else
                dt(i+1)=min([dt(i+1) (cycT(2)-cycT(1))./RT(2)/Ndt(2)]);
             end
         end
    else
        dt(i+1)=dt(1);
    end
    if exist('v_t','var')
        if all(abs(v_t(i,:))<=5e-4) && Timer(i)>tcr(2) && Timer(i)<tcr(3)-(tcr(3)-tcr(2))/15 && CyclicFlag==1
            dt(i+1)=(tcr(3)-tcr(2))/15;
            if Timer(i)+dt(i)+dt(i+1)>tcr(3)
                dt(i+1)=dt(1);
            end
        end
    end

    mm=1;p4=0;p3=0;delta(i)=0;volume=0;
    if ~isempty(Nucleated{i})
    % Restore the coordinates and radious of particle and its mirror
    % partice: an array of [5*numb,4]
    for j=1:length(Nucleated{i}(:,1))
        xx(j,1)=N_p(Nucleated{i}(j,1),2);
        xx(j,2)=N_p(Nucleated{i}(j,1),3);
        xx(j,3)=N_p(Nucleated{i}(j,1),4);
        xx(j,4)=N_p(Nucleated{i}(j,1),7);
        xx(length(Nucleated{i}(:,1))+j,1)=position_MXY(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))+j,2)=position_MXY(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))+j,3)=position_MXY(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))+j,4)=N_p(Nucleated{i}(j,1),7);
        
        xx(length(Nucleated{i}(:,1))*2+j,1)=position_MXZ(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))*2+j,2)=position_MXZ(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))*2+j,3)=position_MXZ(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))*2+j,4)=N_p(Nucleated{i}(j,1),7);
        
        xx(length(Nucleated{i}(:,1))*3+j,1)=position_MYZ(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))*3+j,2)=position_MYZ(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))*3+j,3)=position_MYZ(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))*3+j,4)=N_p(Nucleated{i}(j,1),7);
        
        xx(length(Nucleated{i}(:,1))*4+j,1)=position_MO(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))*4+j,2)=position_MO(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))*4+j,3)=position_MO(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))*4+j,4)=N_p(Nucleated{i}(j,1),7);
    end   
    %%%%%%%% do not calculate the overlapped volume when the particle
    %%%%%%%% impinge with 4 or more other particles
    [Vol Tri impinge4_flag]=impingement34_PBC_modify(xx,length(Nucleated{i}(:,1))*5, ...
        [Nucleated{i}(:,13),Nucleated{i}(:,13),Nucleated{i}(:,13),Nucleated{i}(:,13),Nucleated{i}(:,13)]);% call the impingement34_PBC function
    
    % Calculate the substracted volume for particle i
    for n=1:numb       
        Vol_PBC(n)=Vol(n);
        if impinge4_flag(n)==1||impinge4_flag(n+numb)==1||impinge4_flag(n+numb*2)==1||impinge4_flag(n+numb*3)==1||impinge4_flag(n+numb*4)==1
           impinge4_flag_PBC(n)=1;
        else
           impinge4_flag_PBC(n)=0;
        end
    end  
    while mm<=length(Nucleated{i}(:,1)) %calculate the mean value and the volume fraction
        m=Nucleated{i}(mm,1);
        if N_p(m,7)>0
        N_p(m,8)=(3*Vol_PBC(mm)/(4*pi))^(1/3);%calculate the effective radius
        N_p(m,13)=impinge4_flag_PBC(mm);% Transfer the impinge3_flag to the matrix N_p
        if N_p(m,13)==1 && Timer(i)>N_p(m,6)
           if CyclicFlag==1 && Timer(i)>=tcr(3)
               StableIndex=max(find(Timer-tcr(3)>-4*mean(dt) & Timer<=tcr(3)));
               if isempty(StableIndex)
                   StableIndex=max(find(Timer-tcr(3)>-10*mean(dt) & Timer<=tcr(3)));
               end
               N_p(m,8)=Nucleated{i}(mm,7)*Nucleated{StableIndex}(mm,8)/Nucleated{StableIndex}(mm,7);
               if N_p(m,8)<0
                   N_p(m,8)=0;
               end
           else
               N_p(m,8)=Nucleated{i-1}(mm,8)+0.5*tanh(Nucleated{i-1}(mm,8)/N_p(m,7))*v_t(i,m)*dt(i);% correction on Dec 19
           if N_p(m,8)<Nucleated{i}(mm,8) || N_p(m,8)>1.05*Nucleated{i}(mm,8)
               if ismember(m,Nucleated{i-2}(:,1)) && ismember(m,Nucleated{i-1}(:,1))
                   N_p(m,8)=Nucleated{i-1}(mm,8)+1/2*((Nucleated{i-1}(mm,8)-Nucleated{i-2}(mm,8))/dt(i))*dt(i);
               else
                   N_p(m,8)=N_p(m,8)+v_t(i,m)*dt(i);
               end
           end
           end
           N_p(m,30)=t; % restore the moment when it includes correction of extended volume
        end
        if N_p(m,8)<(N_p(m,7)-eps) && CyclicFlag<=0
            N_p(m,25)=1;
            N_p(m,26)=1;
        end
        
        N_p(m,10)=4/3*pi*N_p(m,8)^3;%actual volume calculate from the effective radius
        N_p(m,11)=N_p(m,10)/(4/3*pi*N_p(m,7)^3);%actual volume/original volume,should always be <=1
        p3=p3+N_p(m,8)^3;
        p4=p4+N_p(m,8)^4;
        volume=volume+4/3*pi*N_p(m,8)^3; %used for the impingement34 function
        R0(i,m)=N_p(m,7);
        Reff(i,m)=N_p(m,8);
        GrDis(i,m)=N_p(m,32);
        end
        mm=mm+1;
    end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update details of each austenite
     for q=1:length(A_P)
         A_P(q,6)=0;
         A_P(q,12)=0;
         A_P(q,14)=0;
        for p=1:length(N_p(:,1))
         if N_p(p,6)>=0 && N_p(p,9)==1 %%% nucleated
             for r=1:length(N_PR{p}(:,1))
                 if N_PR{p}(r,1)==q
                    A_P(q,6)=A_P(q,6)+1; % update the counts of the nucleated ferrite in each austenite
                    A_P(q,12)=A_P(q,12)+N_p(p,10)/length(N_PR{p}(:,1)); % update the ferrite volume in each austenite [um3]
                    if DiffInfo{p}(i,r,7)==1
                        A_P(q,14)=A_P(q,14)+1; % counts of soft impingement
%                         A_P(q,10)=A_P(q,10)+DiffInfo{p}(i,r,1); % C content [mol%]
                    end
                 end
             end
         end
        end
        A_P(q,13)=A_P(q,12)/A_P(q,5); % ferrite fraction
        if A_P(q,13)>1
            A_P(q,13)=1;
        end
     end
     AP_track{i}=A_P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% update details of ferrite
    for s=1:length(site)
        N_p(s,15)=mean(A_P(N_PR{s}(:,1),10)); % average remote C content in austenite [wt.%]
        N_p(s,16)=(100-N_p(s,15))*Ux/(N_p(s,15)+(100-N_p(s,15))*Ux+(100-N_p(s,15))*(1-Ux))*100; % average remote Mn content in austenite [wt.%]
        N_p(s,17)=wC_F_eq(i); % average remote C content in ferrite [wt.%]
        N_p(s,18)=(100-N_p(s,17))*Ux/(N_p(s,17)+(100-N_p(s,17))*Ux+(100-N_p(s,17))*(1-Ux))*100; % average remote Mn content in ferrite [wt.%]
    end
    
     for Ia=1:length(Nucleated{i}(:,1)) %to exlude the nucleation position which has been contained in the ferrite and calculate the C concentration for each site
         a=Nucleated{i}(Ia,1);
        for b=1:length(site)
            if a~=b
                d(a,b)=sqrt((N_p(b,2)-N_p(a,2))^2+(N_p(b,3)-N_p(a,3))^2+(N_p(b,4)-N_p(a,4))^2);
                d_MXY(a,b)=sqrt((N_p(b,2)-position_MXY(a,1))^2+(N_p(b,3)-position_MXY(a,2))^2+(N_p(b,4)-position_MXY(a,3))^2);
                d_MXZ(a,b)=sqrt((N_p(b,2)-position_MXZ(a,1))^2+(N_p(b,3)-position_MXZ(a,2))^2+(N_p(b,4)-position_MXZ(a,3))^2);
                d_MYZ(a,b)=sqrt((N_p(b,2)-position_MYZ(a,1))^2+(N_p(b,3)-position_MYZ(a,2))^2+(N_p(b,4)-position_MYZ(a,3))^2);
                d_MO(a,b)=sqrt((N_p(b,2)-position_MO(a,1))^2+(N_p(b,3)-position_MO(a,2))^2+(N_p(b,4)-position_MO(a,3))^2);
                Rab(a,b)=N_p(a,7);
                if (d(a,b)<=N_p(a,7)&&N_p(b,9)==0)||(d_MXY(a,b)<=N_p(a,7)&&N_p(b,9)==0)|| ...
                        (d_MXZ(a,b)<=N_p(a,7)&&N_p(b,9)==0)||(d_MYZ(a,b)<=N_p(a,7)&&N_p(b,9)==0)|| ...
                        (d_MO(a,b)<=N_p(a,7)&&N_p(b,9)==0)
                   N_p(b,5)=0; % in the range of Rbcc+shield distance
                end
                for c=1:length(N_PR{a}(:,1))
                    if (d(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)||(d_MXY(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)|| ...
                        (d_MXZ(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)||(d_MYZ(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)|| ...
                        (d_MO(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)
                    Lmin=min([d(a,b) d_MXY(a,b) d_MXZ(a,b) d_MYZ(a,b) d_MO(a,b)])-N_p(a,7); % minimum distance from the interface [um]
                    N_p(b,22)=DiffInfo{a}(i,c,1)+(DiffInfo{a}(i,c,2)-DiffInfo{a}(i,c,1))*(1-Lmin/DiffInfo{a}(i,c,3))^2; % carbon concentration [mol%]
                    N_p(b,23)=100*(100-N_p(b,22))*Ux*M_Mn/(N_p(b,22)*M_C+(100-N_p(b,22))*Ux*M_Mn+(100-N_p(b,22))*(1-Ux)*M_Fe); % Mn concentent [wt.%]
                    N_p(b,22)=100*N_p(b,22)*M_C/(N_p(b,22)*M_C+(100-N_p(b,22))*Ux*M_Mn+(100-N_p(b,22))*(1-Ux)*M_Fe); % carbon concentration [wt.%]
                    end
                end
            else
               N_p(b,22)=wC_F_eq(i); % C content [wt.%]
               N_p(b,23)=100*(100-N_p(b,22))*Ux/(N_p(b,22)+(100-N_p(b,22))*Ux+(100-N_p(b,22))*(1-Ux)); % Mn content [wt.%]
            end
     end
        %%%%%%%%%%%%
        N_p(a,29)=0;
        for m=1:length(N_PR{a}(:,1))
            N_p(a,29)=N_p(a,29)+A_P(N_PR{a}(m,1),5)*(1-A_P(N_PR{a}(m,1),13)); % volume that it can maximum grow [um^3]
        end
        N_p(a,29)=N_p(a,29)+N_p(a,10); % the rest surrounding austenite volume+its own volume [um^3]
     end
    
     % minimum distance from a nucleus[um]
    for b=1:length(site)
        if N_p(b,9)==0
        N_p(b,28)=min([min(nonzeros(d(:,b))) min(nonzeros(d_MXY(:,b))) min(nonzeros(d_MXZ(:,b))) min(nonzeros(d_MYZ(:,b))) min(nonzeros(d_MO(:,b)))]);
        N_p(b,12)=0; %ratio of volume fraction to equilibrium volume fraction predicted by phase diagram
        V_Fsurrounding=0;
        if N_p(b,5)==1
             for r=1:length(N_PR{b}(:,1))
                 for rr=1:length(A_PR{N_PR{b}(r,1)})
                     if ~isempty(find(Nucleated{i}(:,1)==A_PR{N_PR{b}(r,1)}(rr))) && A_PR{N_PR{b}(r,1)}(rr)~=b
                        V_Fsurrounding=V_Fsurrounding+N_p(A_PR{N_PR{b}(r,1)}(rr),10)/length(N_PR{A_PR{N_PR{b}(r,1)}(rr)}(:,1));
                     end
                 end
             end
             V_surround(i,b)=sum(A_P(N_PR{b}(:,1),5))-V_Fsurrounding-4/3*pi*N_p(b,8).^3; % surrounding untransformed austenite
             if V_surround(i,b)<0
                 V_surround(i,b)=0;
             end
             N_p(b,32)=min([N_p(b,32) (3*(V_surround(i,b)+N_p(b,10))/(4*pi)).^(1/3)]); % May 10, 2020
        end
        else
            N_p(b,28)=0;
            N_p(b,12)=N_p(b,10)/(4/3*pi*N_p(b,32)^3);
            
            if N_p(b,12)>F_eq(i)
                N_p(b,12)=F_eq(i);
            end
            if ~isempty(Nucleated{i}) && length(Nucleated{i}(:,1))>=min(nonzeros(sumN))
             V_Fsurrounding=0;
             for r=1:length(N_PR{b}(:,1))
                 for rr=1:length(A_PR{N_PR{b}(r,1)})
                     if ~isempty(find(Nucleated{i}(:,1)==A_PR{N_PR{b}(r,1)}(rr))) && A_PR{N_PR{b}(r,1)}(rr)~=b
                        V_Fsurrounding=V_Fsurrounding+N_p(A_PR{N_PR{b}(r,1)}(rr),10)/length(N_PR{A_PR{N_PR{b}(r,1)}(rr)}(:,1));
                     end
                 end
             end
             V_surround(i,b)=sum(A_P(N_PR{b}(:,1),5))-V_Fsurrounding-4/3*pi*N_p(b,8).^3; % surrounding untransformed austenite
             if V_surround(i,b)<0
                 V_surround(i,b)=0;
             end
            if N_p(b,25)==0 && N_p(b,26)==0
                minD_d=d(:,b)-Rab(:,b);
                minD_dMXY=d_MXY(:,b)-Rab(:,b);
                minD_dMXZ=d_MXZ(:,b)-Rab(:,b);
                minD_dMYZ=d_MYZ(:,b)-Rab(:,b);
                minD_dMO=d_MO(:,b)-Rab(:,b);
                minD=min([min(nonzeros(minD_d)) min(nonzeros(minD_dMXY)) min(nonzeros(minD_dMXZ)) ...
                    min(nonzeros(minD_dMYZ)) min(nonzeros(minD_dMO))]); % update the diffusion distance [um]
             
            if ~isempty(find(minD_d==minD))
                N_p(b,31)=find(minD_d==minD);
            else if ~isempty(find(minD_dMXY==minD))
                    N_p(b,31)=find(minD_dMXY==minD);
                else if ~isempty(find(minD_dMXZ==minD))
                        N_p(b,31)=find(minD_dMXZ==minD);
                    else if ~isempty(find(minD_dMYZ==minD))
                            N_p(b,31)=find(minD_dMYZ==minD);
                        else if ~isempty(find(minD_dMO==minD))
                                N_p(b,31)=find(minD_dMO==minD);
                            end
                        end
                    end
                end
            end
             if V_Fsurrounding==0 && N_p(b,31)~=0
                   N_p(b,32)=min([N_p(b,32) (3*(V_surround(i,b)+N_p(b,10))/(4*pi)).^(1/3)]); % May 10, 2020
             else
                   N_p(b,32)=min([N_p(b,32) (3*(V_surround(i,b)+N_p(b,10))/(4*pi)).^(1/3)]); % May 10, 2020
             end
            else if N_p(b,25)==1 && N_p(b,26)==0 && Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),25)==1
                    N_p(b,27)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),27);
                    N_p(b,32)=min([N_p(b,32) (3*(V_surround(i,b)+N_p(b,10))/(4*pi)).^(1/3)]); % May 10, 2020
                else if N_p(b,26)==1
                    N_p(b,27)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),27);
                    N_p(b,32)=min([N_p(b,32) (3*(V_surround(i,b)+N_p(b,10))/(4*pi)).^(1/3)]); % May 10, 2020
                    end
                end
                N_p(b,31)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),31);
            end
               N_p(b,32)=mean(N_PR{b}(:,3)); % average vertice and austenite center distance [um], resume to use on July 17,2020
               N_pD(b,25:32)=N_p(b,25:32);
            end
        end
    end
    
    % Update info for nucleated ferrite
    for b=1:length(Nucleated{i}(:,1))
        for a=1:length(N_p(:,1))
            if N_p(a,1)==Nucleated{i}(b,1)
               Nucleated{i}(b,5:end)=N_p(a,5:end);
               DiffInfo{Nucleated{i}(b,1)}(i,:,4)=Nucleated{i}(b,7);
            end
        end
    end
    end
    
%%%%%%%%% overall statistics
    if ~isempty(Nucleated{i})
         delta(i)=mean(Nucleated{i}(:,8));
         delta_S(i)=sum(Nucleated{i}(:,8).^4)/sum(Nucleated{i}(:,8).^3);
    else
        delta(i)=0;
        delta_S(i)=0;
    end
    F(i)=volume/Lb^3; %volume fraction
    trans=N_p(:,8); % An intermediate array to install the effective radius
    R_sd(i)=std(trans(trans~=0));% Standard deviation of the effective radius 

    %%%%%%%%%%% the extended volume fraction
    if isempty(Nucleated{i})
        Fextend(i)=0;
    else
        Fextend(i)=sum(4/3*pi*Nucleated{i}(:,7).^3)/Lb^3;
    end
  
%%%%%%%%%%%%%% Restore the state of austenite grains
    A_state{i}=A_P;
    A_state{i}(:,length(A_P(1,:))+1)=(3*A_state{i}(:,5).*(1-A_state{i}(:,13))./(4*pi)).^(1/3);% effective austenite radius [um]

%%%%%%%%%%%%%%%%%%%%%%%%%
    if T(i)<=Temp_eq(5,2)||Timer(i)>tcr(end)
        stop_cycle=1;
    end
    t=t+dt(i);
end
%%%% save outputs
save('Simulation_result.mat'); 
% dlmwrite('myfile.txt',[Timer(1:end-1)' T(1:end-1)' N_eff(1:end-1)'/(Lb*1e-6)^3 F' ...
%     F_eq(1:end-1)' delta' R_sd'],'delimiter',' ');
dlmwrite('myfile.txt',[Timer' T' N_eff'/(Lb*1e-6)^3 F' ...
    F_eq' delta' R_sd'],'delimiter',' ');

figure;
subplot(2,3,1);
plot(Timer,T,'ro-');
xlabel('Time (s)');
ylabel('Temperature (K)');
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
grid on;box on;
subplot(2,3,2);
plot(Timer,N_eff,'ro-');
xlabel('Time (s)');
ylabel('N_{\alpha}');
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
grid on;box on;
subplot(2,3,3);
plot(Timer,F,'ro-',Timer,F_eq,'kx-');
legend('Modelling','Equilibrium under PE');
xlabel('Time (s)');
ylabel('f_{\alpha}');
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
grid on;box on;
subplot(2,3,4);
plot(Timer,delta,'ro-',Timer,R_sd,'bx-');
legend('<R_{\alpha}>','\sigma_{R_{\alpha}}');
xlabel('Time (s)');
ylabel('<R_{\alpha}> (\mum)');
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
grid on;box on;
subplot(2,3,5);
plot(T,F,'ro-');
xlabel('Temperature (K)');
ylabel('f_{\alpha}');
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
grid on;box on;
subplot(2,3,6);
plot(Timer,v_t(:,1),'ro-');
xlabel('Time (s)');
ylabel('V_{int} for grain No.1 (\mum/s)');
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
grid on;box on;

%%%%%%%%%% visulize results
figure('Name','Austenite geometry');
subplot(2,2,1);
hold all;
V.plot('alpha', 0.1); % Adjust the transparancy
plot3(A(:,1),A(:,2),A(:,3),'k+','LineWidth',2);% Plot the centroids
plot3(site(:,1),site(:,2),site(:,3),'r.','MarkerSize',20);% Plot the corners
plot3(Nucleated{i-1}(:,2),Nucleated{i-1}(:,3),Nucleated{i-1}(:,4),'b*','LineWidth',2);
grid off;
xlabel('x (micron)','FontSize',12);
ylabel('y (micron)','FontSize',12);
zlabel('z (micron)','FontSize',12);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
box on;
hold off;
% legend('+ centroids of austenite grains','. potential ferrite nucleation sites','* randomly selected nucleation sites');

subplot(2,2,2);
% plot(Timer',N_eff');
plot(T'-273,N_eff');
title('Nucleation at time t');
% xlabel('Time (seconds)');
xlabel('Temperature (^{o}C)');
ylabel('Numble of nuclei');

subplot(2,2,3);
% Plot the ferrite particles
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,7)>0
    ssphere(N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,7),Lb);
    hold on;
    ssphere(position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,7),Lb);
    hold on;
    ssphere(position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,7),Lb);
    hold on;
    ssphere(position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,7),Lb);
    hold on;
    ssphere(position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,7),Lb);
    hold on;
    title('3D visualization-Original');
    hold on;
    end
end
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,7)>0
        sphere0=[N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;         
        sphere0=[position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;
        sphere0=[position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
            hold on;
    end
end
% Plot the original austenite grains
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
V.plot('alpha', 0.4); % Adjust the transparancy
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
grid off;
hold off;
box on;

subplot(2,2,4);
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,8)>0
    ssphere(N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,8),Lb);
        hold on;
    ssphere(position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,8),Lb);
       hold on;
    ssphere(position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,8),Lb);
        hold on;
    ssphere(position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,8),Lb);
        hold on;
    ssphere(position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,8),Lb);
        hold on;
   % axis([0 126 0 126 0 126]);
    title('3D visualization-after substract the overlay');
    hold on;
    end
end
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,7)>0
        sphere0=[N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;         
        sphere0=[position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;
        sphere0=[position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
            hold on;
    end
end
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
V.plot('alpha', 0.05); % Adjust the transparancy
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
grid off;
hold off;
box on;

% % taking a slice
% S = V.slice(2,0);
% S.plot('color','red','alpha',0.2,'linestyle','-','linewidth',5);

%%%%%% Nuclei density and radius
figure('Name','Nuclei density and radius');
subplot(2,2,1);
for b=1:length(Nucleated)
   if isempty(Nucleated{b})
       R_sd0(b)=0;
   else
       R_sd0(b)=std(Nucleated{b}(:,8));
   end
end
plot(T'-273,delta','-bo',T'-273,R_sd','r.-',T'-273,R_sd0','k.-');
legend('<r>','\sigma','\sigma0');
title('Average radius of ferrite grains');
% xlabel('Time (seconds)');
xlabel('Temperature (^{o}C)');
ylabel('Average radius (\mum)');

subplot(2,2,2);
% plot(Timer',F','-bo',Timer',F_eq','r');
plot(T'-273,F','-bo',T'-273,F_eq','r');
legend('modelling','equilibrium');
title('Volume fraction at time t');
% xlabel('Time (seconds)');
xlabel('Temperature (^{o}C)');
ylabel('Volume fraction of ferrite phase');

subplot(2,2,3);
plot(F,N_eff/Lb^3,'ro-');
xlabel('Ferrite fraction');
ylabel('\rho_{\alpha} (\mum^{-3})');

subplot(2,2,4);
plot(F,delta','-bo',F,R_sd','r.-');
legend('<r>','\sigma');
xlabel('Ferrite fraction');
ylabel('\delta_{S},\sigma_{P} (\mum)');

% TrackNo=fix(length(Nucleated{i-1})*rand(1,1));
TrackNo=1;
figure('Name', 'Track status for a particular ferrite');
subplot(2,3,1);
plot(Timer,DiffInfo{TrackNo}(:,1,4),'ko-');
hold all;
for pp=1:length(Timer)
    if ~isempty(Nucleated{pp}) && ~isempty(find(Nucleated{pp}(:,1)==TrackNo))
        plot(Timer(pp),Nucleated{pp}(find(Nucleated{pp}(:,1)==TrackNo),8),'ro-');
        yy(pp)=Nucleated{pp}(find(Nucleated{pp}(:,1)==TrackNo),8);
    end
end
hold off;
xlabel('Time (s)');
ylabel('Ferrite radius (\mum)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

for b=1:length(Timer)
    if ~isempty(Nucleated{b}) && ismember([TrackNo],Nucleated{b})
       GrowSpace(b)=Nucleated{b}(find(Nucleated{b}(:,1)==TrackNo),32);
    end
end
subplot(2,3,2);
plot(Timer,DiffInfo{TrackNo}(:,1,3),'ro-',Timer,GrowSpace,'b+');
xlabel('Time (s)');
ylabel('Length (\mum)');
legend('Diffusion length','Growth length');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,3);
plot(Timer,DiffInfo{TrackNo}(:,1,7),'ro-');
xlabel('Time (s)');
ylabel('Flag of soft impingement');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,4);
plot(Timer,DiffInfo{TrackNo}(:,1,1),'ro-',Timer,DiffInfo{TrackNo}(:,1,2),'b+-',Timer,DiffInfo{TrackNo}(:,1,6),'kd-');
legend('C_{matrix}','C^{i}','C^{eq}');
xlabel('Time (s)');
ylabel('C content (mol%)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,5);
plot(Timer,v_t(:,TrackNo),'ro-');
xlabel('Time (s)');
ylabel('V_{interface} (\mum/s)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,6);
plot(Timer,Gchem(:,TrackNo),'ro-',Timer,Gfriction(:,TrackNo),'b+-', ...
    Timer,Gdiff(:,TrackNo),'md-',Timer,Gfriction(:,TrackNo)+Gdiff(:,TrackNo),'k-');
legend('Gm_{chem}','Gm_{friction}','Gm_{diff}','Gm_{total dissipation}');
xlabel('Time (s)');
ylabel('Gibbs energy (J/mol)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);
for b=1:length(Timer)
    if ~isempty(Nucleated{b}) && ismember([TrackNo],Nucleated{b})
       GrowSpace(b)=Nucleated{b}(TrackNo,32);
       SoftSpace(b)=Nucleated{b}(TrackNo,27);
    end
end
figure;
plot(DiffInfo{TrackNo}(:,1,4),GrowSpace,'ro-',DiffInfo{TrackNo}(:,1,4),SoftSpace,'b+-');
xlabel('Radius (\mum)');
ylabel('Length (\mum)');
legend('GrowSpace','SoftSpace');
set(gca,'fontsize',12);
set(gca,'linewidth',2);
