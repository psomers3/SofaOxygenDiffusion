
%% Variablendefintion
cellarray=cells.empty;
zyklen=0;
positionen=double.empty;
forces=double.empty;
timer=0;
t_hypo=7; % mmHg
t_dead= 0.7; % mmHg
alpha_n=0.3;
alpha_h=0.06;
Zellvolumen=4/3*(5)^3*pi;
timestep=1;
maxradius=5e-6;

entwicklung=[1,1,1,0];
extrapolation=double.empty;

%% Mesh für FEM generieren

model = createpde(1);
gm = multicuboid(500,500,500,'ZOffset',-250); %mit Faktor 10^6
model.Geometry=gm;
generateMesh(model,'Hmax', 10*5,'Hmin', 10, 'GeometricOrder','linear');
TR=triangulation(model.Mesh.Elements',model.Mesh.Nodes');

saveTetraMesh('TetraMesh',model.Mesh.Nodes',model.Mesh.Elements')
oxygen=ones(length(model.Mesh.Elements),1)*30;

%% erste Zelle setzen
startcell=cells;
startcell.pos=[0 0 0];
startcell.radius=5e-6;
cellarray=addtoarray(startcell,cellarray);
figure
