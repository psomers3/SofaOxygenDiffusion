%% Parameter aktualtiserung
TR=triangulation(model.Mesh.Elements',model.Mesh.Nodes');
timer=0;
numcells=length(cellarray);
zyklenalt=zyklen;

for c=1:length(cellarray)
    [ID,~]=pointLocation(TR,cellarray(c).pos*10^6);
    cellarray(c).oxygen=oxygen(ID);
 
end
for i=1:length(cellarray)
    if cellarray(i).oxygen<7 && cellarray(i).phenotype==0
        cellarray(i).phenotype=1;%hypoxic
    end

    if cellarray(i).oxygen<0.7
        cellarray(i).phenotype=2;%necrotic
    end
end
%% Simulation
while length(cellarray)<numcells+5 && zyklen<zyklenalt+1000
    tic
    zyklen=zyklen+1;
    %% Zellwachstum für jede Zell
    for c = 1:length(cellarray)
        if cellarray(c).phenotype==0
            cellarray(c)=cellarray(c).grow;
        end
    end
    
    %% Infos über Nachbarn für jede Zelle des cellarray
    for c=1:length(cellarray)
        cellarray=getNeighbour(c,cellarray);
        if length(cellarray(c).nachbarn)>16
            cellarray(c).teilbar=false;
        end

    end
    
    
    %% Zellbewegung
    for c=1:length(cellarray)
        cellarray(c)=cellarray(c).move3D;
    end
    
    
    %% Zellteilung
    cellarraykopie=cellarray;
    for c=1:length(cellarraykopie)
        if cellarray(c).teilbar==true 
            if cellarray(c).phenotype==0
                if rand<1/1000
                   cellarray=celldivide3D(c,maxradius,cellarray);
                end  
            end
        end
    end
    

%% Kräfte für jede Zelle berechnen
    for c=1:length(cellarray)
        cellarray(c).F_rep=0;
        cellarray(c).F_adh=0;
        cellarray(c).F=[0 0 0];
            for vglc=1:length(cellarray(c).nachbarn)
                [F_rep,F_adh,F]=computeForces(cellarray(c),cellarray(c).nachbarn(vglc));
                cellarray(c).F_rep=cellarray(c).F_rep+F_rep;
                cellarray(c).F_adh=cellarray(c).F_adh+F_adh;
                cellarray(c).F=cellarray(c).F+F;
            end
        if cellarray(c).F_rep>4.1954e-05
            cellarray(c).teilbar=false;
        end
    end
%% phenotypwechsel
    if cellarray(i).phenotype==1 && cellarray(i).oxygen>7
        if rand<1/(24*60)
            cellarray(i).phenotype=0;   %back to normoxic
        end
    end
      

    %% Auswertungsstuff
    extrapolation(zyklen)=length(cellarray);
    timer=timer+toc;
    if mod(zyklen,10)==0
        disp("Anzahl Zellen: "+length(cellarray)+" in Timestep "+zyklen+" in "+timer+" s");
    end
end


%% Wie viele Zellen sind in welchem Tetraeder?
for c=1:length(cellarray)
    [ID,~]=pointLocation(TR,cellarray(c).pos*10^6);
    cellarray(c).ID=ID;
end

nbHypoCellsInTetra=zeros(length(TR.ConnectivityList),1);
nbNormCellsInTetra=zeros(length(TR.ConnectivityList),1);
for c=1:length(cellarray)
    if cellarray(c).phenotype==0
       nbNormCellsInTetra(cellarray(c).ID)=nbNormCellsInTetra(cellarray(c).ID)+1;
    end
    if cellarray(c).phenotype==1
       nbHypoCellsInTetra(cellarray(c).ID)=nbHypoCellsInTetra(cellarray(c).ID)+1;
    end
end
 %% TetraUptakeCoefficient berechnen
[V,VolumeEveryTetra] = volume(model.Mesh);

TetraederNorm = [TR.ConnectivityList,VolumeEveryTetra',nbNormCellsInTetra, nbNormCellsInTetra*Zellvolumen]; 
densitieNorm = TetraederNorm(:,7)./TetraederNorm(:,5); 


TetraederHypo = [TR.ConnectivityList,VolumeEveryTetra',nbHypoCellsInTetra, nbHypoCellsInTetra*Zellvolumen]; 
densitieHypo = TetraederHypo(:,7)./TetraederHypo(:,5); 

tetraUptakeCoefficient=alpha_n*densitieNorm+alpha_h*densitieHypo;

%% Asuwertungsstuff
entwicklung(length(entwicklung(:,1))+1,:)=[zyklen,length(cellarray),sum(nbNormCellsInTetra),sum(nbHypoCellsInTetra)];
disp("Anzahl normotische Zellen: "+sum(nbNormCellsInTetra))
disp("Anzahl hypotische Zellen: "+sum(nbHypoCellsInTetra))
disp("Anzahl tote Zellen: "+(length(cellarray)-sum(nbNormCellsInTetra)-sum(nbHypoCellsInTetra)))
clf
plotSphereOxygen(cellarray,5e-6)
%% Datensicherung
save('workspace')
