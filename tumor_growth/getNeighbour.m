function [cellarray] = getNeighbour(c,cellarray)
%sets all surrounding cells as neigbour of the cell c
cellarray(c).nachbarn=cells.empty;
        for vglc=1:length(cellarray)
            cellkopie=cellarray(vglc);
            cellkopie.nachbarn=[];
            dist=norm(cellarray(c).pos-cellarray(vglc).pos);
            if 0<dist & dist<10e-6
                cellarray(c).nachbarn=addtoarray(cellkopie,cellarray(c).nachbarn);
            end
        end
end

