function [cellarray] = celldivide11112(c,maxradius,cellarray)

    %% Mutterzelle bearbeiten
    cellarray(c).radius=0.5*cellarray(c).radius;
    cellarray(c).teilbar=false;
    mutterpos=cellarray(c).pos;
    dist=cellarray(c).radius;
    theta=rand*180;
    phi=rand*360;
    x=dist*sin(theta)*cos(phi);
    y=dist*sin(theta)*sin(phi);
    z=dist*cos(theta);

    cellarray(c).pos=[mutterpos(1)+x mutterpos(2)+y mutterpos(3)+z];

    %% Erzeugung neuer Zelle
    newcell=cellarray(c);
    newcell.pos=[mutterpos(1)-x mutterpos(2)-y mutterpos(3)-z];
    newcell.oxygen=cellarray(c).oxygen;
    cellarray =addtoarray(newcell, cellarray);
end

