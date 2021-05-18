function plotSphere(cellarray,maxradius)

for i=1:length(cellarray)
    rad=cellarray(i).radius;
    pos=cellarray(i).pos;
    [X Y Z]=sphere(16);
    X2 = X * rad;
    Y2 = Y * rad;
    Z2 = Z * rad;

    plot=surf(X2+pos(1),Y2+pos(2),Z2+pos(3));
    axis equal
    if cellarray(i).phenotype==0
        set(plot,'FaceColor',[0 1 0], 'FaceAlpha',0.5,'EdgeAlpha', 0.5); 
    end
    if cellarray(i).phenotype==1
        set(plot,'FaceColor',[0 0 1], 'FaceAlpha',0.5,'EdgeAlpha', 0.5); 
    end
    if cellarray(i).phenotype==2
        set(plot,'FaceColor',[1 0 0], 'FaceAlpha',0.5,'EdgeAlpha', 0.5); 
    end
    hold on
end
hold off
end

