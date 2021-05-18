function saveTetraMesh(filename,nodes,elements)
% first column is 1..#nodes
nodes(:,2:4)=nodes(:,1:3);
nodes(:,1)=1:size(nodes,1);
elements(:,2:5)=elements(:,1:4);
elements(:,1)=1:size(elements,1);
fileID = fopen(filename,'w');
fprintf(fileID,'$NOD\n');
fprintf(fileID,'%d\n',size(nodes,1));
fprintf(fileID,'%d %3.3f %3.3f %3.3f\n',nodes');
fprintf(fileID,'$ENDNOD\n');
fprintf(fileID,'$ELM\n');
fprintf(fileID,'%d\n',size(elements,1));
fprintf(fileID,'%d 4 1 1 4 %d %d %d %d\n',elements(:,1:5)');
fprintf(fileID,'$ENDELM\n');
fclose(fileID);
end