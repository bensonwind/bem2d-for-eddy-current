function [p,f] = nastranmat_comsol(FileName)
%   Read a mesh in NASTRAN ASCII format. 

    fileID = fopen(FileName, 'r');
    p = []; f = [];
    tline = fgetl(fileID);
    while ischar(tline)
        Nodes = strfind(tline, 'GRID');
        if ~isempty(Nodes)
            tmp = sscanf(tline, 'GRID%d%f%f');
            p = [p; tmp([2 3]).'];
        end
        
        Tri  = strfind(tline, 'CTRIA3');
        if ~isempty(Tri)
            tmp = sscanf(tline, 'CTRIA3%d%d%d%d%d');
            f = [f; tmp([3 4 5]).'];
        end
        
        tline = fgetl(fileID);
    end
  
% Restore outer normals
%     normals = zeros(size(t)); 
%     for m = 1:length(t)
%         Vertexes        = p(t(m, 1:3)', :)';
%         r1              = Vertexes(:, 1);
%         r2              = Vertexes(:, 2);
%         r3              = Vertexes(:, 3);
%         tempv           = cross(r2-r1, r3-r1);  %   clockwise normal
%         temps           = sqrt(tempv(1)^2 + tempv(2)^2 + tempv(3)^2);
%         normals(m, :)   = tempv/temps; 
%     end
    
    fclose(fileID);

end





