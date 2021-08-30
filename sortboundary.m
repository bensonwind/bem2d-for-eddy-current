function [new_edge] = sortboundary(edge,dir)
    %UNTITLED3 此处显示有关此函数的摘要
    %   此处显示详细说明
    [n,m] = size(edge);
    new_edge = zeros(n,m);
    flip_flag = ones(n,1);
    currentP_idx = [];
    for i = 1:n
       if i == 1
           new_edge(i,:) = edge(i,:);
           if dir==true
               currentP_idx = new_edge(i,1);
               flip_flag(i) = -1;
           else
               currentP_idx = new_edge(i,2);
           end
       else
            idx = find( edge(:,1) == currentP_idx | edge(:,2) == currentP_idx);
            if length(idx)~=2
                return
            end
            
            if new_edge(i-1,:) == edge(idx(1),:)
                new_edge(i,:) = edge(idx(2),:);
            else
                new_edge(i,:) = edge(idx(1),:);
            end
            if new_edge(i,1)==currentP_idx
                currentP_idx = new_edge(i,2);
            else
                currentP_idx = new_edge(i,1);
                flip_flag(i) = -1;
            end
       end
    end
    
    for i = 1:n
        if flip_flag(i) == -1
            new_edge(i,:) = fliplr(new_edge(i,:));
        end
    end
end
