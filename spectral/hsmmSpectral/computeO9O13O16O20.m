function  computeO9O13O16O20(A1, A, D1, D, O)


[o x] = size(O);
d = size(D1,1);


%O1O5D6X6
M = zeros(x, x, x, d);



% for ind_d6 = 1:d
% for ind_x6 = 1:x
% for ind_x1 = 1:x
% for ind_x5 = 1:x
%    
%     for ind_d1 = 1:d
%     for ind_d2 = 1:d
%     for ind_d3 = 1:d    
%     for ind_d4 = 1:d
%     for ind_d5 = 1:d
%     for ind_x2 = 1:x    
%     for ind_x3 = 1:x
%     for ind_x4 = 1:x
%         
%         
%         M(ind_x1, ind_x5, ind_x6, ind_d6) = M(ind_x1, ind_x5, ind_x6, ind_d6) + ...
%             A(ind_x6, ind_x5, ind_d5) * D(ind_d6, ind_x6, ind_d5) *  A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4) * A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2) * A(ind_x2, ind_x1, ind_d1) * D(ind_d2, ind_x2, ind_d1) *D1(ind_d1, ind_x1) * A1(ind_x1);
%         
%         
%     end    
%     end
%     end
%     end
%     end
%     end
%     end
%     end
%     
%         
% end
% end
% end
% end
% 
% 
% Mm = zeros(x*x, d*x);
% 
% for ind_d6 = 1:d
% for ind_x6 = 1:x
% for ind_x5 = 1:x
% for ind_x1 = 1:x
%     
%     Mm(ind_x1 + x*(ind_x5), ind_x6 + d*(ind_d6-1)) = M(ind_x1, ind_x5, ind_x6, ind_d6);
%     
% end
% end
% end
% end



%X7X11D6X6
X = zeros(o,o,d,x);

for ind_o7 = 1:o
for ind_o11 = 1:o    
for ind_x7 = 1:x
for ind_x11 = 1:x
for ind_d6 = 1:d
for ind_x6 = 1:x
    
    for ind_d10 = 1:d
    for ind_d9 = 1:d    
    for ind_d8 = 1:d    
    for ind_d7 = 1:d    
    for ind_x10 = 1:x
    for ind_x9 = 1:x    
    for ind_x8 = 1:x
        
        X(ind_o7, ind_o11, ind_d6, ind_x6) = X(ind_o7, ind_o11, ind_d6, ind_x6) + ...
            O(ind_o7, ind_x7) * O(ind_o11, ind_x11) * A(ind_x11, ind_x10, ind_d10) * A(ind_x10, ind_x9, ind_d9) * D(ind_d10, ind_x10, ind_d9) * A(ind_x9, ind_x8, ind_d8) * D(ind_d9, ind_x9, ind_d8) * A(ind_x8, ind_x7, ind_d7) * D(ind_d8, ind_x8, ind_d7) * A(ind_x7, ind_x6, ind_d6) * D(ind_d7, ind_x7, ind_d6);
        
        
    end
    end
    end
    end
    end
    end
    end
end
end
end
end
end
end



Xm = zeros(o*o, x*d);

for ind_o11 = 1:o
for ind_o7 = 1:o
for ind_d6 = 1:d
for ind_x6 = 1:x
    
    Xm(ind_o7 + o*(ind_o11-1), ind_x6 + d*(ind_d6-1)) = X(ind_o7, ind_o11, ind_d6, ind_x6);
    
end
end
end
end



Y = zeros(o,o);

for ind_o2 = 1:o
for ind_o3 = 1:o
for ind_x2 = 1:x
    
for ind_x3 = 1:x
for ind_x1 = 1:x
for ind_d2 = 1:d
for ind_d1 = 1:d
    
    Y(ind_o3, ind_o2) = Y(ind_o3, ind_o2) + ...
        O(ind_o2, ind_x2) * O(ind_o3, ind_x3) * A(ind_x3, ind_x2, ind_d2) * A(ind_x2, ind_x2, ind_d2) * D(ind_d2, ind_x2, ind_d2) * D1(ind_d1, ind_x1) * A1(ind_x1);
    
end
end
end
end
end
end
end



Y;









    