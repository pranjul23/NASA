%function to test P(Oi,Oj|x2,d2)

function computeOiOjX2D2(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


% M = zeros(o, o, x, d);
% 
% %O2O3X2D2
% for ind_o2 = 1:o
%     for ind_o3 = 1:o
%         for ind_d2 = 1:d
%             for ind_x2 = 1:x
%                 
%                  for ind_d3 = 1:d
% %                 for ind_d4 = 1:d
% %                 for ind_d5 = 1:d  
%                 for ind_x3 = 1:x
%                  for ind_x4 = 1:x
% %                 for ind_x5 = 1:x
% %                 for ind_x6 = 1:x    
%                             
%                     M(ind_o2, ind_o3, ind_x2, ind_d2) = M(ind_o2, ind_o3, ind_x2, ind_d2) + ...
%                         O(ind_o2, ind_x2) * O(ind_o3, ind_x4)  * A(ind_x4, ind_x3, ind_d3) * A(ind_x3, ind_x2, ind_d2) * D(ind_d3, ind_x3, ind_d2);
%                         
%                             
% %                 end
% %                  end
%                  end
%                 end
%                 end
% %                 end
% %                 end
%                 
%             end
%         end
%     end
% end
% 
% 
% %normalize
% for ind_x2 = 1:x
% for ind_d2 = 1:d
%      M(:, :, ind_x2, ind_d2) = M(:, :, ind_x2, ind_d2) / sum(sum(M(:, :, ind_x2, ind_d2)));
% end
% end
% 
% F = zeros(o*o, x*d);
% for ind_o3 = 1:o
% for ind_o2 = 1:o
% for ind_d2 = 1:d
% for ind_x2 = 1:x 
%     
%     F(ind_o2 + o*(ind_o3-1), ind_x2 + x*(ind_d2-1)) = M(ind_o2, ind_o3, ind_x2, ind_d2);
%     
% end
% end
% end
% end




M2 = zeros(o, o, x, d);
X = zeros(x,x,x,d);
W = zeros(x,x,x,d);
Y = zeros(d,d,x,d);

D4 = zeros(d,x,d);
X4 = zeros(x,x,d);

Tr = zeros(x,x,x,x,d);

%O2O5X2D2
for ind_o5 = 1:o
    for ind_o2 = 1:o
        for ind_d2 = 1:d
            for ind_x2 = 1:x
                
                for ind_d3 = 1:d
                for ind_d4 = 1:d
                for ind_d5 = 1:d  
                for ind_x3 = 1:x
                for ind_x4 = 1:x
                for ind_x5 = 1:x
                for ind_x6 = 1:x    
                            
                    M2(ind_o2, ind_o5, ind_x2, ind_d2) = M2(ind_o2, ind_o5, ind_x2, ind_d2) + ...
                        O(ind_o2, ind_x5) * O(ind_o5, ind_x6) * A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
                    
                    X(ind_x4, ind_x6, ind_x2, ind_d2) = X(ind_x4, ind_x6, ind_x2, ind_d2) + ...
                        A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
                    
                    W(ind_x4, ind_x5, ind_x2, ind_d2) = W(ind_x4, ind_x5, ind_x2, ind_d2) + ...
                        A(ind_x5, ind_x4, ind_d4) * A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);                                        
                    
%                     Y(ind_d4, ind_d5, ind_x2, ind_d2) = Y(ind_d4, ind_d5, ind_x2, ind_d2) + ...
%                         A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);                                                           
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


  
for ind_x4 = 1:x
    for ind_x5 = 1:x
        for ind_x6 = 1:x
            
             for ind_d2 = 1:d
            for ind_x2 = 1:x
                
                for ind_x3 = 1:x
                for ind_d3 = 1:d
                for ind_d4 = 1:d
                for ind_d5 = 1:d
                
Tr(ind_x4, ind_x5, ind_x6, ind_x2, ind_d2) = Tr( ind_x4, ind_x5, ind_x6, ind_x2, ind_d2) + ...
                        A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4)* A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);

                end
                end
                end
                end
            end
             end
        end
    end
end


% for ind_x4 = 1:x
% for ind_d2 = 1:d
% for ind_x2 = 1:x
%                 
%     for ind_d3 = 1:d
%     for ind_x3 = 1:x
%                                 
%         X4(ind_x4, ind_x2, ind_d2) = X4(ind_x4, ind_x2, ind_d2) + ...
%                             A(ind_x4, ind_x3, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
%     end
%     end
% end
% end
% end
   


% for ind_d4 = 1:d
% for ind_d2 = 1:d
% for ind_x2 = 1:x
%                 
%     for ind_d3 = 1:d
%     for ind_x3 = 1:x
%     for ind_x4 = 1:x    
%                                 
%          D4(ind_d4, ind_x2, ind_d2) = D4(ind_d4, ind_x2, ind_d2) + ...
%                         A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
%     end
%     end
%     end
% end
% end
% end
                   


%normalize
for ind_x2 = 1:x
for ind_d2 = 1:d
     M2(:, :, ind_x2, ind_d2) = M2(:, :, ind_x2, ind_d2) / sum(sum(M2(:, :, ind_x2, ind_d2)));
     X(:, :, ind_x2, ind_d2) = X(:, :, ind_x2, ind_d2) / sum(sum(X(:, :, ind_x2, ind_d2)));
     W(:, :, ind_x2, ind_d2) = W(:, :, ind_x2, ind_d2) / sum(sum(W(:, :, ind_x2, ind_d2)));
     Tr(:, :, :, ind_x2, ind_d2) = Tr(:, :, :, ind_x2, ind_d2) / sum(sum(sum(Tr(:, :, :,ind_x2, ind_d2))));
     
     %Y(:, :, ind_x2, ind_d2) = Y(:, :, ind_x2, ind_d2) / sum(sum(Y(:, :, ind_x2, ind_d2)));
%      X4(:, ind_x2, ind_d2) = X4(:, ind_x2, ind_d2)/sum(X4(:, ind_x2, ind_d2));
%      D4(:, ind_x2, ind_d2) = D4(:, ind_x2, ind_d2)/sum(D4(:, ind_x2, ind_d2));
end
end



QQ = zeros(x,x,x,d);

for ind_x4 = 1:x
for ind_x6 = 1:x
for ind_x2 = 1:x
for ind_d2 = 1:d
   
    for ind_x5 = 1:x
        QQ(ind_x4, ind_x6, ind_x2, ind_d2) = QQ(ind_x4, ind_x6, ind_x2, ind_d2) + ...
        Tr(ind_x4, ind_x5, ind_x6, ind_x2, ind_d2);
    end
    
end
end
end
end


for ind_x2 = 1:x
for ind_d2 = 1:d
    
    QQ(:, :, ind_x2, ind_d2) = QQ(:, :, ind_x2, ind_d2)/sum(sum(QQ(:, :, ind_x2, ind_d2)));
    
end
end


QQm = zeros(x*x, x*d);

for ind_x6 = 1:x
for ind_x5 = 1:x
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    QQm(ind_x5 + x*(ind_x6-1), ind_x2 + x*(ind_d2-1)) = QQ(ind_x5, ind_x6, ind_x2, ind_d2);
    
end
end
end
end


% X4X6X2D2 = zeros(x,x,x,d);
% 
% for ind_x6 = 1:x
%     for ind_x4 = 1:x
%         for ind_d2 = 1:d
%             for ind_x2 = 1:x
%                 
%                 for ind_d4 = 1:d
%                 for ind_d5 = 1:d                
%                 for ind_x5 = 1:x
%                 
%                     X4X6X2D2(ind_x4, ind_x6, ind_x2, ind_d2) = X4X6X2D2(ind_x4, ind_x6, ind_x2, ind_d2) + ...
%                         A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4) * D4(ind_d4, ind_x2, ind_d2) * X4(ind_x4, ind_x2, ind_d2);
%                 end
%                 end
%                 end
%                 
%             end            
%         end
%     end
% end
% 
% 
% 
% X4X5X2D2 = zeros(x,x,x,d);
% 
% for ind_x5 = 1:x
%     for ind_x4 = 1:x
%         for ind_d2 = 1:d
%             for ind_x2 = 1:x
%                 
%                 for ind_d4 = 1:d
%                 
%                     X4X5X2D2(ind_x4, ind_x5, ind_x2, ind_d2) = X4X5X2D2(ind_x4, ind_x5, ind_x2, ind_d2) + ...
%                         A(ind_x5, ind_x4, ind_d4) * X4(ind_x4, ind_x2, ind_d2) * D4(ind_d4, ind_x2, ind_d2);
%                 end                
%             end            
%         end
%     end
% end


F2 = zeros(o*o, x*d);
S = zeros(x*x, x*d);
SW = zeros(x*x, x*d);
H = zeros(d*d, x*d);

XX4 = zeros(x, x*d);
DD4 = zeros(d, x*d);


X4X6X2D2_m = zeros(x*x,x*d);
X4X5X2D2_m = zeros(x*x,x*d);



% for ind_x6 = 1:x
% for ind_x5 = 1:x
% for ind_d2 = 1:d
% for ind_x2 = 1:x 
%     
%     X4X6X2D2_m(ind_x5 + x*(ind_x6-1), ind_x2 + x*(ind_d2-1)) = X4X6X2D2(ind_x5, ind_x6, ind_x2, ind_d2);
%     
% end
% end
% end
% end
% 
% 
% for ind_x6 = 1:x
% for ind_x5 = 1:x
% for ind_d2 = 1:d
% for ind_x2 = 1:x 
%     
%     X4X5X2D2_m(ind_x5 + x*(ind_x6-1), ind_x2 + x*(ind_d2-1)) = X4X5X2D2(ind_x5, ind_x6, ind_x2, ind_d2);
%     
% end
% end
% end
% end


Trm = zeros(x*x*x,x*d);

for ind_x6 = 1:x
for ind_x5 = 1:x
for ind_x4 = 1:x 
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    Trm(ind_x4 + x*(ind_x5-1) + x*x*(ind_x6-1), ind_x2 + x*(ind_d2-1)) = Tr(ind_x4, ind_x5, ind_x6, ind_x2, ind_d2);
    
end
end
end
end
end


for ind_o5 = 1:o
for ind_o2 = 1:o
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    F2(ind_o2 + o*(ind_o5-1), ind_x2 + x*(ind_d2-1)) = M2(ind_o2, ind_o5, ind_x2, ind_d2);
    
end
end
end
end


for ind_x6 = 1:x
for ind_x5 = 1:x
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    S(ind_x5 + x*(ind_x6-1), ind_x2 + x*(ind_d2-1)) = X(ind_x5, ind_x6, ind_x2, ind_d2);
    
end
end
end
end

for ind_x6 = 1:x
for ind_x5 = 1:x
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    SW(ind_x5 + x*(ind_x6-1), ind_x2 + x*(ind_d2-1)) = W(ind_x5, ind_x6, ind_x2, ind_d2);
    
end
end
end
end


for ind_d5 = 1:d
for ind_d4 = 1:d
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    H(ind_d4 + d*(ind_d5-1), ind_x2 + x*(ind_d2-1)) = Y(ind_d4, ind_d5, ind_x2, ind_d2);
    
end
end
end
end


for ind_x4 = 1:x
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    XX4(ind_x4, ind_x2 + x*(ind_d2-1)) = X4(ind_x4, ind_x2, ind_d2);
    
end
end
end


for ind_d4 = 1:d
for ind_d2 = 1:d
for ind_x2 = 1:x 
    
    DD4(ind_d4, ind_x2 + x*(ind_d2-1)) = D4(ind_d4, ind_x2, ind_d2);
    
end
end
end


r


%% ========================================================================

% L = zeros(o, o, o, x, d);
% 
% %O2O6X2D2
%     for ind_o4 = 1:o
%     for ind_o5 = 1:o
%     for ind_o6 = 1:o
%         for ind_d2 = 1:d
%             for ind_x2 = 1:x
%                 
%                 for ind_d3 = 1:d
%                 for ind_d4 = 1:d
%                 for ind_d5 = 1:d  
%                 for ind_x3 = 1:x
%                 for ind_x4 = 1:x
%                 for ind_x5 = 1:x
%                 for ind_x6 = 1:x    
%                             
%                     L(ind_o4, ind_o5, ind_o6, ind_x2, ind_d2) = L(ind_o4, ind_o5, ind_o6, ind_x2, ind_d2) + ...
%                         O(ind_o4, ind_x4) * O(ind_o5, ind_x5) * O(ind_o6, ind_x6) * A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4) * A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
%                             
%                 end
%                 end
%                 end
%                 end
%                 end
%                 end
%                 end
%                 
%             end
%         end
%     end
%     end
%     end
% 
% 
% %normalize
% for ind_x2 = 1:x
% for ind_d2 = 1:d
%      L(:, :, :, ind_x2, ind_d2) = L(:, :, :, ind_x2, ind_d2) / sum(sum(sum(L(:, :, :, ind_x2, ind_d2))));
% end
% end
% 
% H = zeros(o*o*o, x*d);
% for ind_o6 = 1:o
% for ind_o5 = 1:o    
% for ind_o4 = 1:o
% for ind_d2 = 1:d
% for ind_x2 = 1:x 
%     
%     H(ind_o4 + o*(ind_o5-1) + o^2*(ind_o6-1), ind_x2 + x*(ind_d2-1)) = L(ind_o4, ind_o5, ind_o6, ind_x2, ind_d2);
%     
% end
% end
% end
% end
% end
% 
% 
% 
% %% ========================================================================
% 
% % O6X2D2 = zeros(o, x, d);
% % 
% % for ind_o6 = 1:o
% %     for ind_d2 = 1:d
% %        for ind_x2 = 1:x 
% %         
% %         for ind_x3 = 1:x
% %         for ind_x4 = 1:x
% %         for ind_x5 = 1:x
% %         for ind_x6 = 1:x
% %         for ind_d3 = 1:d
% %         for ind_d4 = 1:d
% %         for ind_d5 = 1:d    
% %            
% %             O6X2D2(ind_o6, ind_x2, ind_d2) = O6X2D2(ind_o6, ind_x2, ind_d2) + ...
% %                 O(ind_o6, ind_x6) * A(ind_x6, ind_x5, ind_d5) * D(ind_d5, ind_x5, ind_d4) * A(ind_x5, ind_x4, ind_d4) * D(ind_d4, ind_x4, ind_d3) * A(ind_x4, ind_x3, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
% %         end
% %         end
% %         end
% %         end
% %         end
% %         end
% %         end
% %        end
% %     end
% % end
% % 
% % 
% % %normalize
% % for ind_x2 = 1:x
% %     for ind_d2 = 1:d
% %         O6X2D2(:, ind_x2, ind_d2) = O6X2D2(:, ind_x2, ind_d2)/sum(O6X2D2(:, ind_x2, ind_d2));
% %     end
% % end
% % 
% % 
% % 
% % 
% % Q = zeros(o, x*d);
% % for ind_o6 = 1:o
% % for ind_d2 = 1:d 
% % for ind_x2 = 1:x
% %     
% %     Q(ind_o6, ind_x2 + x*(ind_d2-1)) = O6X2D2(ind_o6, ind_x2, ind_d2);
% %     
% % end
% % end
% % end
% 
% %% ========================================================================
% 
% % O3X2D2 = zeros(o,x,d);
% % 
% % for ind_o3 = 1:o
% %     for ind_d2 = 1:d
% %         for ind_x2 = 1:x
% %                         
% %             for ind_x3 = 1:x
% %                 
% %                 O3X2D2(ind_o3, ind_x2, ind_d2) = O3X2D2(ind_o3, ind_x2, ind_d2) + ...
% %                     O(ind_o3, ind_x3) * A(ind_x3, ind_x2, ind_d2);
% %             end
% %             
% %         end        
% %     end
% % end
% % 
% % %normalize
% % for ind_x2 = 1:x
% %     for ind_d2 = 1:d        
% %         O3X2D2(:, ind_x2, ind_d2) = O3X2D2(:, ind_x2, ind_d2)/sum(O3X2D2(:, ind_x2, ind_d2));        
% %     end
% % end
% % 
% % K = zeros(o, x*d);
% % for ind_o3 = 1:o
% % for ind_d2 = 1:d
% % for ind_x2 = 1:x 
% %     
% %     K(ind_o3, ind_x2 + x*(ind_d2-1)) = O3X2D2(ind_o3, ind_x2, ind_d2);
% %     
% % end
% % end
% % end
% 
% 
% %% ========================================================================
% 
% % W = zeros(o, o, x, d);
% % 
% % %O2O3X2D2
% % for ind_o2 = 1:o
% %     for ind_o3 = 1:o
% %         for ind_d2 = 1:d
% %             for ind_x2 = 1:x
% %                 
% %                 for ind_x3 = 1:x
% %                             
% %                     W(ind_o2, ind_o3, ind_x2, ind_d2) = W(ind_o2, ind_o3, ind_x2, ind_d2) + ...
% %                         O(ind_o3, ind_x3) * O(ind_o2, ind_x2) * A(ind_x3, ind_x2, ind_d2);
% %                             
% %                 end
% %                 
% %             end
% %         end
% %     end
% % end
% % 
% % 
% % %normalize
% % for ind_x2 = 1:x
% % for ind_d2 = 1:d
% %      W(:, :, ind_x2, ind_d2) = W(:, :, ind_x2, ind_d2) / sum(sum(W(:, :, ind_x2, ind_d2)));
% % end
% % end
% % 
% % S = zeros(o*o, x*d);
% % for ind_o3 = 1:o
% % for ind_o2 = 1:o
% % for ind_d2 = 1:d
% % for ind_x2 = 1:x 
% %     
% %     S(ind_o2 + o*(ind_o3-1), ind_x2 + x*(ind_d2-1)) = W(ind_o2, ind_o3, ind_x2, ind_d2);
% %     
% % end
% % end
% % end
% % end
% 
% %% ========================================================================
% 
% 
% V = zeros(o, o, x, d);
% 
% %O5O6X2D2
% for ind_o5 = 1:o
%     for ind_o6 = 1:o
%         for ind_d2 = 1:d
%             for ind_x2 = 1:x
%                 
%                 for ind_d3 = 1:d
%                 for ind_d4 = 1:d
%                 for ind_d5 = 1:d  
%                 for ind_x3 = 1:x
%                 for ind_x4 = 1:x
%                 for ind_x5 = 1:x
%                 for ind_x6 = 1:x    
%                             
%                     V(ind_o5, ind_o6, ind_x2, ind_d2) = V(ind_o5, ind_o6, ind_x2, ind_d2) + ...
%                         O(ind_o6, ind_x6) * O(ind_o5, ind_x5) * A(ind_x6, ind_x5, ind_d5) * A(ind_x5, ind_x4, ind_d4) * D(ind_d5, ind_x5, ind_d4) * A(ind_x4, ind_x3, ind_d3) * D(ind_d4, ind_x4, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
%                             
%                 end
%                 end
%                 end
%                 end
%                 end
%                 end
%                 end
%                 
%             end
%         end
%     end
% end
% 
% 
% %normalize
% for ind_x2 = 1:x
% for ind_d2 = 1:d
%      V(:, :, ind_x2, ind_d2) = V(:, :, ind_x2, ind_d2) / sum(sum(V(:, :, ind_x2, ind_d2)));
% end
% end
% 
% G = zeros(o*o, x*d);
% for ind_o6 = 1:o
% for ind_o5 = 1:o
% for ind_d2 = 1:d
% for ind_x2 = 1:x 
%     
%     G(ind_o5 + o*(ind_o6-1), ind_x2 + x*(ind_d2-1)) = V(ind_o5, ind_o6, ind_x2, ind_d2);
%     
% end
% end
% end
% end












