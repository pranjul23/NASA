function comboMatrix = createAllCombos(numValArray)

    totalVals = prod(numValArray);
    comboMatrix = ones(totalVals, length(numValArray));
    for n=2:1:totalVals
       comboMatrix(n,:) = comboMatrix(n-1,:); 
       index = length(numValArray);
       while(1)
        comboMatrix(n,index) = comboMatrix(n,index) + 1;
        if (comboMatrix(n,index) > numValArray(index))
           comboMatrix(n,index) = 1; 
           index = index - 1;
        else
            break
        end           
       end
    end
end
