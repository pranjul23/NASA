%using parametric distribution,
%calculate duration in state
function d = duration(i, Nstates)

persistent shape scale;

%initialize
if isempty(shape)    
    %Gamma distribution
    shape = 10*rand(Nstates,1);
    scale = 5*rand(Nstates,1);        
end

d = gamrnd(shape(i), scale(i));