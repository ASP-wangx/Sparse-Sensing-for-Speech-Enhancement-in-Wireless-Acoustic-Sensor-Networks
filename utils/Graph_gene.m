%% Graph generation
function A = Graph_gene(thr,mic)
%% setup
M = size(mic,1);% mic number

%% distance based communication graph
D0 = pdist2(mic,mic,'Euclidean');
A = zeros(M,M);
for i = 1:M
    for j = 1:M
       if(D0(i,j)<=thr) 
           A(i,j) = 1;
       end
    end
end
A = A - diag(diag(A));   % Adjacency matrix


end
