function [K_nearest_neighbors, Distances] = Dijkstra_with_stopping2_LLPD(IDX,D,s)
%        [Distances, K_nearest_neighbors] = Dijkstra_with_stopping2(IDX,D,s,k)
% this function implements a dijkstra-like algorithm for finding the k
% nearest neighbors of a source vertex s, with respect to the longest leg
% path distance. Designed to run in a k-nearest neighbors graph
%
% INPUT
% ======================================================================= %
% IDX ............... n-by-k matrix. j-th row contains k nearest neighbors
% of vertex j.
% D ................. n-by-k matrix. Contains the distances corresponding
% to IDX.
% s ................. source vertex
%
% OUTPUT
% ======================================================================= %
% K_nearest_neighbors ........ k nearest neighbors in path distance
% Distances .................. Corresponding path lengths.
%
% Daniel Mckenzie
% December 2018

[n,k] = size(IDX);

% ====== find the neighbors of s, initialize min-priority queue ====== %
% NB: Q is set up to be a max priority queue, instead of a min priority
% queue as we really need. To work around this all distances are keyed into
% the queue as their negatives. So that the closest node has the highest
% priority.

Neighbors_of_s = IDX(s,:);
Dists_to_neighbors = D(s,:);
Q = pq_create(n);
List_of_Dist = -inf(n,1);  %list keeping track of key-index pairs in queue

% ======== Populate Q and List_of_Dist with neighbors of s ========== %
% Note that s is stored, by default, as the nearest neighbor of itself at a
% distance of zero
List_of_Dist(Neighbors_of_s) = -Dists_to_neighbors;
for i = 1:k
    pq_push(Q,Neighbors_of_s(i),-Dists_to_neighbors(i));
end

% ======================= Initialize output arrays  ==================== %
K_nearest_neighbors = zeros(k,1);
Distances = zeros(k,1);


% ============ Now perform k iterations of Dijkstra ============== % 
 for i = 1:k
     [curridx,currdist] = pq_pop(Q);
     K_nearest_neighbors(i) = curridx;
     Distances(i) = -currdist;
     Neighbors_to_examine = IDX(curridx,:);
     Dists_to_neighbors = D(curridx,:);
     for j = 2:k 
         % Examine all k-nearest neighbors of current vertex (curridx). 
         % If the length of path from s via curridx to j is shorter than
         % currently saved path, update it.
         tempdist = min(currdist,-Dists_to_neighbors(j));
         if tempdist > List_of_Dist(Neighbors_to_examine(j))
             pq_push(Q,Neighbors_to_examine(j),tempdist);
             List_of_Dist(Neighbors_to_examine(j)) = tempdist;
         end
     end
 end
 
 pq_delete(Q);
end
    
