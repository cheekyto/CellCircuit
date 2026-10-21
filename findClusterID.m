function clusterID = findClusterID(spikeTimestamp, id, n_neighbor, option)

% option = 1;         % =1, find cluster position using the field "cluster"
                      % =2, find cluster position using the field "position"
                      % =3, find neighboring clusters' position using the field "cluster"

MATRIX_SIZE = 64;
N=2*n_neighbor+1;

n_clusters = size(spikeTimestamp,2);
data = nan(n_clusters,1);

switch option
    case 1
        for i=1:n_clusters
            data(i) = spikeTimestamp(i).cluster;
        end
    case 2
        for i=1:n_clusters
            data(i) = spikeTimestamp(i).position;
        end
    case 3
        for i=1:n_clusters
            data(i) = spikeTimestamp(i).cluster;
        end
end

clusterID = [];
for i=1:size(spikeTimestamp,2)
    if data(i) == id
        clusterID = vertcat(clusterID, i);
    end
end

if numel(clusterID)<1
    clusterID = nan;
end

if option ==3
    idx_center = clusterID;
    clusterID = nan(N*N,1);
    contador=0;
    for x=(-1*n_neighbor):1:n_neighbor
        for y=(-1*n_neighbor):1:n_neighbor
            contador = contador + 1;
            if x==0 && y==0
                clusterID(contador) = idx_center;
            else
                aux_idx = id + MATRIX_SIZE*y + x; 
                for i=1:size(spikeTimestamp,2)
                    if data(i) == aux_idx
                        clusterID(contador) = i;
                    end
                end
            end
        end
     end
end

end
               