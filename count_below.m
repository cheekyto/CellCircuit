function n_below = count_below(values, queries, inclusive)

% function to count, for each query, how many elements of "values" are below it (binary-search equivalent, for sorted or unsorted values)
% It replaces loops like "for i=1:numel(queries), n(i) = sum(values < queries(i)); end", which take n_values*n_queries operations
% (revised October 2026: it uses "histcounts", ~4x faster than sorting values and queries together)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% values:               array of values (e.g., spike times of a channel)
% queries:              array of queries (e.g., spike times of another channel)
% inclusive:            =1, count the values <= query;  =0, count the values < query (strict)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_below(i):           number of elements of "values" below queries(i) (column array)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

values = values(:);     queries = queries(:);
n_values = numel(values);   n_queries = numel(queries);
n_below = zeros(n_queries,1);
if n_queries == 0
    return
end

if inclusive                                                                % values <= query  <=>  NOT (-values < -query): use the strict count on the negated arrays
    n_below = sum(~isnan(values)) - count_below(-values, -queries, 0);      % (NaN values are never counted, as in the strict count)
    return
end

if ~any(isnan(queries))                                                     % (revised October 2026) usual case: histogram of the values, with the queries as bin edges
    [edges, ~, loc] = unique(queries);                                      % sorted, strictly increasing edges; loc(i): edge of the query i
    counts = histcounts(values, [-Inf; edges; Inf]);                        % counts(1): values < edges(1); counts(j): values in [edges(j-1), edges(j)) (NaN values are ignored)
    below = cumsum(counts(:));                                              % below(j): number of values < edges(j)
    n_below = below(loc);
    n_below = n_below(:);
    return
end

% Queries with NaN: sort values and queries together (slower); the row index breaks the ties, so that equal values go after the query (strict)
all_data = [queries, (1:n_queries)'; values, n_queries+(1:n_values)'];      % queries first: in a tie, the value is NOT counted
[~, order] = sortrows(all_data, [1 2]);                                     % order(k): row of "all_data" at the position k of the sorted list (sorted by value, then by row)
is_value = order > n_queries;
cum_values = cumsum(is_value);                                              % number of values up to each position of the sorted list
is_query = ~is_value;
n_below(order(is_query)) = cum_values(is_query);                            % for each query, the number of values sorted before it

end
