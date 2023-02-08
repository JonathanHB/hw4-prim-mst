import numpy as np
import heapq as hq
from typing import Union

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):

        #arbitrary
        init_node = 0

        #number of nodes
        n = self.adj_mat.shape[0]
        #nodes_visited = np.zeros(n) #one hot encoding
        nodes_visited = [init_node]
        #unused_edges = [[x, init_node, w] for x, w in enumerate(self.adj_mat[init_node])]
        #unused_edges.sort(key=lambda x: x[2])
        h = []
        [hq.heappush(h, (w, init_node, x)) for x, w in enumerate(self.adj_mat[init_node]) if w != 0]

        self.mst = np.zeros(self.adj_mat.shape)

        while True:
            #get the lowest-weighted edge
            try:
                edge = hq.heappop(h)
            except IndexError as e:
                raise IndexError("graph is probably disconnected, which this implementation does not support.") from e

            #this check is necessary for cyclic graphs
            if edge[2] not in nodes_visited:
                #add the lowest-weighted edge's destination node's edges to the heap
                [hq.heappush(h, (w, edge[2], x)) for x, w in enumerate(self.adj_mat[edge[2]]) if w != 0]
                self.mst[edge[1]][edge[2]] = edge[0]
                self.mst[edge[2]][edge[1]] = edge[0]
                nodes_visited.append(edge[2])
                #print(nodes_visited)
                #print(edge)
                if len(nodes_visited) == n:
                    break
            #else:
            #    print(f"node {edge[2]} already in heap; skipping it")

        #print(self.adj_mat)
        #print(self.mst)


        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """

        #self.mst = None

# file_path = '../data/small.csv'
# g = Graph(file_path)
# g.construct_mst()