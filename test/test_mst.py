import pytest
import numpy as np
from mst import Graph
from sklearn.metrics import pairwise_distances
import networkx as nx


def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001):
    """
    
    Helper function to check the correctness of the adjacency matrix encoding an MST.
    Note that because the MST of a graph is not guaranteed to be unique, we cannot 
    simply check for equality against a known MST of a graph. 

    Arguments:
        adj_mat: adjacency matrix of full graph
        mst: adjacency matrix of proposed minimum spanning tree
        expected_weight: weight of the minimum spanning tree of the full graph
        allowed_error: allowed difference between proposed MST weight and `expected_weight`

    TODO: Add additional assertions to ensure the correctness of your MST implementation. For
    example, how many edges should a minimum spanning tree have? Are minimum spanning trees
    always connected? What else can you think of?

    """

    mstg = nx.from_numpy_matrix(adj_mat, parallel_edges=False, create_using=None)
    assert nx.is_connected(mstg), "MST is not connected"


    def approx_equal(a, b):
        return abs(a - b) < allowed_error

    n_edges = 0
    total = 0
    for i in range(mst.shape[0]):
        for j in range(i+1):
            if mst[i,j] != 0:
                n_edges+=1
            total += mst[i, j]

    print(total)
    print(expected_weight)
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'
    assert n_edges == mst.shape[0]-1, f"MST has {n_edges} edges but should have {mst.shape[0]-1}."


def test_mst_small():
    """
    
    Unit test for the construction of a minimum spanning tree on a small graph.
    
    """
    file_path = '../data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 8)


def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = '../data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)


def test_mst_student():

    n_nodes = 4
    thresh = .2

    #generate an adjacency matrix for a connected graph with n_nodes nodes, each connected to a fraction 1-thresh of their neighbors
    adjout = []

    while True:
        print("iter")

        testadj = np.random.rand(n_nodes,n_nodes)
        for x in range(n_nodes):
            for y in range(0,x):
                #make a fraction thresh of the edges essentially nonexistent
                if testadj[x][y] < thresh:
                    testadj[x,y] = 10**100
                    testadj[y,x] = 10**100
                #rescale the remaining edges to [0,1]
                else:
                    rescaled = (testadj[x][y]-thresh)/(1-thresh)
                    testadj[x,y] = rescaled
                    testadj[y,x] = rescaled


        mstg = nx.from_numpy_matrix(testadj, parallel_edges=False, create_using=None)
        if nx.is_connected(mstg):
            adjout = testadj
            break

        print("regenerating disconnected graph")


    #test
    g = Graph(adjout)
    g.construct_mst()

    #networkx benchmark
    g2 = nx.from_numpy_matrix(adjout, parallel_edges=False, create_using=None)
    expected_mst_sum = sum([e[2]['weight'] for e in nx.minimum_spanning_edges(g2)])

    check_mst(g.adj_mat, g.mst, expected_mst_sum)

    """
    
    TODO: Write at least one unit test for MST construction.
    
    """
    pass
