from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    #check that the distance between A and B is the same as between B and A
    assert cluster.compute_similarity(activesite_a, activesite_b) == cluster.compute_similarity(activesite_b, activesite_a)
    #check that the distance between A and A is 0
    assert cluster.compute_similarity(activesite_a, activesite_a) == 0.0
    #check that the distance is always positive
    assert cluster.compute_similarity(activesite_a, activesite_a) >= 0.0


def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # clusters produce k number of final labels
    assert len(cluster.cluster_by_partitioning(active_sites)) == k

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # clusters produce k number of final labels
    assert len(cluster.cluster_hierarchically(active_sites)) == k
