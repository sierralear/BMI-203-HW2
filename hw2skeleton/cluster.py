from .utils import Atom, Residue, ActiveSite
import numpy as np
import pandas as pd
from collections import Counter
from itertools import product
from hw2skeleton import io
import os

my_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
         'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

path = os.path.join("data")
active_site = io.read_active_sites(path);

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0
    
    #creating list comprehension of all AAs for site a residues
    list_res_a = [r.type for r in site_a.residues]
    
    #creating list comprehension of all AAs for site b residues
    list_res_b = [r.type for r in site_b.residues]
    
    #create two histogram count lists for site a and b from the list comprehensions
    count_a = Counter(list_res_a)
    count_b = Counter(list_res_b)
    
    #initialize two 20-element AA dictionaries: one for site A, one for site B
    a_dict = {aa:0 for aa in my_aa}
    b_dict = {aa:0 for aa in my_aa}
    
    #convert normalized histogram count into the 20-D vector AA dictionary
    #a site
    for aa, count in count_a.items():
        a_dict[aa] = count
    a_vector = np.array(list(a_dict.values()))
    a_vector = a_vector/np.sum(a_vector) #to return percentage of amino acid identity
    #b site
    for aa, count in count_b.items():
        b_dict[aa] = count
    b_vector = np.array(list(b_dict.values()))
    b_vector = b_vector/np.sum(b_vector)
        
    #calculate Euclidian distance between the two sites' 20-D vectors
    similarity = np.sqrt(np.sum((a_vector - b_vector)**2))
    
    return similarity

def compute_similarity_partitioning(site_a, centroid):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: one ActiveSite instance, one "centroid" consisting of a 20-element vector with numnbers between 0 and 1
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0
    
    #creating list comprehension of all AAs for site a residues
    list_res_a = [r.type for r in site_a.residues]
    
    #create two histogram count lists for site a and b from the list comprehensions
    count_a = Counter(list_res_a)
    
    #initialize two 20-element AA dictionaries: one for site A, one for site B
    a_dict = {aa:0 for aa in my_aa}
    
    #convert normalized histogram count into the 20-D vector AA dictionary
    #a site
    for aa, count in count_a.items():
        a_dict[aa] = count
    a_vector = np.array(list(a_dict.values()))
    a_vector = a_vector/np.sum(a_vector) #to return percentage of amino acid identity
        
    #calculate Euclidian distance between the two sites' 20-D vectors
    similarity = np.sqrt(np.sum((a_vector - centroid)**2))
    
    return similarity, a_vector

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
        #number of clusters
    k = 3 #based off histogram--while an elbow plot would be a better way to pick clusters, I chose
              #this for the sake of brevity
    
    #creating random centroids
        # my centroid will be a 20-element vector with randomized numbers between 0 and 1
    np.random.seed(3) #seed
    centroid_list = [None] * k #initializing my list of k centroids
    for i in range(len(centroid_list)): 
        centroid_list[i] = np.random.rand(20) #creatign random numbers to fill 20 elements of each centroid vector
        centroid_list[i] /= centroid_list[i].sum() #to ensure all elements sum to 1 to mimic amino acid compositionality

    
    #k-means algorithm    
    partition_combos = list(product(range(0,len(active_site)), range(0,len(centroid_list)))) #create list of tuples
                                                                                          #with all combinations of
                                                                                #active site instances and centroids
    for iter in range(100):    #running KNN 100 times before stopping
        df = pd.DataFrame(columns = ['ActiveSite', 'Centroids', 'Similarity']) #initialize dataframe
        a_list = [] #initialize list of the compositionality of all amino acids
        for a, b in partition_combos:
            similarity, a_vector = compute_similarity_partitioning(active_site[a], centroid_list[b]) #calculate distances between all centroids and all active sites
            df = df.append({'ActiveSite': a, 'Centroids': b, 'Similarity': similarity}, ignore_index=True) #create df listing distance to each centroid for each amino acid
            a_list.append(a_vector) #create list of compositionality of all amino acids
    
        #labeling each ActiveSite to a centroid based on having the lowest similarity score
        df_new = df.sort_values("Similarity").groupby("ActiveSite").first() #assigns active site to closest centroid
    
        #recalculate centroid based on active site assigned to it
        clusterings = [] #initialize clusterings list
        for centroid in df_new["Centroids"].unique(): #looping through each centroid 
            as_in_centroids = df_new.loc[df_new["Centroids"] == centroid].index #find all active site labelled for that cluster
            clusterings.append(as_in_centroids.tolist()) #add that to my clustering list
            centroid_activesite_comps = [] #intialize list of compositionality for a given centroid
            for activesite in as_in_centroids: #lopping through each active site for a given centroid
                centroid_activesite_comps.append(a_list[int(activesite)]) #adding given active site's 20-element vector to list
            centroid_list[i] = np.mean(centroid_activesite_comps, axis = 0) #recalculating centroid as the row element-wise mean of the new amino acids
                
    return clusterings

as_num = list(range(0,len(active_site)))

def min_similarity(distance_matrix, group_labels):
    '''
    This function determines the the two active sites that are closest to each other and then changes it so those
    points and their associated cluster are converted into one single cluster. It also flags that minimum distance
    so it is never used again in future iterations.
    
    Input: distance_matrix -- a matrix returning all the distances from one active site to another active site
            group_labels -- the current clusters/labels all active sites currently belong to
    
    
    Output: row_i, col_i: the two active sites that are closest to each other
            distance: the minimum distance between those two active sites
    
    '''
    
    flat_i = np.nanargmin(distance_matrix) #finds the "flattened index" of the minimum distance
    row_i = int(flat_i/len(as_num)) #calculates what the row_index or one of the active sites from "flattened index"
    col_i = flat_i - (row_i * len(as_num)) #calculates the col_index or the other active site from "flattened index"
    distance = distance_matrix[row_i, col_i] #is the distance between those two active sites
        
    min_label = np.min([group_labels[row_i], group_labels[col_i]]) #calculates what the lower of the labels of those two active sites are
    max_label = np.max([group_labels[row_i], group_labels[col_i]]) #calculates what the higher of the labels of those two active sites are
    group_labels[np.where(group_labels == max_label)] = min_label #changes all labels of the maxmimum label to the minimum label
                                                                    #basically: combines two clusters into one cluster label
    
    distance_matrix[row_i, col_i] = group_labels[row_i] + 1000 #flagging that distance as "used up"
    
    return row_i, col_i, distance

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    group = np.array(as_num) #creating cluster labels
    dist_matrix = np.zeros((len(as_num), len(as_num))) # create upper triangular distance matrix recording similarities for each pair of active sites
    for row in as_num: #going through rows of distance matrix
        for col in range(row,len(as_num)): #going through columns of distance matrix
            similarity = compute_similarity(active_site[row], active_site[col]) #calculate similarity
            dist_matrix[row, col] = similarity #input the similarity into the proper location in matrix
        
    for row in as_num: #NaN-ing redundant values in matrix to make upper triangular matrix
        for col in range(0,row+1):
            dist_matrix[row, col] = np.NaN
            
    k = 3 #final number of clusters--set at 3 to better compare to the 3 clusters from Objective 2
    while len(np.unique(group)) > k: #the algorithm will run until I've reached my desired number of clusters
        min_similarity(dist_matrix, group) #calculating
        
    cluster_labels = np.unique(group) #returns a list of arrays of my clusters
    clusterings = []
    for cluster in cluster_labels:
        clusterings.append(np.where(group == cluster))
        ##Note: doing this somehow gave me a list returning tuples that then returns an array
        
        
    #instead turns the previous weird list of tuples of arrays into a list of list of lists (couldn't figure out how to mkae it just a list of lists)
    NN_clusters_list = [] 
    for cluster in clusterings:
        clustering = [l.tolist() for l in cluster]
        NN_clusters_list.append(clustering)
    


    return NN_clusters_list
