import numpy as np
import pandas as pd
from kneebow.rotor import Rotor
from sklearn.cluster import AgglomerativeClustering

def computeClustering(dataNorm, matNeighbors, nbClust=None):
    """
    Function to compute 2D agglomerative clustering with neighboring constraint
    
    INPUTS:
    ----------
    dataNorm: input normalized data
    matNeighbors: connectivity matrix 
    nbClust: number of clusters to obtain, optional
    
    OUTPUTS:
    ----------
    merge: output of AgglomerativeClustering, children of each non-leaf node
    height: output of AgglomerativeClustering, distance between nodes at each step
    nbClust: number of clusters selected for clustering
    clustering: clustering of dataNorm elements
    """
    
    ##Select only relevant columns
    dataNorm = dataNorm.iloc[:,2:]
    # Data is log-transformed
    data = np.log(dataNorm + 1)

    if nbClust is None:
      ward = AgglomerativeClustering(
          n_clusters=1, connectivity=matNeighbors, linkage="ward", 
          compute_distances=True
      ).fit(data)
    else:
      ward = AgglomerativeClustering(
          n_clusters=nbClust, connectivity=matNeighbors, linkage="ward", 
          compute_distances=True, compute_full_tree=True
      ).fit(data)
    
    merge = pd.DataFrame(ward.children_)
    height = pd.DataFrame(ward.distances_)
    
    if nbClust is None:
      nbClust = getElbow(height)
      ward = AgglomerativeClustering(
          n_clusters=nbClust, connectivity=matNeighbors, linkage="ward", 
          compute_distances=True, compute_full_tree=True
      ).fit(data)

    clustering = ward.labels_ + 1
    
    return merge, height, nbClust, clustering

def getElbow(height):
    """
    Obtain hierarchy level corresponding to elbow of increase of total within-cluster variance (not monotonous)
    
    INPUTS:
    ----------
    height: output of AgglomerativeClustering, distance between nodes at each step
    
    OUTPUTS:
    ----------
    nbOptClust: optimal number of clusters, corresponding to height elbow
    """
    crit = np.column_stack((range(1, len(height) + 1), height))
    rotor = Rotor()
    rotor.fit_rotate(crit)
    elbowIndex = rotor.get_elbow_index()
    
    nbInt = len(height)+1
    nbOptClust = nbInt - (elbowIndex + 1)
    return nbOptClust
  
