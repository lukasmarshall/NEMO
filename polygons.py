# Copyright (C) 2014, 2015, 2016 The University of New South Wales
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

"""Support code for the 43 polygons of the AEMO study."""

import json
import numpy as np
import dijkstra
from latlong import LatLong
import regions
import networkx as nx
import networkGraph
from geopy.distance import vincenty

# Useful for testing (in polygon 42 region model)
wildcard = 31

# Get the networkX-based graph model of the Transmission Network, overlayed with generator information. 
G = networkGraph.getNetworkAndGeneratorGraph()
# Get a sorted reference list of nodes in G.
node_keys = sorted(list(G.nodes()))
total_len = 0
for edge in G.edges():
    total_len += G.edge[edge[0]][edge[1]]["length"]
print "total Length"
print total_len

# These dictionary objects will hold the fraction of a region's load in each polygon.
regions.nsw.polygons = {}
regions.qld.polygons = {}
regions.sa.polygons = {}
regions.snowy.polygons = {}
regions.tas.polygons = {}
regions.vic.polygons = {}

# Add node key indices to the region polygons lists. Eventually going to make these node lists. Stick with polygons for the moment. 
# Also - what is the 'snowy' region? Left it blank for now. ACT? At the moment ACT is put in as part of nsw for simplicity. 
for idx, node in enumerate(node_keys):
    if G.node[node]['state'] == "New South Wales" or  G.node[node]['state'] == "Australian Capital Territory":
        regions.nsw.polygons[idx] = 0
    elif G.node[node]['state'] == "Queensland":
        regions.qld.polygons[idx] = 0
    elif G.node[node]['state'] == "South Australia":
        regions.sa.polygons[idx] = 0
    elif G.node[node]['state'] == "Victoria":
        regions.vic.polygons[idx] = 0
    elif G.node[node]['state'] == "Tasmania":
        regions.tas.polygons[idx] = 0
    else:
        print G.node[node]['state']

# Arbitrarily spread demand across every node in a region. This needs to change (badly) but I just want to get the code working. 
# TODO: Get populations in the network graph so we can have realistic demand. 
for region in regions.All:
    for key in list(region.polygons):
        region.polygons[key] = float(1) / float(len(list(region.polygons)))


# Ensure all weights sum to one.
for r in regions.All:
    if len(r.polygons) > 0:
        assert round(sum(r.polygons.values())) == 1


_polygons = {}

# Vertices of the closed polygons (nb. from Ben - must be closed)
# Here we instead add single coordinate tuples as the polygons. Note that the code is expecting a tuple of tuples, hence the ((double brackets)).  
for idx, node in enumerate(node_keys):
    _polygons[idx] = ((node),)

numpolygons = len(_polygons)

# Table mapping polygon number to region.
_region_table = [None] * (numpolygons + 1)
for rgn in [regions.nsw, regions.qld, regions.sa, regions.tas, regions.vic]:
    for poly in rgn.polygons:
        _region_table[poly] = rgn


def region(poly):
    """Return the region a polygon resides in.

    >>> region(1)
    QLD1
    >>> region(40)
    TAS1
    """
    return _region_table[poly]

# These map a limit to an index (polygon number is used as index.)
# Make the wind, PV, CST limits 1 MW for all regions. 
wind_limit = [None].append([1] * (len(node_keys) + 1))
pv_limit = [None].append([1] * (len(node_keys) + 1))
cst_limit = [None].append([1] * (len(node_keys) + 1))



def dumps():
    """Dump the polygon data in GeoJSON format.

    >>> assert len(dumps()) > 0
    """
    polys = []
    for i in range(1, numpolygons):
        geometry = {'type': 'Polygon', 'coordinates': [_polygons[i]]}
        properties = {'Polygon #': i, 'CST limit (GW)': cst_limit[i], 'PV limit (GW)': pv_limit[i], 'Wind limit (GW)': wind_limit[i]}
        feature = {'type': 'Feature', 'geometry': geometry, 'properties': properties}
        polys.append(feature)
    top = {'type': 'FeatureCollection', 'features': polys}
    return json.dumps(top)


def _centroid(vertices):
    """Find the centroid of a polygon."""

    # Ensure the polygon is closed
    if len(vertices) > 1:
        assert vertices[0] == vertices[-1]
        thesum = 0
        vsum = (0, 0)
        for i in range(len(vertices) - 1):
            v1 = vertices[i]
            v2 = vertices[i + 1]
            cross = v1[0] * v2[1] - v1[1] * v2[0]
            thesum += cross
            vsum = (((v1[0] + v2[0]) * cross) + vsum[0], ((v1[1] + v2[1]) * cross) + vsum[1])
            z = 1. / (3. * thesum)
        return (vsum[0] * z, vsum[1] * z)
    else:
        return vertices[0]

centroids = {}
for i, vertices in _polygons.iteritems():
    a, b = _centroid(vertices)
    centroids[i] = LatLong(b, a)


def path(poly1, poly2):
    """
    Return a path from polygon 1 to polygon 2.

    >>> path(1, 30)
    [(1, 4), (4, 10), (10, 16), (16, 23), (23, 30)]
    >>> path(23, 43)
    [(23, 30), (30, 35), (35, 38), (38, 41), (41, 43)]
    """
    return connections[(poly1, poly2)]


def subset(path, polysuperset):
    """
    Are all polygons in path present in superset?

    >>> subset([(1,2), (2,3)], [1,2,3])
    True
    >>> subset([(1,4), (4,3)], [1,2,3])
    False
    """
    # Flatten the list of pairs into one long list.
    polylist = [i for sub in path for i in sub]
    # Now for a simple set operation.
    return set(polylist) <= set(polysuperset)


def direct_p(poly1, poly2):
    """
    Return True if region A and B are directly connected.

    >>> direct_p(1, 2)
    True
    >>> direct_p(1, 40)
    False
    """
    return len(path(poly1, poly2)) <= 1


def dist(poly1, poly2):
    """Return the distance between two polygon centroids.

    >>> dist(1,1)
    0
    >>> dist(1,43)
    2910
    >>> dist(1,43) == distances[1,43]
    True
    """
    # print "Getting dist: ("+str(poly1)+" , "+str(poly2)+")"

    # Get the nx node objects that correspond to the two poly ids that were passed. 
    nx_node_1 = node_keys[poly1 - 1]
    nx_node_2 = node_keys[poly2 - 1]

    distance = vincenty(nx_node_1, nx_node_2).kilometers

    # print "Getting dist: ("+str(nx_node_1)+" , "+str(nx_node_2)+")"
    # print distance

    return distance


def pathlen(path):
    """Return the total length of a path

    >>> pathlen([])
    0
    >>> x = dist(1,4)
    >>> assert pathlen([(1, 4)]) == x
    >>> y = dist(4,7)
    >>> assert pathlen([(1, 4), (4, 7)]) == x + y
    """
    # Get the networkX node ids (tuples) from the path. 
    path_nx_nodes = [node_keys[n-1] for n in path]
    # iterate to find the distance. 
    dist = 0
    for i in range(len(path)-1):
        dist += dist(path[i], path[i+1])
    return dist


# A proposed transmission network.

# net = {1: {2: dist(1, 2), 3: dist(1, 3), 4: dist(1, 4)},
#        2: {1: dist(2, 1), 3: dist(2, 3), 5: dist(2, 5)},
#        3: {1: dist(3, 1), 2: dist(3, 2), 4: dist(3, 4), 6: dist(3, 6)},
#        4: {1: dist(4, 1), 3: dist(4, 3), 6: dist(4, 6), 7: dist(4, 7), 10: dist(4, 10)},
#        5: {2: dist(5, 2), 6: dist(5, 6), 8: dist(5, 8)},
#        6: {3: dist(6, 3), 4: dist(6, 4), 5: dist(6, 5), 9: dist(6, 9)},
#        7: {4: dist(7, 4), 10: dist(7, 10), 11: dist(7, 11), 16: dist(7, 16)},
#        8: {5: dist(8, 5), 9: dist(8, 9), 14: dist(8, 14)},
#        9: {6: dist(9, 6), 8: dist(9, 8), 10: dist(9, 10), 15: dist(9, 15)},
#        10: {4: dist(10, 4), 7: dist(10, 7), 9: dist(10, 9), 16: dist(10, 16)},
#        11: {7: dist(7, 11), 16: dist(16, 11), 17: dist(11, 17)},
#        12: {13: dist(12, 13), 19: dist(12, 19)},
#        13: {12: dist(13, 12), 14: dist(13, 14), 20: dist(13, 20)},
#        14: {8: dist(14, 8), 13: dist(14, 13), 15: dist(14, 15), 21: dist(14, 21)},
#        15: {9: dist(15, 9), 14: dist(15, 14), 16: dist(15, 16), 22: dist(15, 22)},
#        16: {7: dist(16, 7), 10: dist(16, 10), 11: dist(16, 11), 15: dist(16, 15), 17: dist(16, 17), 23: dist(16, 23), 24: dist(16, 24)},
#        17: {11: dist(17, 11), 16: dist(17, 16), 24: dist(17, 24)},
#        18: {19: dist(18, 19), 25: dist(18, 25)},
#        19: {12: dist(19, 12), 18: dist(19, 18), 20: dist(19, 20), 26: dist(19, 26)},
#        20: {13: dist(20, 13), 19: dist(20, 19), 21: dist(20, 21), 27: dist(20, 27)},
#        21: {14: dist(21, 14), 20: dist(21, 20), 22: dist(21, 22), 29: dist(21, 29)},
#        22: {15: dist(22, 15), 21: dist(22, 21), 23: dist(22, 23), 29: dist(22, 29)},
#        23: {16: dist(23, 16), 22: dist(23, 22), 24: dist(23, 24), 30: dist(23, 30)},
#        24: {16: dist(24, 16), 17: dist(24, 17), 23: dist(24, 23), 31: dist(23, 31)},
#        25: {18: dist(25, 18), 26: dist(25, 26)},
#        26: {19: dist(26, 19), 25: dist(26, 25), 27: dist(26, 27)},
#        27: {20: dist(27, 20), 26: dist(27, 26), 28: dist(27, 28), 32: dist(27, 32), 33: dist(27, 33)},
#        28: {21: dist(28, 21), 27: dist(28, 27), 29: dist(28, 29), 33: dist(28, 33)},
#        29: {22: dist(29, 22), 28: dist(29, 28), 30: dist(29, 30), 34: dist(29, 34)},
#        30: {23: dist(30, 23), 29: dist(30, 29), 31: dist(30, 31), 35: dist(30, 35)},
#        31: {24: dist(31, 24), 30: dist(31, 30), 36: dist(31, 36)},
#        32: {27: dist(32, 27), 33: dist(32, 33), 37: dist(32, 37)},
#        33: {27: dist(33, 27), 28: dist(33, 28), 32: dist(33, 32), 34: dist(33, 34), 37: dist(33, 37), 39: dist(33, 39)},
#        34: {29: dist(34, 29), 33: dist(34, 33), 35: dist(34, 35), 39: dist(34, 39)},
#        35: {30: dist(35, 30), 34: dist(35, 34), 36: dist(35, 36), 38: dist(35, 38), 39: dist(35, 39)},
#        36: {31: dist(36, 31), 35: dist(36, 35), 38: dist(36, 38)},
#        37: {32: dist(37, 32), 33: dist(37, 33), 39: dist(37, 39)},
#        38: {35: dist(38, 35), 36: dist(38, 36), 39: dist(38, 39), 41: dist(38, 41)},
#        39: {33: dist(39, 33), 34: dist(39, 34), 35: dist(39, 35), 37: dist(39, 37), 38: dist(39, 38), 40: dist(39, 40)},
#        40: {39: dist(40, 39), 41: dist(40, 41), 42: dist(40, 42)},
#        41: {38: dist(41, 38), 40: dist(41, 40), 43: dist(41, 43)},
#        42: {40: dist(42, 40), 43: dist(42, 43)},
#        43: {41: dist(43, 41), 42: dist(43, 42)}}

distances = networkGraph.getFromPickle('./pickles/distances.pkl')
if distances is None:
    distances = np.zeros((numpolygons + 1, numpolygons + 1))
    # mark row 0 and column 0 as unused (there is no polygon #0)
    print "Creating distances matrix."
    distances[0] = np.nan
    distances[::, 0] = np.nan
    for p1 in range(1, distances.shape[0]):
        for p2 in range(1, distances.shape[0]):
            distances[p1, p2] = dist(p1, p2)
    networkGraph.saveToPickle(distances, './pickles/distances.pkl')

print "Creating Limit matrix."
existing_net = np.zeros((numpolygons + 1, numpolygons + 1))
# mark row 0 and column 0 as unused (there is no polygon #0)
# Existing_net appears to store transfer limits between network nodes. 
# Leaving it as-is for the moment - not sure what to do with this until we make it into a timeseries thing based on transmission data. 
existing_net[0] = np.nan
existing_net[::, 0] = np.nan

for (p1, p2, limit) in  []:
    # Ensure connection between these nodes actually exists in the graph. 
    # G[] is the adjacency dictionary 

    assert node_keys[p1-1] in list(G[node_keys[p2-1]])
    existing_net[p1, p2] = limit


print "Getting Connections "

connections = networkGraph.getFromPickle('./pickles/connections.pkl')
if connections is None:
    print "Pickle not found. Generating connections. This can take ages."
    connections = {}
    for dest in range(1, numpolygons + 1):
        for src in range(1, numpolygons + 1):
            src_node = node_keys[src-1]
            dest_node = node_keys[dest-1]
            shortest_path_nodes = nx.shortest_path(G, source=src_node, target=dest_node, weight="length")
            shortest = [node_keys.index(n) for n in shortest_path_nodes]
            pairs = []
            for i in range(len(shortest) - 1):
                pairs.append((shortest[i], shortest[i + 1]))
            connections[(src, dest)] = pairs
    networkGraph.saveToPickle(connections, './pickles/connections.pkl')




print "Finished setting up electricity network model. "
# If run as a script in itself, print some info about the gens. 
if __name__ == '__main__':
    print "Testing Polygons."