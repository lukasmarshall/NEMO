import networkx as nx
import matplotlib.pyplot as plt
import geojson
import geopy
from geopy.distance import vincenty

import pickle
import os



def saveToPickle(my_object, fileName):
    print("Pickling my_object to file: "+str(fileName)+"...")
    pickle.dump(my_object, open(fileName, "wb"))
    print ("Saved.")

def getFromPickle(fileName):
    if os.path.isfile(fileName):
        my_object = pickle.load(open(fileName, "rb"))
        return my_object
    else:
        return None


# G = nx.Graph()

# # Add a single node
# G.add_node(1)

# # Add a list of nodes
# G.add_nodes_from([2,3])


# # Add an edge between 1 and 2 and associate it with an object
# G.add_edge(1,2, object={'number':2})
# # You can do this for nodes and whole graphs also.
# # Weight is a special attribute and is used by graph algos. 

# # add an 'nbunch' of nodes (iterable container of nodes)
# H = nx.path_graph(10)
# G.add_nodes_from(H)

# # Add an edge from 1 to 2
# G.add_edge(1,2)

# # Add an edge tuple
# e = (2,3)
# G.add_edge(*e) # the * unpacks the edge tuple. cool. 

# # Add edges from a list
# G.add_edges_from([(1,2), (1,3)])

# # Add eges from an 'ebunch' iterable container of edges
# G.add_edges_from(H.edges())

# # Remove some nodes
# G.remove_nodes_from([2,3])

# # Print the nodes:
# print G.nodes()

# # Print the edges
# print G.edges()


unconnected_features = ['Heron Creek Tee', 'Middleback Tee', 'Steeple Flat', 'Bocco Rock','Nevertire', 'Dubbo Tee']


def constructGraph(features):

	G = nx.Graph()
	seen = []
	for feature in features:
		name = feature['properties']['NAME']
		state = feature['properties']['STATE']
		capacityKV = feature['properties']['CAPACITYKV']
		nodes = name.split(' to ')
		origin = nodes[0].strip()
		destination = nodes[len(nodes) - 1].strip()

		if not origin.isspace() and not origin == "" and not destination.isspace() and not destination == "" and not state == "Western Australia" and not state == "Northern Territory" and not origin in unconnected_features and not destination in unconnected_features:
			print name
			# The maximum distance 2 points can be apart to be considered the same node (in meters) 
			COLOCATION_DISTANCE = 150

			origin += " "+feature['properties']['STATE']
			destination += " "+feature['properties']['STATE']
			path = feature['geometry']['coordinates'][0] # LIST OF LAT/LONGS FROM START TO FINISH
			# We're rounding to 3 decimal places here as it gives us the approximate width of a substation as the error  
			# so roughly geolocated points will be treated as one node, which is what we want.
			origin_coords = {'lat':float(path[0][0]), 'lon':float(path[0][1])}
			dest_coords = {'lat':float(path[len(path) - 1][0]), 'lon':float(path[len(path) - 1][1])}
			# Edge information - shape length and labels, some other things.
			edge_info = {
				'length':float(feature['properties']['SHAPE_Length']), 
				'object_id':feature['properties']['OBJECTID'], 
				'path_name':name, 
				'node_names':[origin, destination],
				'capacityKV':capacityKV,
			}
			
			# Create node name tuples
			origin_node = (origin_coords['lat'], origin_coords['lon'])
			dest_node = (dest_coords['lat'], dest_coords['lon'])

			# Check distance against all other nodes to find effectively colocated spots.
			# This slows things down A LOT because we make the complexity n^2
			for node in list(G.nodes()):
				origin_distance =  vincenty((origin_coords['lat'], origin_coords['lon']), node).meters
				dest_distance = vincenty((dest_coords['lat'], dest_coords['lon']), node).meters
				# If the origin or the destination are close enough to an existing 
				if origin_distance < COLOCATION_DISTANCE:
					origin_node = node
				if dest_distance < COLOCATION_DISTANCE:
					dest_node = node
			
			# We create two nodes.
			G.add_node(
				origin_node, 
				pos=(origin_node[1], origin_node[0])
				)
			G.add_node(
				dest_node, 
				pos=(dest_node[1], dest_node[0])
				)
			# We create an edge using the two nodes
			G.add_edge(origin_node, dest_node, attr_dict=edge_info)
	
	return G









def getNetworkGraph():
	json_file =  open("networkMap/features.geojson")
	json_data = geojson.load(json_file)

	G = getFromPickle('./pickles/networkGraph.pkl')
	if not G:
		G = constructGraph(json_data.features)
		saveToPickle(G, './pickles/networkGraph.pkl')
	return G


def getNetworkAndGeneratorGraph():
	G = getNetworkGraph()
	G = getFromPickle('./pickles/networkAndGeneratorGraph.pkl')
	if not G:
		G = _addGenerators(G)
		saveToPickle(G, './pickles/networkAndGeneratorGraph.pkl')
	
	return G


def _addGenerators(G):
	print "Adding generators."
	# Get the geojson containing all gens in australia (2017)
	json_file =  open("generatorMap/generators.geojson")
	json_data = geojson.load(json_file)
	features = json_data.features
	# remove features in non-NEM states
	non_nem_states = ["Western Australia", "Northern Territory"]
	features[:] = [feature for feature in features if feature.get('properties').get('STATE') not in non_nem_states]

	# For every feature, find a corresponding nod ein the network graph to associate it with. 
	for feature in features:
		state = feature.properties['STATE']
		capacity = float(feature.properties['GENERATIONMW']) if feature.properties['GENERATIONMW'] else 0

		# Associate each feature with a graph node.
		# Check distance against all nodes to find probably corresponding node..
		# This is slow because we make the complexity n^2.
		nodes = list(G.nodes())
		min_distance = vincenty((feature.geometry.coordinates[0],feature.geometry.coordinates[1]), nodes[0]).meters
		closest_node = nodes[0]
		for node in list(G.nodes()):
			# Get the vincenty distance between a given node and the generator
			distance =  vincenty((feature.geometry.coordinates[0],feature.geometry.coordinates[1]), node).meters
			closest_node = node if distance < min_distance else closest_node
			min_distance = distance if distance < min_distance else min_distance
		feature.graph_info = { 'node':closest_node, 'node_distance': min_distance} #keep some info about which part of the graph the feature was assigned to, for future analysis


		# If the node is not currently a 'generators'type node, make it one.  
		if G.node[closest_node]['node_type'] != 'generator':
			G.node[closest_node]['node_type'] = 'generator'
			G.node[closest_node]['generators'] = []
			G.node[closest_node]['generation_capacity'] = 0
			G.node[closest_node]['generation_type'] = ""
			G.node[closest_node]['generation_label'] = ""

		# Add the current generator feature to the generators node
		G.node[closest_node]['generators'].append(feature)
		G.node[closest_node]['num_generators'] = len(G.node[closest_node]['generators'])
		
		if feature['properties']['GENERATIONMW']:
			G.node[closest_node]['generation_capacity'] += float(feature['properties']['GENERATIONMW'])

		# Give the node a 'fuel type'for graphing purposes. 
		# At the moment I'm just choosing whichever was the last primary fuel type seen at the node. Could do better but I want a single type label for each one for coloring purposes.
		if feature['properties']['PRIMARYFUELTYPE']:
			fuelType = feature['properties']['PRIMARYFUELTYPE']
			G.node[closest_node]['generation_type'] = fuelType
			if not fuelType in seen_gen_types:
				seen_gen_types.append(fuelType)
			G.node[closest_node]['generation_type_color_code'] = seen_gen_types.index(fuelType)

		if feature['properties']['NAME']:
			G.node[closest_node]['generation_label'] += feature['properties']['NAME']+" "


	# print features
	for feature in features:
		state = feature.properties['STATE']
		print state
		print feature
		print feature['graph_info']['node_distance']
	return G


def examineConnectedSubgraphs(G):
	counter = 0
	unconnected = []
	for comp in nx.connected_components(G):
		print str(len(comp)) + " " + str(comp)+"\n"
		counter += 1
		if len(comp) < 100:
			for node in comp:
				unconnected.append(node)

	# Loop through all edges, if in an unconnected set then print.
	for n1, n2, names in G.edges(data='node_names'):
		if n1 in unconnected or n2 in unconnected:
			print names
			print n1
			print n2

	print "Number of Graphs: "+str(counter)
	# Draw the graph
	