import pandas as pd
from csv import reader
import numpy as np

def writeExELEMFile(data_file, nodes,ele_name):
    """
        Write out an ex Node file with the given coords.
        :param filename: Filename to write to.
        :param coords: List of coordinate lists.
        :return: None
    """
    with open( ele_name, 'w') as f:
        f.write(' Group name: Converted_tree\n')
        f.write(' Shape.  Dimension=1\n')
        f.write(' #Scale factor sets=1\n')
        f.write('  l.lagrange, #Scale factors= 2\n')
        f.write(' #Nodes= 2\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('    x.  l.Lagrange, no modify, standard node based.\n')
        f.write('      #Nodes= 2\n')
        f.write('       1.  #Values=1\n')
        f.write('        Value indices:   1\n')
        f.write('        Scale factor indices:   1\n')
        f.write('       2.  #Values=1\n')
        f.write('        Value indices:   1\n')
        f.write('        Scale factor indices:   2\n')
        f.write('    y.  l.Lagrange, no modify, standard node based.\n')
        f.write('      #Nodes= 2\n')
        f.write('      1.  #Values=1\n')
        f.write('        Value indices:   1\n')
        f.write('        Scale factor indices:   1\n')
        f.write('      2.  #Values=1\n')
        f.write('        Value indices:   1\n')
        f.write('        Scale factor indices:   2\n')
        f.write('    z.  l.Lagrange, no modify, standard node based.\n')
        f.write('      #Nodes= 2\n')
        f.write('      1.  #Values=1\n')
        f.write('        Value indices:   1\n')
        f.write('        Scale factor indices:   1\n')
        f.write('      2.  #Values=1\n')
        f.write('        Value indices:   1\n')
        f.write('        Scale factor indices:   2\n')

        with open(data_file, 'r') as read_obj:
            # pass the file object to reader() to get the reader object
            csv_reader = reader(read_obj)
            # Iterate over each row in the csv using reader object
            count=1
            for row in csv_reader:
                in_node=np.where(np.all(nodes[:,0:3]==row[0:3],axis=1))
                out_node=np.where(np.all(nodes[:,0:3]==row[3:6],axis=1))
                f.write(' Element:         %.d 0 0\n' % (count))
                f.write('   Nodes:         \n       %.d\t' % (nodes[in_node,3]))
                f.write('    %d\n' % (nodes[out_node,3]))
                f.write('   Scale factors:        \n     0.10000E+01  0.10000E+01\n')
                count=count+1
def filter_nodes(data_file):

    with open(data_file, 'r') as read_obj:
        # pass the file object to reader() to get the reader object
        csv_reader = reader(read_obj)
        # Iterate over each row in the csv using reader object

        nodes = [[0 for x in range(3)] for y in range(3140)]  # CIP data


        count = 0
        for row in csv_reader:
            # row variable is a list that represents a row in csv
            nodes[count] = row[0:3]
            count = count + 1
            nodes[count] = row[3:6]
            count = count + 1
    nodes = np.array(nodes)
    nodes_re, indices = np.unique(nodes, axis=0, return_index=True) # remove duplicate nodes. store value in nodes_re and store index in indices. axis=0 just show the unique rows
    indices = sorted(indices)# sort indices
    nodes = nodes[indices]# allocate nodes based on the sorted indices
    add_col= [i for i in range(1,len(nodes)+1)]
    add_col=np.array([add_col])
    add_col=add_col.T #transpose vector as a coloum
    nodes=np.append(nodes,add_col,axis=1)
    return nodes
def writeExNodeFile( nodes,node_file):
    """
        Write out an ex Node file with the given coords.
        :param filename: Filename to write to.
        :param coords: List of coordinate lists.
        :return: None
    """
    with open(node_file, 'w') as f:
        f.write(' Group name : converted_tree\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 0\n')
        f.write('   y.  Value index= 2, #Derivatives= 0\n')
        f.write('   z.  Value index= 3, #Derivatives= 0\n')
        for i in range(len(nodes)):
            f.write(' Node:         %d\n' % (i+1))
            f.write('    %s\n' % (nodes[i][0]))
            f.write('    %s\n' % (nodes[i][1]))
            f.write('    %s\n' % (nodes[i][2]))


def filter_radius(data_file):

    with open(data_file, 'r') as read_obj:
        # pass the file object to reader() to get the reader object
        csv_reader = reader(read_obj)

        radius = [[0 for x in range(1)] for y in range(421)] # Graph_edges

        count = 0
        for row in csv_reader:
            # row variable is a list that represents a row in csv
            radius[count] = row[6:7]  # radii
            # radius[count] = row[7:8] #radii
            # radius[count] = row[8:9] #diameter
            # radius[count] = row[9:10] # ratio diameter
            # radius[count] = row[10:11] # ratio len to diameter
            # radius[count] = row[6:7] # len

            count = count + 1

    radius = np.array(radius)
    # radius_re, indices = np.unique(radius, axis=1, return_index=True) # remove duplicate row. store value in nodes_re and store index in indices. axis=0 just show the unique rows
    # indices = sorted(indices)# sort indices
    # radius = radius[indices]# allocate nodes based on the sorted indices
    return radius


def writeExRadiusFile(radius, file_name):
    """
        Write out an ex Node file with the given coords.
        :param filename: Filename to write to.
        :param coords: List of coordinate lists.
        :return: None
    """
    with open(file_name, 'w') as f:
        f.write('Group name: converted_tree\n')
        f.write(' Shape.  Dimension=1\n')
        f.write(' #Scale factor sets=0\n')
        f.write(' #Nodes= 0\n')
        f.write(' #Fields= 1\n')
        f.write(' 1) radius, field, rectangular cartesian, #Components=1\n')
        f.write('  radius.  l.Lagrange, no modify, grid based.\n')
        f.write('  #xi1=1\n')
        # with open(data_file, 'r') as read_obj:
        #     csv_reader = reader(read_obj)
        # count = 1
        # for row in csv_reader:
        #     find_radius = np.where(np.all(radius[:, 7] == row[7:8], axis=1))
        #     f.write(' Element:         %.d 0 0\n' % (count))
        #     f.write(' Values:         \n        %d' % (radius[find_radius, 7]))
        #     f.write('    %d\n' % (radius[find_radius, 7]))
        for i in range(len(radius)):
            f.write(' Element:         %.d 0 0\n' % (i+1))
            f.write(' Values:       \n  %s\t' % (radius[i][0]))
            f.write('    %s\n' % (radius[i][0]))

def writeExRadiusGRAPHEDGESFile(radius, file_name):
    """
        Write out an ex Node file with the given coords.
        :param filename: Filename to write to.
        :param coords: List of coordinate lists.
        :return: None
    """
    with open(file_name, 'w') as f:
        f.write('CMISS Version 1.21 ipfiel File Version 3\n')
        f.write(' Heading: test\n')
        f.write(' \n')
        f.write(' The number of elements is [  935]: 20\n')
        f.write(' \n')
        # f.write(' Element number:        1\n')
        # f.write(' 1) radius, field, rectangular cartesian, #Components=1\n')
        # f.write('  radius.  l.Lagrange, no modify, grid based.\n')
        # f.write('  #xi1=1\n')
        # with open(data_file, 'r') as read_obj:
        #     csv_reader = reader(read_obj)
        # count = 1
        # for row in csv_reader:
        #     find_radius = np.where(np.all(radius[:, 7] == row[7:8], axis=1))
        #     f.write(' Element:         %.d 0 0\n' % (count))
        #     f.write(' Values:         \n        %d' % (radius[find_radius, 7]))
        #     f.write('    %d\n' % (radius[find_radius, 7]))
        for i in range(len(radius)):
            f.write(' Element number:         %.d \n' % (i+1))
            # f.write(' The field variable value is [ 0.00000D+00]:        %s\t' % (radius[i][0]))
            f.write(' The field variable value is [ 0.00000D+00]:   %s\n' % (radius[i][0]))
            f.write(' \n')




