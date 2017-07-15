#Eric Bochsler and Daniel Stromecki
#Homework 7: Final Project
#Globular Cluster Characterization

import numpy as np
import matplotlib.pylab as plt
import urllib2
from scipy.stats import pearsonr

       
def print_file_lines(filename):
    '''Opens a file of the specified filename and prints each line of the file
    on a numbered line.
    '''
    my_file = urllib2.urlopen(filename, 'r')
    n = 0
    for line in my_file:
        n += 1
        print str(n), repr(line.strip()) 
    my_file.close()
    
    
def read_cluster_dat_file(url):
    '''Reads data in from a specified url pathway. The data is stored as a
    list of strings where each string corresponds to a line of data. Only data
    from the New General Catalogue is stored. The list of strings is returned.
    '''
    my_file = urllib2.urlopen(url, 'r')
    cluster_list = []
    label_list = []
    for line in my_file:
        line = line.strip()
        if test_if_ngc(line):
            cluster_list.append(line)
        elif test_if_label(line):
            label_list.append(line)
    my_file.close()
    return cluster_list, label_list
    
    
def test_index_zero(line, string):
    '''Returns True if index zero in a given list is equal to a specified
    string.
    '''
    if line != '':
        columns = line.split()
        if columns[0] == string:
            return True
    return False
    

def test_if_ngc(line):
    '''Determines if a line of data is from the New General Catalogue. If the 
    ID of the data begins with NGC the function returns True, otherwise returns 
    False.
    '''
    return test_index_zero(line, 'NGC')
    
assert test_if_ngc('NGC 101') == True
assert test_if_ngc('Pal 1000') == False
    
    
def test_if_label(line):
    '''Determines if a line of data is a header label for the data below. 
    Returns true if the line begins with ID, otherwise returns False.
    '''
    return test_index_zero(line, 'ID')
    
    
assert test_if_label('ID 101 4 AD') == True
assert test_if_label('NGC 10') == False
    
def find_id_set(str_list):
    '''Takes a list of strings that represent lines of data. Finds and returns
    a set of all globular cluster IDs within the given list.
    '''
    id_set = set([])
    for line in str_list:
        columns = line.split()
        id_set.add(columns[0] + ' ' + columns[1])
    return id_set
    
assert find_id_set(['a 1', 'b 2', 'c 3', 'd 4 5']) == \
set(['a 1', 'b 2', 'c 3', 'd 4'])
    
    
def make_cluster_nested_dict(id_set):
    '''Given a set of strings, returns a dictionary where each string is a key
    that maps to another dictionary.
    '''
    cluster_dict = {}
    for name in id_set:
        cluster_dict[name] = {}
    return cluster_dict


def update_1(data, cluster_dict):
    '''Given data in a list of strings and a nested dictionary, returns an 
    updated version of the nested dictionary. This update adds distance from the
    Sun to the nested dictionary for each NGC globular cluster.
    '''
    for line in data:
        columns = line.split()
        name = columns[0] + ' ' + columns[1]
        if len(columns) == 17:
            cluster_dict[name]['R_Sun'] = float(columns[12])
        else:
            cluster_dict[name]['R_Sun'] = float(columns[10])
    return cluster_dict
    
    
def update_2(data, cluster_dict):
    '''Given data in a list of strings and a nested dictionary, returns an 
    updated version of the nested dictionary. This update adds metallicity, 
    metallicity weight, absolute visual magnitude, and color magnitude values 
    to the nested dictionary for each NGC globular cluster.
    '''
    for line in data:
        columns = line.split()
        name = columns[0] + ' ' + columns[1]
        cluster_dict[name]['[Fe/H]'] = float(columns[2])
        cluster_dict[name]['wt'] = int(columns[3])
        cluster_dict[name]['V'] = float(columns[8])
        if name == 'NGC 5466' or name == 'NGC 6380':
            cluster_dict[name]['B-V'] = float(columns[9])
        elif name != 'NGC 6540':
            cluster_dict[name]['B-V'] = float(columns[10])
            cluster_dict[name]['U-B'] = float(columns[9])
    return cluster_dict
    
    
def update_3(data, cluster_dict):
    '''Given data in a list of strings and a nested dictionary, returns an 
    updated version of the nested dictionary. This update adds radial velocity
    values to the nested dictionary for each NGC globular cluster.
    '''
    for line in data:
        columns = line.split()
        name = columns[0] + ' ' + columns[1]
        cluster_dict[name]['v_r'] = float(columns[2])
    return cluster_dict
    
    
def add_data_to_dict(data_list, cluster_dict):
    '''Given a list of lists that contain strings of data and a dictionary
    where clusters are keys that map to empty dictionaries, returns a nested
    dictionary containing labeled data.
    '''
    cluster_dict = update_1(data_list[0], cluster_dict)
    cluster_dict = update_2(data_list[1], cluster_dict)
    cluster_dict = update_3(data_list[2], cluster_dict)
    return cluster_dict
    

def split_data_set(str_list, size):
    '''Splits up a list of strings into a list of lists. Each list in the list
    of lists is of equal length and contains subsequent lines.
    
    Parameters:
        str_list: A list of strings.
        size: Number of equal length lists to be contained in list of lists.
            
    Returns:
        list_of_lists: A list of lists of strings.
    '''
    list_of_lists = []
    length = len(str_list)
    for num in range(size):
        lst = []
        for i in range((num) * (length / size), (num + 1) * (length / size)):
            lst.append(str_list[i])
        list_of_lists.append(lst)
    return list_of_lists
            
assert split_data_set(['a', 'b', 'c', 'd', 'e', 'f'], 3) == \
[['a', 'b'], ['c', 'd'], ['e', 'f']] 
assert split_data_set(['a', 'b', 'c', 'd', 'e', 'f'], 2) == \
[['a', 'b', 'c'], ['d', 'e', 'f']]
            
    
def metallicity_weights_histogram(weights_list):
    '''Takes a list of metallicities and returns list integers where each index
    of the list corresponds to number of clusters with that weight.
    '''
    histogram = [0] * 15
    for metal_value in weights_list:
        histogram[metal_value] += 1
    return histogram
    
    
def find_mean(val_list):
    '''Given a list of ints or floats, returns the mean float of the list.
    '''
    length = len(val_list)
    list_sum = 0
    for val in val_list:
        list_sum += val
    return float(list_sum) / length
    
assert find_mean([4.0, 1.0, 3.0, 8.0]) == 4.0


def print_key():
    '''Prints a key that explains abbreviations and gives units for each.
    '''
    print 'Key to Data Values:'
    print 'v_r: Radial velocity in km/s'
    print 'V: Absolute visual magnitude'
    print 'R_Sun: Distance from cluster to the Sun in Kiloparsecs'
    print 'wt: Number of measurements in [Fe/H] value'
    print 'B-V: B-V color magnitude'
    print 'U-B: U-B color magnitude'
    print '[Fe/H]: Metallicity ratio value'
    print ''


def print_mean_vals(nested_dict):
    '''Given a nested dictionary, prints the mean values for each quantity
    averaged over all included NGC globular clusters.
    '''
    clusters = nested_dict.keys()
    label_list = nested_dict[clusters[0]].keys()
    for label in label_list:
        val_list = []
        for cluster in clusters:
            if label in nested_dict[cluster].keys():
                val_list.append(nested_dict[cluster][label])
        mean = find_mean(val_list)
        print label + ':', mean
        
        
def percent_clusters(nested_dict):
    '''Given a nested dictionary, this function returns the decimal percent of
    globular clusters in the dictionary that are moving toward Earth.
    '''
    clusters = nested_dict.keys()
    toward_earth = 0
    label = 'v_r'
    for cluster in clusters:
        v = nested_dict[cluster][label]
        if v < 0:
            toward_earth += 1
    return float(toward_earth) / len(clusters)
    
        
def extract_color_mag(nested_dict):
    '''Given a nested dictionary, returns lists of U-B color magnitudes, 
    B-V color magnitudes, and absolute visual magnitudes (V) for the NGC 
    globular clusters.
    '''
    clusters = nested_dict.keys()
    v_list = []
    bv_color_list = []
    ub_color_list = []
    for cluster in clusters:
        if 'B-V' in nested_dict[cluster].keys() and \
        'U-B' in nested_dict[cluster].keys():
            bv_color_list.append(nested_dict[cluster]['B-V'])
            v_list.append(nested_dict[cluster]['V'])
            ub_color_list.append(nested_dict[cluster]['U-B'])
    return ub_color_list, bv_color_list, v_list
    
    
def extract_distance_velocity_wt(nested_dict):
    '''Given a nested dictionary, returns lists of distances from the Sun,
    radial velocities, and metallicity weights for the NGC globular clusters.
    '''
    clusters = nested_dict.keys()
    r_sun = []
    v_r = []
    wt = []
    for cluster in clusters:
        r_sun.append(nested_dict[cluster]['R_Sun'])
        v_r.append(nested_dict[cluster]['v_r'])
        wt.append(nested_dict[cluster]['wt'])
    return r_sun, v_r, wt
    
   
def plot_color_mag_diagram(color_list, mag_list, save_file='colormag.png'):
    '''Creates and saves a scatter plot of the absolute visual magnitude 
    plotted against the B-V color of NGC globular clusters.
    
    Parameters:
        color_list: List of floats representing B-V colors of NGC 
        globular clusters.
        mag_list: List of floats representing absolute visual magnitudes of NGC
        globular clusters
        save_file: String representing save path for plot. 
        Defaults to colormag.png
    '''
    plt.clf()
    plt.scatter(color_list, mag_list)
    plt.xlabel('Color (B-V)')
    plt.ylabel('Magnitude (V)')
    plt.title('Color Magnitude Diagram for NGC Globular Clusters')
    plt.savefig(save_file)
    plt.clf()
    
    
def plot_color_color_diagram(ub_color, bv_color, save_file='colorcolor.png'):
    '''Creates and saves  plot of the U-B color magnitude against the B-V color 
    magnitude for NGC globular clusters. Default saves to colorcolor.png.
    '''
    plt.clf()
    plt.scatter(ub_color, bv_color)
    plt.xlabel('Color (U-B)')
    plt.ylabel('Color (B-V)')
    plt.title('Color Diagram for NGC Globular Clusters')
    plt.savefig(save_file)
    plt.clf()
    
    
def plot_distance_velocity(distance, v_r, save_file='velocity.png'):
    '''Creates and saves  plot of the distance from the Sun (R_Sun) against
    the radial velocity of the cluster (v_r). Default saves to velocity.png.
    '''
    plt.clf()
    plt.scatter(distance, v_r)
    plt.xlabel('Distance from Sun (Kiloparsecs)')
    plt.ylabel('Radial Velocity (km/s)')
    plt.title('Relation between NGC Radial Velocities and Distances')
    plt.savefig(save_file)
    plt.clf()
    
    
def plot_histogram(numbers, title, x_label, y_label, save_file='histogram.png'):
    '''Creates and saves a plot of a histogram to default image histogram.png.
    
    Parameters:
        title: String representing title of histogram.
        x_label: String representing label of x-axis of histogram.
        y_label: String representing label of y-axis of histogram.
        save_file: String representing save path for image of histogram. 
        Defaults to histogram.png.
    '''
    plt.clf()
    plt.plot(numbers)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.savefig(save_file)
    plt.clf()

    
def main():
    cluster_list, label_list = \
    read_cluster_dat_file('http://physwww.physics.mcmaster.ca/~harris/mwgc.dat')
    id_set = find_id_set(cluster_list)
    split_data_list = split_data_set(cluster_list, 3)
    nested_dict = make_cluster_nested_dict(id_set)
    nested_dict = add_data_to_dict(split_data_list, nested_dict)
    print_key()
    print 'Mean Values for Globular Clusters from NGC:'
    print_mean_vals(nested_dict)
    print ''
    percent = percent_clusters(nested_dict)
    print 'Statistics:'
    print 'Fraction of clusters moving toward Earth:', percent
    ub_color, bv_color, v_mag = extract_color_mag(nested_dict)
    r_sun, v_r, wt = extract_distance_velocity_wt(nested_dict)
    histogram = metallicity_weights_histogram(wt)
    plot_color_mag_diagram(bv_color, v_mag)
    plot_color_color_diagram(ub_color, bv_color)
    plot_distance_velocity(r_sun, v_r)
    plot_histogram(histogram, 'Metallicity Measurement Weights',\
    'Number of Independent Measurements', 'Number of NGC Globular Clusters')
    print 'Correlation and p-value between V and B-V:', \
    pearsonr(v_mag, bv_color)
    print 'Correlation and p-value between B-V and U-B:', \
    pearsonr(bv_color, ub_color)
    print 'Correlation and p-value between R_Sun and v_r:', pearsonr(r_sun, v_r)
    
        
if __name__ == '__main__':
    main()
