try:
    import structure_analysis_utils as utils
except ImportError:
    from . import structure_analysis_utils as utils

import math

filename = input("Enter filename: ")
molecule = utils.open_file(filename)
torsion_array = utils.array_of_carbon_nitrogen_torsions(molecule)

avg_torsion = utils.find_average_torsion(torsion_array)
#print("Mean: ", avg_torsion)

torsion_variance = utils.find_torsion_circular_variance(torsion_array)
#print("Variance: ", torsion_variance)

array_nitrogens = utils.array_of_imine_nitrogens(molecule)
imine_torsions = utils.array_of_imine_torsions(array_nitrogens, molecule)
imine_torsion_mean = utils.find_torsion_circular_mean(imine_torsions)
imine_torsion_variance = utils.find_torsion_circular_variance(imine_torsions)
#print(imine_torsions)
print("All layers")
print("Mean: ", imine_torsion_mean)
print("Variance: ", imine_torsion_variance)
print("")

# utils.plot_torsion_angles(torsion_array)
separate_molecules = utils.array_of_separate_molecules(molecule)
fst_list = ["First", "Second", "Third"]
for key, each in enumerate(separate_molecules):
    array_imine_nitrogens = utils.array_of_imine_nitrogens(each)
    array_imine_torsions = utils.array_of_imine_torsions(array_imine_nitrogens, each)
    imine_torsion_mean = utils.find_torsion_circular_mean(array_imine_torsions)
    imine_torsion_variance = utils.find_torsion_circular_variance(array_imine_torsions)
    print(fst_list[key], "Layer")
    print("Mean: ", imine_torsion_mean)
    print("Variance: ", imine_torsion_variance)
    print("")
    
array_carbons = utils.array_of_imine_carbons(molecule)

in_layer_NC_imine_distances = utils.get_in_layer_NC_imine_distances(array_nitrogens, array_carbons)
in_layer_imine_distance_mean = utils.get_interlayer_imine_distance_mean(in_layer_NC_imine_distances)
in_layer_imine_distance_var = utils.get_interlayer_imine_distance_variance(in_layer_NC_imine_distances)
print("C-N In-layer Distance")
print("Mean: ", in_layer_imine_distance_mean)
print("Variance: ", in_layer_imine_distance_var)
print("")

interlayer_NC_imine_distances = utils.get_interlayer_NC_imine_distances(array_nitrogens, array_carbons)
interlayer_NC_imine_distance_mean = utils.get_interlayer_imine_distance_mean(interlayer_NC_imine_distances)
interlayer_NC_imine_distance_var = utils.get_interlayer_imine_distance_variance(interlayer_NC_imine_distances)
print("C-N Interlayer Distance")
print("Mean: ", interlayer_NC_imine_distance_mean)
print("Variance: ", interlayer_NC_imine_distance_var)
print("")

interlayer_NN_imine_distances = utils.get_interlayer_NN_imine_distances(array_nitrogens, array_nitrogens)
interlayer_NN_imine_distance_mean = utils.get_interlayer_imine_distance_mean(interlayer_NN_imine_distances)
interlayer_NN_imine_distance_var = utils.get_interlayer_imine_distance_variance(interlayer_NN_imine_distances)

print("N-N Interlayer Distance")
print("Mean: ", interlayer_NN_imine_distance_mean)
print("Variance: ", interlayer_NN_imine_distance_var)
print("")

interlayer_CC_imine_distances = utils.get_interlayer_CC_imine_distances(array_carbons, array_carbons)
interlayer_CC_imine_distance_mean = utils.get_interlayer_imine_distance_mean(interlayer_CC_imine_distances)
interlayer_CC_imine_distance_var = utils.get_interlayer_imine_distance_variance(interlayer_CC_imine_distances)

print("C-C Interlayer Distance")
print("Mean: ", interlayer_CC_imine_distance_mean)
print("Variance: ", interlayer_CC_imine_distance_var)
print("")

framework_offset, framework_separation = utils.get_framework_offset(in_layer_imine_distance_mean, \
                            interlayer_NN_imine_distance_mean, interlayer_NC_imine_distance_mean)
framework_skew = utils.get_framework_skew(molecule, array_nitrogens, array_carbons)
framework_skew_mean = utils.get_framework_skew_mean(framework_skew)
framework_skew_var = utils.get_framework_skew_var(framework_skew)

print("Framework Separation ", framework_separation)
print("Framework Offset ", framework_offset)
print("Framework Skew Mean ", framework_skew_mean)
print("Framework Skew Variance ", framework_skew_var)
print("")
print(framework_skew)

    