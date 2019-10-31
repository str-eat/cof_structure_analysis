from openbabel import openbabel, pybel
import matplotlib.pyplot as plt
from scipy.stats import circvar, circmean
import numpy as np
from math import sqrt


def open_file(filename):
    molecule = pybel.readfile("xyz", filename).__next__().OBMol
    return molecule

def find_average_torsion(array_of_torsions):
    torsion_array = array_of_torsions
    average_torsion = torsion_array.mean()
    return average_torsion

def find_torsion_circular_variance(array_of_torsions):
    torsion_array = array_of_torsions
    torsion_variance = circvar(torsion_array)
    return torsion_variance

def find_torsion_circular_mean(array_of_torsions):
    torsion_mean = circmean(array_of_torsions)
    return torsion_mean

def plot_torsion_angles(array_of_torsions):
    torsion_array = array_of_torsions
    print(torsion_array[:10])
    plt.hist(torsion_array, bins=50)
    plt.show()

def array_of_separate_molecules(OBMolecule):
    molecule = OBMolecule
    separate_molecules = list(molecule.Separate())
    return separate_molecules

def array_of_carbon_nitrogen_torsions(OBMolecule):
    molecule = OBMolecule
    torsion_list = list()
    for torsion in openbabel.OBMolTorsionIter(molecule):
        if 0 in torsion:
            continue
        atom1 = molecule.GetAtom(torsion[0])
        atom2 = molecule.GetAtom(torsion[1])
        atom3 = molecule.GetAtom(torsion[2])
        atom4 = molecule.GetAtom(torsion[3])
        if {atom1.GetAtomicNum(), atom2.GetAtomicNum(), atom3.GetAtomicNum(), atom4.GetAtomicNum()} == {6, 7}:
            angle = molecule.GetTorsion(atom1, atom2, atom3, atom4)            
            torsion_list.append(angle)
    torsion_array = np.array(torsion_list, dtype='float')
    return torsion_array

def array_of_imine_nitrogens(OBMolecule):
    molecule = OBMolecule
    imine_nitrogens = list()
    for atom in openbabel.OBMolAtomIter(molecule):
        atomic_num = atom.GetAtomicNum()
        atom_in_ring_size = atom.MemberOfRingSize()
        if atomic_num == 7:
            if atom_in_ring_size == 0:
                imine_nitrogens.append(atom)
    return imine_nitrogens

def array_of_imine_carbons(OBMolecule):
    molecule = OBMolecule
    imine_carbons = list()
    for atom in openbabel.OBMolAtomIter(molecule):
        atomic_num = atom.GetAtomicNum()
        atom_in_ring_size = atom.MemberOfRingSize()
        if atomic_num == 6:
            if atom_in_ring_size == 0:
                imine_carbons.append(atom)
    return imine_carbons

def array_of_imine_torsions(imine_nitrogens, molecule):
    
    def beta_atoms(atom, nitrogen, imine, in_ring, begin_end):
        for bond in openbabel.OBAtomBondIter(atom):
            
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            if begin_atom != nitrogen and begin_atom != atom:
                if in_ring == True and begin_end is "begin":
                    imine.insert(0, begin_atom.GetIdx())
                elif in_ring == True and begin_end is "end":
                    imine.append(begin_atom.GetIdx())
                elif in_ring == False and begin_end is "begin":
                    imine.append(begin_atom.GetIdx())                    
                else:
                    imine.insert(0, begin_atom.GetIdx())
                return imine

            if end_atom != nitrogen and end_atom != atom:
                if in_ring == True and begin_end is "begin":
                    imine.insert(0, end_atom.GetIdx())
                elif in_ring == True and begin_end is "end":
                    imine.append(end_atom.GetIdx())
                elif in_ring == False and begin_end is "begin":
                    imine.append(end_atom.GetIdx())
                else:
                    imine.insert(0, begin_atom.GetIdx())
                return imine
        return imine
    
    def get_torsions(imine_struct):                  
        torsion_list = list()
        for torsion in imine_struct:                   
            atom1 = molecule.GetAtom(torsion[0])
            atom2 = molecule.GetAtom(torsion[1])
            atom3 = molecule.GetAtom(torsion[2])
            atom4 = molecule.GetAtom(torsion[3])        
            angle = molecule.GetTorsion(atom1, atom2, atom3, atom4)            
            torsion_list.append(angle)
        torsion_array = np.array(torsion_list, dtype='float')
        return torsion_array        

    imine_struct = list()
    for nitrogen in imine_nitrogens:        
        imine = list()
        imine.append(nitrogen.GetIdx())

        for bond in openbabel.OBAtomBondIter(nitrogen):            
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            if begin_atom != nitrogen:                
                in_ring = begin_atom.MemberOfRingSize() != 0
                begin_end = "begin"
                if in_ring == True:
                    imine.insert(0, begin_atom.GetIdx())
                    imine = beta_atoms(begin_atom, nitrogen, imine, in_ring, begin_end)
                else:
                    imine.append(begin_atom.GetIdx())

            if end_atom != nitrogen:                
                in_ring = end_atom.MemberOfRingSize() != 0
                begin_end = "end"                
                if in_ring == True:
                    imine.append(end_atom.GetIdx())
                    imine = beta_atoms(end_atom, nitrogen, imine, in_ring, begin_end)
                else:
                    imine.insert(0, end_atom.GetIdx())                    

        imine_struct.append(imine)
    torsion_list = get_torsions(imine_struct)
    return torsion_list

def get_in_layer_NC_imine_distances(imine_nitrogens, imine_carbons):
    distance_list = []
    for idx, nitrogen in enumerate(imine_nitrogens):
        distance_loop = []
        for carbon in imine_carbons:            
            distance = nitrogen.GetDistance(carbon)
            distance_loop.append(distance)
        distance_loop.sort()
        distance_list.append(distance_loop[0])
    return distance_list

def get_interlayer_NC_imine_distances(imine_nitrogens, imine_carbons):
    distance_list = []
    for idx, nitrogen in enumerate(imine_nitrogens):
        distance_loop = []
        for carbon in imine_carbons:            
            distance = nitrogen.GetDistance(carbon)
            distance_loop.append(distance)
        distance_loop.sort()
        distance_list.append(distance_loop[1])
    return distance_list

def get_interlayer_NN_imine_distances(imine_nitrogens1, imine_nitrogens2):
    distance_list = []
    for idx, nitrogen1 in enumerate(imine_nitrogens1):
        distance_loop = []
        for nitrogen2 in imine_nitrogens2:            
            if nitrogen1 != nitrogen2:
                distance = nitrogen1.GetDistance(nitrogen2)
                distance_loop.append(distance)
        distance_loop.sort()
        distance_list.append(distance_loop[0])
    return distance_list

def get_interlayer_CC_imine_distances(imine_carbons1, imine_carbons2):
    distance_list = []
    for idx, carbon1 in enumerate(imine_carbons1):
        distance_loop = []
        for carbon2 in imine_carbons2:
            if carbon1 != carbon2:
                distance = carbon1.GetDistance(carbon2)
                distance_loop.append(distance)
        distance_loop.sort()
        distance_list.append(distance_loop[0])
    return distance_list

def get_interlayer_imine_distance_mean(distance_list):
    distance_mean = sum(distance_list)/len(distance_list)
    return distance_mean

def get_interlayer_imine_distance_variance(distance_list):
    distance_variance = np.var(distance_list)
    return distance_variance

def get_framework_offset(in_layer_CN_distance, interlayer_NN_distance, interlayer_CN_distance):
    half_perimeter = (in_layer_CN_distance + interlayer_NN_distance + interlayer_CN_distance) / 2
    area = sqrt(half_perimeter * (half_perimeter - in_layer_CN_distance) * (half_perimeter - interlayer_NN_distance) \
                        * (half_perimeter - interlayer_CN_distance))
    height = (area * 2) / in_layer_CN_distance
    offset_distance = sqrt((interlayer_NN_distance * interlayer_NN_distance) - (height * height))
    return offset_distance, height

def get_framework_skew(molecule, array_nitrogens, array_carbons):
    skew = list()
    NC_distance_list = list()
    for nitrogen in array_nitrogens:
        distance_loop = list()
        for carbon in array_carbons:
            distance = nitrogen.GetDistance(carbon)
            distance_loop.append((nitrogen.GetIdx(), carbon.GetIdx(), distance))
        NC_distance_list.append((min(distance_loop, key=lambda x:x[2])))

    NN_distance_list = list()
    for nitrogen1 in array_nitrogens:
        distance_loop = list()
        for nitrogen2 in array_nitrogens:
            if nitrogen1 != nitrogen2:
                distance = nitrogen1.GetDistance(nitrogen2)
                distance_loop.append((nitrogen1.GetIdx(), nitrogen2.GetIdx(), distance))
        NN_distance_list.append((min(distance_loop, key=lambda x:x[2])))

    for nc in NC_distance_list:
        for nn in NN_distance_list:
            if nc[0] == nn[0]:
                for match in NC_distance_list:
                    if nn[1] == match[0]:
                        skew.append(molecule.GetTorsion(nc[1], nc[0], match[0], match[1]))
    return skew

def get_framework_skew_mean(skew_list):
    skew_mean = circmean(skew_list)
    return skew_mean

def get_framework_skew_var(skew_list):
    skew_variance = circvar(skew_list)
    return skew_variance


