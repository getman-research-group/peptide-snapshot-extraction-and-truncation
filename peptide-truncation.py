# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 16:34:14 2023

@author: biswas.161
"""

from ase import io, geometry, Atoms
from ase.io import gromacs, proteindatabank
from ase.build import sort
from ase.geometry.analysis import Analysis
import numpy as np
import math
import os

####Input snapshot number
snapshot = 1 

####Input Ln ion
Ln = "La"

LanM1_Ln_pdb = gromacs.read_gromacs("md{}.gro".format(snapshot)) #read gro file based on snapshot

#Creating lists to read ASE Atoms object
atom_symbols = LanM1_Ln_pdb.get_chemical_symbols() #list of symbols
atom_coordinates = LanM1_Ln_pdb.get_positions() #array of coordinates
atom_coordinates_list = atom_coordinates.tolist() #converting array to list of coordinates
resname = LanM1_Ln_pdb.get_array('residuenames')
resname_list = resname.tolist()
atomtype = LanM1_Ln_pdb.get_array('atomtypes')
atomtype_list = atomtype.tolist()
resnum = LanM1_Ln_pdb.get_array('residuenumbers')
resnum_list = resnum.tolist()

gromacs.write_gromacs("md_test.gro", LanM1_Ln_pdb)

atom_list = zip(atom_symbols, atomtype_list, resname_list, atom_coordinates_list, resnum_list) #zipping together symbol and coordinates list
zipped_list = list(atom_list) #converting zip to list format

#Creating empty lists to populate with truncated lists
new_zipped_list_pos = []
new_zipped_list_sym = []
new_zipped_list_resname = []
new_zipped_list_atomtype = []
new_zipped_list_resnum = []

zl_pos = []
zl_sym = []
zl_resname = []
zl_atomtype = []
zl_resnum = []

#Position/coordinates of the lanthanide (last position in pdb file)
Ln_zl = zipped_list[-1][3]

####Input cutoff distance
water_cutoff = 6

#Loop through entire list of atoms in pdb file and add waters within input water cutoff               
for i in range((len(zipped_list)-1)): 
    
    if zipped_list[i][0] != 'O': 
        
        if zipped_list[i][0] != 'H':
        
            new_zipped_list_pos.append(zipped_list[i][3])
            new_zipped_list_sym.append(zipped_list[i][0])
            new_zipped_list_resname.append(zipped_list[i][2])
            new_zipped_list_atomtype.append(zipped_list[i][1])
            new_zipped_list_resnum.append(zipped_list[i][4])
            
        if zipped_list[i][0] == 'H':
            
            if zipped_list[i-1][0] == 'O' and zipped_list[i+1][0] == 'H': #if O 'H' H
                
                if math.dist(zipped_list[i-1][3], Ln_zl) <= water_cutoff:
                    
                    new_zipped_list_pos.append(zipped_list[i][3])
                    new_zipped_list_sym.append(zipped_list[i][0])
                    new_zipped_list_resname.append(zipped_list[i][2])
                    new_zipped_list_atomtype.append(zipped_list[i][1])
                    new_zipped_list_resnum.append(zipped_list[i][4])
                    zl_pos.append(zipped_list[i][3])
                    zl_sym.append(zipped_list[i][0])
                    zl_resname.append(zipped_list[i][2])
                    zl_atomtype.append(zipped_list[i][1])
                    zl_resnum.append(zipped_list[i][4])
                    
            elif zipped_list[i-2][0] == 'O' and zipped_list[i-1][0] == 'H': #if O H 'H'
                
                if math.dist(zipped_list[i-2][3], Ln_zl) <= water_cutoff:
                    
                    new_zipped_list_pos.append(zipped_list[i][3])
                    new_zipped_list_sym.append(zipped_list[i][0])
                    new_zipped_list_resname.append(zipped_list[i][2])
                    new_zipped_list_atomtype.append(zipped_list[i][1])
                    new_zipped_list_resnum.append(zipped_list[i][4])
                    zl_pos.append(zipped_list[i][3])
                    zl_sym.append(zipped_list[i][0])
                    zl_resname.append(zipped_list[i][2])
                    zl_atomtype.append(zipped_list[i][1])
                    zl_resnum.append(zipped_list[i][4])
                    
            elif (zipped_list[i-1][0] != 'O' and zipped_list[i+1][0] != 'H') or (zipped_list[i-2][0] != 'O' and zipped_list[i-1][0] != 'H'): #if non OHH 'O' or 'H'    
                
                new_zipped_list_pos.append(zipped_list[i][3])
                new_zipped_list_sym.append(zipped_list[i][0])
                new_zipped_list_resname.append(zipped_list[i][2])
                new_zipped_list_atomtype.append(zipped_list[i][1])
                new_zipped_list_resnum.append(zipped_list[i][4])

                
    if zipped_list[i][0] == 'O':
        
        if (zipped_list[i+1][0] != 'H' or zipped_list[i+2][0] != 'H'): 
            
            new_zipped_list_pos.append(zipped_list[i][3])
            new_zipped_list_sym.append(zipped_list[i][0])
            new_zipped_list_resname.append(zipped_list[i][2])
            new_zipped_list_atomtype.append(zipped_list[i][1])
            new_zipped_list_resnum.append(zipped_list[i][4])
            
        elif math.dist(zipped_list[i][3], Ln_zl) <= water_cutoff:
            
            new_zipped_list_pos.append(zipped_list[i][3])
            new_zipped_list_sym.append(zipped_list[i][0])
            new_zipped_list_resname.append(zipped_list[i][2])
            new_zipped_list_atomtype.append(zipped_list[i][1])
            new_zipped_list_resnum.append(zipped_list[i][4])
            
            if (zipped_list[i+1][0] == 'H' and zipped_list[i+2][0] == 'H'):
                
                zl_pos.append(zipped_list[i][3])
                zl_sym.append(zipped_list[i][0])
                zl_resname.append(zipped_list[i][2])
                zl_atomtype.append(zipped_list[i][1])
                zl_resnum.append(zipped_list[i][4])
    
#Adding Lns to list
new_zipped_list_pos.append(zipped_list[-1][3])
new_zipped_list_sym.append(zipped_list[-1][0]) 
new_zipped_list_resname.append(zipped_list[-1][2])
new_zipped_list_atomtype.append(zipped_list[-1][1])
new_zipped_list_resnum.append(zipped_list[-1][4])

#Update arrays with Ln
resname_array = np.array(new_zipped_list_resname)
atomtype_array = np.array(new_zipped_list_atomtype)
resnum_array = np.array(new_zipped_list_resnum)

LanM1_Ln_pdb_new = Atoms(new_zipped_list_sym, new_zipped_list_pos)
LanM1_Ln_pdb_new.set_array('residuenames', resname_array)
LanM1_Ln_pdb_new.set_array('atomtypes', atomtype_array)
LanM1_Ln_pdb_new.set_array('residuenumbers', resnum_array)

#Convert pdb file to gromacs type file
gromacs.write_gromacs("LanM1_pdb_new.gro", LanM1_Ln_pdb_new)

LanM1_Ln_pdb_new = gromacs.read_gromacs("LanM1_pdb_new.gro")
charge = LanM1_Ln_pdb_new.get_initial_charges()
#Creating list to read ASE Atoms object
atom_symbols_new = LanM1_Ln_pdb_new.get_chemical_symbols() #list of symbols
atom_coordinates_new = LanM1_Ln_pdb_new.get_positions() #array of coordinates
atom_coordinates_list_new = atom_coordinates_new.tolist() #converting array to list of coordinates

resname_new = LanM1_Ln_pdb_new.get_array('residuenames')
resname_list_new = resname_new.tolist()
atomtype_new = LanM1_Ln_pdb_new.get_array('atomtypes')
atomtype_list_new = atomtype_new.tolist()
resnum_new = LanM1_Ln_pdb_new.get_array('residuenumbers')
resnum_list_new = resnum_new.tolist()

atom_list_new = zip(atom_symbols_new, atomtype_list_new, resname_list_new, atom_coordinates_list_new, resnum_list_new) #zipping together symbol and coordinates list
zipped_list_new = list(atom_list_new) #converting zip to list format

#Creating empty lists to populate with truncated lists
zipped_list_pos_new = []
zipped_list_sym_new = []
zipped_list_resname_new = []
zipped_list_atomtype_new = []
zipped_list_resnum_new = []

#Position/coordinates of the lanthanide (last position in pdb file)
Ln_zl_new = zipped_list_new[-1][3]

#### Input peptide cutoff distance
cutoff = 6 #4.8

ana = Analysis(LanM1_Ln_pdb_new)

#Get number of bonds, lengths, max and min lengths for C-C, C-N, C-O, C-H and N-H
C_C_bonds = ana.get_bonds('C', 'C', unique=True)
C_N_bonds = ana.get_bonds('C', 'N', unique=True)
C_O_bonds = ana.get_bonds('C', 'O', unique=True)
C_H_bonds = ana.get_bonds('C', 'H', unique=True)
C_H_length = ana.get_values(C_H_bonds)
C_H_length_max = max(C_H_length[0])
C_H_length_avg = np.average(C_H_length[0])
N_H_bonds = ana.get_bonds('N', 'H', unique=True)
O_H_bonds = ana.get_bonds('O', 'H', unique=True)

C_C_x = [] #'C' in 'C'-C bond
C_C_y = [] #'C' in C-'C' bond

#Separate ('C'==x, 'C'==y) from C_C_bonds list into 2 different lists
for i1 in range((len(C_C_bonds[0]))):
    
    C_C_x.append(C_C_bonds[0][i1][0])
    C_C_y.append(C_C_bonds[0][i1][1]) 
    
C_N_x = [] #C in C-N bond
C_N_y = [] #N in C-N bond

#Separate ('C'==x, 'N'==y) from C_N_bonds list into 2 different lists    
for i2 in range((len(C_N_bonds[0]))):
    
    C_N_x.append(C_N_bonds[0][i2][0]) #all the Cs connected to Ns
    C_N_y.append(C_N_bonds[0][i2][1]) #all the Ns
    
C_O_x = []
C_O_y = []

#Separate ('C'==x, 'O'==y) from C_O_bonds list into 2 different lists  
for i3 in range((len(C_O_bonds[0]))):
    
    C_O_x.append(C_O_bonds[0][i3][0]) #all the Cs connected to Os
    C_O_y.append(C_O_bonds[0][i3][1]) #all the Os

C_H_x = []
C_H_y = []

#Separate ('C'==x, 'H'==y) from C_H_bonds list into 2 different lists    
for i4 in range((len(C_H_bonds[0]))):
    
    C_H_x.append(C_H_bonds[0][i4][0]) #all the Cs connected to Hs
    C_H_y.append(C_H_bonds[0][i4][1]) #all the Hs
    
N_H_x = []
N_H_y = []

#Separate ('N'==x, 'H'==y) from N_H_bonds list into 2 different lists    
for i5 in range((len(N_H_bonds[0]))):
    
    N_H_x.append(N_H_bonds[0][i5][0]) #all the Cs connected to Hs
    N_H_y.append(N_H_bonds[0][i5][1]) #all the Hs
    
O_H_x = []
O_H_y = []

#Separate ('O'==x, 'H'==y) from O_H_bonds list into 2 different lists    
for i8 in range((len(O_H_bonds[0]))):
    
    O_H_x.append(O_H_bonds[0][i8][0]) #all the Os connected to Hs
    O_H_y.append(O_H_bonds[0][i8][1]) #all the Hs
 
N_new = [] #Account for the Ns included in the structure
C_C_new = [] #Add C of C-C bonds to this list if the attached N makes the cutoff
C_N_new = []
C_O_new = []
C_H_new = []
N_H_new = []
O_new = []

#Checking distance between all C atoms and Ln and adding to list if distance is <=cutoff - no duplicates   
for i in range(len(C_C_x)):
    
    if math.dist(zipped_list_new[C_C_x[i]][3], Ln_zl_new) <= cutoff: #checking 'C'-C
        
        if zipped_list_new[C_C_x[i]][3] not in zipped_list_pos_new:
            
            zipped_list_pos_new.append(zipped_list_new[C_C_x[i]][3])
            zipped_list_sym_new.append(zipped_list_new[C_C_x[i]][0])
            zipped_list_resname_new.append(zipped_list_new[C_C_x[i]][2])
            zipped_list_atomtype_new.append(zipped_list_new[C_C_x[i]][1])
            zipped_list_resnum_new.append(zipped_list_new[C_C_x[i]][4])
            C_C_new.append(C_C_x[i])
        
    if math.dist(zipped_list_new[C_C_y[i]][3], Ln_zl_new) <= cutoff: #checking C-'C'
        
        if zipped_list_new[C_C_y[i]][3] not in zipped_list_pos_new:
            
            zipped_list_pos_new.append(zipped_list_new[C_C_y[i]][3])
            zipped_list_sym_new.append(zipped_list_new[C_C_y[i]][0])
            zipped_list_resname_new.append(zipped_list_new[C_C_y[i]][2])
            zipped_list_atomtype_new.append(zipped_list_new[C_C_y[i]][1])
            zipped_list_resnum_new.append(zipped_list_new[C_C_y[i]][4])
            C_C_new.append(C_C_y[i])

for k in range(len(C_N_x)): 
    
    if C_N_x[k] in C_C_new: #checking if C of C-N is part of C of C-C making cutoff - if yes, include the N
        
        if C_N_y[k] not in N_new:
            
            zipped_list_pos_new.append(zipped_list_new[C_N_y[k]][3])
            zipped_list_sym_new.append(zipped_list_new[C_N_y[k]][0])
            zipped_list_resname_new.append(zipped_list_new[C_N_y[k]][2])
            zipped_list_atomtype_new.append(zipped_list_new[C_N_y[k]][1])
            zipped_list_resnum_new.append(zipped_list_new[C_N_y[k]][4])
            C_N_new.append(C_N_y[k])
            N_new.append(C_N_y[k])

    if C_N_x[k] in C_C_x or C_C_y: #checking if C of C-N is part of original C-bonds - if yes, next condition
             
             if math.dist(zipped_list_new[C_N_y[k]][3], Ln_zl_new) <= (cutoff): #Add Ns within cutoff region, attached to atleast 1 C of C-C
                 
                 if C_N_y[k] not in N_new:
                     
                     zipped_list_pos_new.append(zipped_list_new[C_N_y[k]][3])
                     zipped_list_sym_new.append(zipped_list_new[C_N_y[k]][0])
                     zipped_list_resname_new.append(zipped_list_new[C_N_y[k]][2])
                     zipped_list_atomtype_new.append(zipped_list_new[C_N_y[k]][1])
                     zipped_list_resnum_new.append(zipped_list_new[C_N_y[k]][4])
                     N_new.append(C_N_y[k])
                 
                 if math.dist(zipped_list_new[C_N_x[k]][3], Ln_zl_new) <= (cutoff+1.6): #Add Cs bonded to Ns 1 A beyond cutoff region (to make sure only C-C truncations)
                     
                     if C_N_x[k] not in C_C_new:
                         
                         zipped_list_pos_new.append(zipped_list_new[C_N_x[k]][3])
                         zipped_list_sym_new.append(zipped_list_new[C_N_x[k]][0])
                         zipped_list_resname_new.append(zipped_list_new[C_N_x[k]][2])
                         zipped_list_atomtype_new.append(zipped_list_new[C_N_x[k]][1])
                         zipped_list_resnum_new.append(zipped_list_new[C_N_x[k]][4])
                         C_C_new.append(C_N_x[k])
                 
for j in range(len(C_O_x)): #include all the Os where 'C'-O is included already

    if C_O_x[j] in C_C_new:
        
         zipped_list_pos_new.append(zipped_list_new[C_O_y[j]][3])
         zipped_list_sym_new.append(zipped_list_new[C_O_y[j]][0])
         zipped_list_resname_new.append(zipped_list_new[C_O_y[j]][2])
         zipped_list_atomtype_new.append(zipped_list_new[C_O_y[j]][1])
         zipped_list_resnum_new.append(zipped_list_new[C_O_y[j]][4])
         O_new.append(C_O_y[j])
         C_O_new.append(C_O_x[j])
            
for l in range(len(C_H_x)): #include all the Hs where 'C'-H is included already

    if C_H_x[l] in C_C_new:
          
          zipped_list_pos_new.append(zipped_list_new[C_H_y[l]][3])
          zipped_list_sym_new.append(zipped_list_new[C_H_y[l]][0])  
          zipped_list_resname_new.append(zipped_list_new[C_H_y[l]][2])
          zipped_list_atomtype_new.append(zipped_list_new[C_H_y[l]][1])
          zipped_list_resnum_new.append(zipped_list_new[C_H_y[l]][4])

for m in range(len(N_H_x)): #include all the Hs where 'N'-H is included already

    if N_H_x[m] in N_new:
          
          zipped_list_pos_new.append(zipped_list_new[N_H_y[m]][3])
          zipped_list_sym_new.append(zipped_list_new[N_H_y[m]][0])
          zipped_list_resname_new.append(zipped_list_new[N_H_y[m]][2])
          zipped_list_atomtype_new.append(zipped_list_new[N_H_y[m]][1])
          zipped_list_resnum_new.append(zipped_list_new[N_H_y[m]][4])
          N_H_new.append(N_H_x[m])
          
for s in range(len(O_H_x)): #include all the Hs where 'N'-H is included already

    if O_H_x[s] in O_new:
        
        if zipped_list_new[O_H_y[s]][3] not in zl_pos: 
            
          zipped_list_pos_new.append(zipped_list_new[O_H_y[s]][3])
          zipped_list_sym_new.append(zipped_list_new[O_H_y[s]][0]) 
          zipped_list_resname_new.append(zipped_list_new[O_H_y[s]][2])
          zipped_list_atomtype_new.append(zipped_list_new[O_H_y[s]][1])
          zipped_list_resnum_new.append(zipped_list_new[O_H_y[s]][4])

for p in range(len(C_C_x)): 
    
    if (C_C_x[p] in C_C_new) or (C_C_y[p] in C_C_new): #if 'C'-C or C-'C' is included in the added Cs
        
        if C_C_x[p] not in C_C_new: #
            
            zipped_list_pos_new.append(zipped_list_new[C_C_x[p]][3])
            zipped_list_sym_new.append('H')
            zipped_list_resname_new.append(zipped_list_new[C_C_x[p]][2])
            zipped_list_atomtype_new.append('H')
            zipped_list_resnum_new.append(zipped_list_new[C_C_x[p]][4])
            
        if C_C_y[p] not in C_C_new:
            
            zipped_list_pos_new.append(zipped_list_new[C_C_y[p]][3])
            zipped_list_sym_new.append('H') 
            zipped_list_resname_new.append(zipped_list_new[C_C_y[p]][2])
            zipped_list_atomtype_new.append('H')
            zipped_list_resnum_new.append(zipped_list_new[C_C_y[p]][4])
 

#Adding Lns to list
zipped_list_pos_new.append(zipped_list_new[-1][3])
zipped_list_sym_new.append(zipped_list_new[-1][0])
zipped_list_resname_new.append(zipped_list_new[-1][2])
zipped_list_atomtype_new.append(zipped_list_new[-1][1])
zipped_list_resnum_new.append(zipped_list_new[-1][4])

zipped_list_pos_new = zipped_list_pos_new + zl_pos
zipped_list_sym_new = zipped_list_sym_new + zl_sym
zipped_list_resname_new = zipped_list_resname_new + zl_resname
zipped_list_atomtype_new = zipped_list_atomtype_new + zl_atomtype
zipped_list_resnum_new = zipped_list_resnum_new + zl_resnum


resname_array_new = np.array(zipped_list_resname_new)
atomtype_array_new = np.array(zipped_list_atomtype_new)
resnum_array_new = np.array(zipped_list_resnum_new)

LanM1_Ln_pdb_trunc = Atoms(zipped_list_sym_new, zipped_list_pos_new)


#Print list to pdb file          
#LanM1_Ln_pdb_trunc = Atoms(zipped_list_sym_new, zipped_list_pos_new)


ana_trunc = Analysis(LanM1_Ln_pdb_trunc)
C_H_bonds_new = ana_trunc.get_bonds('C', 'H', unique=True)
C_H_length_new = ana_trunc.get_values(C_H_bonds_new)

C_H_x_new = []
C_H_y_new = []
    
for i6 in range((len(C_H_bonds_new[0]))):
    
    C_H_x_new.append(C_H_bonds_new[0][i6][0]) #all the Cs connected to Hs
    C_H_y_new.append(C_H_bonds_new[0][i6][1]) #all the Hs

print(zipped_list_pos_new[[C_H_x_new[0]][0]])
print(zipped_list_pos_new[[C_H_y_new[0]][0]])

for r in range(len(C_H_x_new)):
    
    if math.dist(zipped_list_pos_new[[C_H_x_new[r]][0]], zipped_list_pos_new[[C_H_y_new[r]][0]]) > C_H_length_max:
        
        LanM1_Ln_pdb_trunc.set_distance([C_H_x_new[r]][0], [C_H_y_new[r]][0], C_H_length_avg, fix=0) 
        
ana_trunc_2 = Analysis(LanM1_Ln_pdb_trunc)
C_H_bonds_new_2 = ana_trunc.get_bonds('C', 'H', unique=True)

C_H_x_new_2 = []
C_H_y_new_2 = []
    
for i7 in range((len(C_H_bonds_new_2[0]))):
    
    C_H_x_new_2.append(C_H_bonds_new[0][i7][0]) #all the Cs connected to Hs
    C_H_y_new_2.append(C_H_bonds_new[0][i7][1]) #all the Hs
    
for r in range(len(C_H_x_new_2)):
    
    if math.dist(zipped_list_pos_new[[C_H_x_new[r]][0]], zipped_list_pos_new[[C_H_y_new[r]][0]]) > C_H_length_max:
        
        LanM1_Ln_pdb_trunc.set_distance([C_H_x_new[r]][0], [C_H_y_new[r]][0], C_H_length_avg, fix=0)  
        
#Calculating formal charge        
fc = 0

#Calculating Ln charge
fc = fc + 3

#Calculating O-C-O charge of -1
C_O_O = []
for q in range((len(C_O_new))):
    
    if (C_O_new.count(C_O_new[q]) == 2):
        
        if C_O_new[q] not in C_O_O:
        
            C_O_O.append(C_O_new[q])
            fc = fc - 1

#Calculating N-H3 charge             
N_H3 = []
print(len(N_H_new))
for t in range((len(N_H_new))):
    
    if (N_H_new.count(N_H_new[t]) == 3):
        
        if N_H_new[t] not in N_H3:
            
            N_H3.append(N_H_new[t])
            fc = fc + 1
            
#Create list of variables to freeze
freeze = []
for u in range((len(zipped_list_sym_new))):
    
    if math.dist(zipped_list_pos_new[u], zipped_list_new[-1][3]) >= 4.7:
        
        freeze.append(u+1)

#Write the list into pdb file
io.write("md{}_LanM1_{}.pdb".format(snapshot, Ln), LanM1_Ln_pdb_trunc, "proteindatabank")

#Add residue name, type, and number to the pdb file 
LanM1_Ln_pdb_trunc.set_array('residuenames', resname_array_new)
LanM1_Ln_pdb_trunc.set_array('atomtypes', atomtype_array_new)
LanM1_Ln_pdb_trunc.set_array('residuenumbers', resnum_array_new)

#Sorted into pdb order
LanM1_Ln_trunc_sort = sort(LanM1_Ln_pdb_trunc, tags=resnum_array_new)
gromacs.write_gromacs("snapshot_truncated_ns.gro", LanM1_Ln_trunc_sort)                                           