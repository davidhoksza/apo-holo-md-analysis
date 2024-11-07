from __future__ import annotations
from typing import Callable

import numpy as np
from Bio.PDB.Polypeptide import three_to_one

class Atom:
    def __init__(self, index = None, name = None, chain = None, resn = None, resi = None, coord = None) -> None:        
        self.index = int(index)
        self.name = name
        self.chain = chain
        self.resn = resn
        self.resi = int(resi)
        self.coord = np.array(coord)
    
    def __repr__(self) -> str:
        return f'index: {self.index}, name:{self.name}, chain:{self.chain}, resn:{self.resn}, resi:{self.resi}, coord:{self.coord}'
    
    def __str__(self) -> str:
        return self.__repr__(self)
    
    def __lt__(self, other):
        return self.resi < other.resi or (self.resi == other.resi and self.name < other.name)


class Residue:
    def __init__(self, chain = None, resn = None, resi = None) -> None:
        self.chain = chain
        self.resn = resn
        self.resi = resi

    def __repr__(self) -> str:
        return f'chain: {self.chain}, resn:{self.resn}, resi:{self.resi}'
    
    def __str__(self) -> str:
        return self.__repr__(self)

    @classmethod
    def from_atom(cls, atom:Atom):
        return cls(atom.chain, atom.resn, atom.resi)


def get_coords(atoms:list[Atom]) -> np.array[np.array]:
    return np.array([a.coord for a in atoms])


class Structure:
    def __init__(self, model, definition:str, state:int=1) -> None:
                
        self.__model = model
        self.__definition = definition
        self.__state = state

        self.__atoms:list[Atom] = []
        self.__atoms_sorted = False        
        self.__get_atoms_from_model()

        self.__sasa = None

    @property
    def atoms(self):
        return self.__atoms
    
    @property
    def atoms_sorted(self):
        return self.__atoms_sorted
    

    def __get_atoms_from_model(self):
        self.__atoms = [Atom(index=a.index, name=a.name, chain=a.chain, resn=a.resn, resi=a.resi, coord=a.coord) for a in self.__model.atom]
        self.__atoms_sorted = False
    

    def get_residues(self, sort=True) -> list[Residue]:
        sa:list[Atom] = sorted(self.__atoms) if sort else self.__atoms
        last_resi = None
        residues = []
        for a in sa:
            if a.resi != last_resi:
                residues.append(Residue.from_atom(a))
                last_resi = a.resi
        return residues
    

    def get_residues_atoms(self, sort=True) -> list[list[Atom]]:        
        sa:list[Atom] = sorted(self.__atoms) if sort else self.__atoms        
        residues = []
        residue = []
        for a in sa:
            if len(residue) > 0 and a.resi != residue[-1].resi:
                residues.append(residue)
                residue=[]
            residue.append(a)
        residues.append(residue)
                
        return residues
    
    def get_sequence(self, one_letter=True) -> str:
        if one_letter:
            return ''.join([three_to_one(r.resn) for r in self.get_residues()])
        else:
            return '-'.join([r.resn for r in self.get_residues()])

    
    def sort_atoms(self):
        if not self.__atoms_sorted:
            self.__atoms = sorted(self.atoms)
            self.__atoms_sorted = True

    def get_coords(self) -> np.ndarray:
        return get_coords(self.atoms)
    

    def get_distances(self, other:Structure, sort=True) -> np.ndarray:
        
        if sort:
            self.sort_atoms()
            other.sort_atoms()

        return np.sqrt(np.sum((self.get_coords() - other.get_coords())**2, axis=1))
    
    
    def rmsd_fast(self, other:Structure):
        assert(self.get_sequence() == other.get_sequence())

        a1 = np.array([a.coord for a in self.atoms])
        a2 = np.array([a.coord for a in other.atoms])

        assert(len(a1) == len(a2))

        return np.sqrt((((a1 - a2)**2).sum())/len(a1))

    
    def rmsd(self, other:Structure, sort=True, atom_names=None, get_coords_from_residue:Callable[[list[Atom]], np.ndarray]=None) -> float:
        
        assert(self.get_sequence() == other.get_sequence())

        if sort:
            self.sort_atoms()
            other.sort_atoms()        
        
        if get_coords_from_residue:
            coords1 = np.array([get_coords_from_residue(residue_atoms) for residue_atoms in self.get_residues_atoms()]).flatten()
            coords2 = np.array([get_coords_from_residue(residue_atoms) for residue_atoms in other.get_residues_atoms()]).flatten()

        else:            
            atoms1 = [a for a in atoms1 if a.name in atom_names] if atom_names else self.atoms
            atoms2 = [a for a in atoms2 if a.name in atom_names] if atom_names else other.atoms

            coords1 = get_coords(atoms1)
            coords2 = get_coords(atoms2)

        assert(len(coords1) == len(coords2))


        return np.sqrt((((coords1 - coords2)**2).sum())/len(coords1))
    
    def get_sasa(self, cmd):
        if self.__sasa is None:
            self.__sasa = cmd.get_area(self.__definition, self.__state)

        return self.__sasa
        
    


        
    