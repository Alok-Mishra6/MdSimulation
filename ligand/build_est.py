#!/usr/bin/env python3
"""Build 17-beta-estradiol with hydrogens, in the AlphaFold 3 predicted pose.

AlphaFold returns heavy atoms only, and a force field needs hydrogens. The
chemistry (which atom bonds to which, and with what order) is NOT guessed from
interatomic distances here -- it is read from the PDB Chemical Component
Dictionary entry for EST, which is the authoritative definition. Only the
*coordinates* come from the prediction.

The result is checked against the CCD's own InChIKey, so a wrong stereoisomer
cannot pass silently: 17-beta-estradiol and its 17-alpha epimer differ by one
chiral centre and would give different keys.

Usage:  build_est.py <EST.cif> <holo.pdb> <chain> <out.pdb>
"""
import re
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

REFERENCE_INCHIKEY = "VOXZDWNPVJITMN-ZBRFXRBCSA-N"  # from the CCD entry itself

BOND_ORDER = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "AROM": Chem.BondType.AROMATIC,
}


def _loop(text, prefix):
    """Rows of the first cif loop_ whose columns all start with `prefix`."""
    cols, rows, in_loop = [], [], False
    for line in text.splitlines():
        s = line.strip()
        if s == "loop_":
            cols, in_loop = [], True
            continue
        if in_loop and s.startswith(prefix):
            cols.append(s.split(".", 1)[1])
            continue
        if cols and s and not s.startswith("_"):
            if s.startswith("#"):
                break
            rows.append(s.split())
            continue
        if cols and rows:
            break
        if s.startswith("_"):
            in_loop = False
    return cols, rows


def read_ccd(path):
    """(name -> element, [(name1, name2, rdkit bond type)]) from a CCD entry."""
    text = open(path).read()
    acols, arows = _loop(text, "_chem_comp_atom.")
    ia, ie = acols.index("atom_id"), acols.index("type_symbol")
    elements = {r[ia].strip('"'): r[ie] for r in arows}

    bcols, brows = _loop(text, "_chem_comp_bond.")
    i1 = bcols.index("atom_id_1")
    i2 = bcols.index("atom_id_2")
    io = bcols.index("value_order")
    bonds = [(r[i1].strip('"'), r[i2].strip('"'), BOND_ORDER[r[io]])
             for r in brows]
    return elements, bonds


def read_pose(path, chain):
    """name -> (x, y, z) for the HETATM records of one chain."""
    out = {}
    for line in open(path):
        if not line.startswith("HETATM"):
            continue
        if line[21] != chain:
            continue
        out[line[12:16].strip()] = (float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54]))
    return out


def main():
    ccd_path, pdb_path, chain, out_path = sys.argv[1:5]
    elements, bonds = read_ccd(ccd_path)
    pose = read_pose(pdb_path, chain)

    heavy = [n for n, e in elements.items() if e != "H"]
    missing = sorted(set(heavy) - set(pose))
    extra = sorted(set(pose) - set(heavy))
    if missing or extra:
        sys.exit("atom-name mismatch against the CCD: missing %s, unexpected %s"
                 % (missing, extra))

    # Heavy-atom skeleton: CCD chemistry, predicted coordinates.
    mol = Chem.RWMol()
    index = {}
    for name in heavy:
        index[name] = mol.AddAtom(Chem.Atom(elements[name]))
    n_heavy_bonds = 0
    for a, b, order in bonds:
        if a in index and b in index:
            mol.AddBond(index[a], index[b], order)
            n_heavy_bonds += 1

    conf = Chem.Conformer(mol.GetNumAtoms())
    for name, i in index.items():
        conf.SetAtomPosition(i, Point3D(*pose[name]))
    mol.AddConformer(conf)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    # Stereochemistry is read back off the predicted 3D coordinates, not assumed.
    Chem.AssignStereochemistryFrom3D(mol)

    withh = Chem.AddHs(mol, addCoords=True)

    # Relax hydrogens only: the heavy-atom frame is the prediction and is kept.
    ff_ok = True
    mp = AllChem.MMFFGetMoleculeProperties(withh)
    if mp is None:
        ff_ok = False
    else:
        ff = AllChem.MMFFGetMoleculeForceField(withh, mp)
        for atom in withh.GetAtoms():
            if atom.GetAtomicNum() > 1:
                ff.AddFixedPoint(atom.GetIdx())
        ff.Minimize(maxIts=2000)

    key = Chem.MolToInchiKey(withh)
    print("heavy atoms      : %d" % mol.GetNumAtoms())
    print("heavy bonds      : %d" % n_heavy_bonds)
    print("with hydrogens   : %d atoms" % withh.GetNumAtoms())
    print("formula          : %s" % Chem.rdMolDescriptors.CalcMolFormula(withh))
    print("net charge       : %d" % Chem.GetFormalCharge(withh))
    print("hydrogen relax   : %s" % ("MMFF94, heavy atoms fixed" if ff_ok
                                     else "SKIPPED (no MMFF params)"))
    print("InChIKey         : %s" % key)
    print("CCD reference    : %s" % REFERENCE_INCHIKEY)

    if key != REFERENCE_INCHIKEY:
        if key.split("-")[0] == REFERENCE_INCHIKEY.split("-")[0]:
            sys.exit("STEREOCHEMISTRY MISMATCH: same skeleton, different "
                     "stereoisomer. The predicted pose is not 17-beta-estradiol.")
        sys.exit("CONSTITUTION MISMATCH: this is not estradiol at all.")
    print("stereochemistry  : MATCHES 17-beta-estradiol")

    for atom in withh.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("EST")
        info.SetResidueNumber(1)
        info.SetChainId(chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    Chem.MolToPDBFile(withh, out_path, flavor=4)
    print("wrote            : %s" % out_path)

    # acpype takes MDL or MOL2, never PDB (PDB input would make it shell out to
    # obabel). An MDL molfile also carries the bond orders explicitly, so
    # antechamber never has to re-guess the chemistry from distances.
    mdl_path = out_path.rsplit(".", 1)[0] + ".mdl"
    withh.SetProp("_Name", "EST")
    Chem.MolToMolFile(withh, mdl_path)
    print("wrote            : %s" % mdl_path)


if __name__ == "__main__":
    main()
