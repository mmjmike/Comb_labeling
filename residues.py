class Residue:

    def __init__(self, name, one_letter_code, three_letters_code,
                 mol_weight, pack_weight, brutto, proline=False):
        self.name = name
        self.one_let = one_letter_code
        self.three_let = three_letters_code
        self.MW = mol_weight
        self.PW = pack_weight
        self.usage = 1
        self.brutto = brutto
        self.proline = proline


class ResidueTypes:

    def __init__(self):
        self.residues = []
        self.one_let = {}
        self.three_let = {}

    def add_residue(self, residue):
        self.residues.append(residue)
        self.one_let[residue.one_let] = residue
        self.three_let[residue.three_let] = residue

def make_rt():
    rt = ResidueTypes()
    rt.add_residue(Residue("Alanine", "A", "ALA", 89.09, 89.09, "C3H7NO2"))
    rt.add_residue(Residue("Arginine", "R", "ARG", 174.20, 174.20+35.45+1.01, "C6H14N4O2"))
    rt.add_residue(Residue("Asparagine", "N", "ASN", 132.12, 132.12+18, "C4H8N2O3"))
    rt.add_residue(Residue("Aspartic acid", "D", "ASP", 133.10, 133.10, "C4H7NO4"))
    rt.add_residue(Residue("Glutamic acid", "E", "GLU", 147.13, 147.13, "C5H9NO4"))
    rt.add_residue(Residue("Glutamine", "Q", "GLN", 145.15, 145.15, "C5H10N2O3"))
    rt.add_residue(Residue("Glycine", "G", "GLY", 75.07, 75.07, "C2H5NO2"))
    rt.add_residue(Residue("Histidine", "H", "HIS", 155.16, 155.16+35.45+1.01+18, "C6H9N3O2"))
    rt.add_residue(Residue("Isoleucine", "I", "ILE", 131.17, 131.17, "C6H13NO2"))
    rt.add_residue(Residue("Leucine", "L", "LEU", 131.17, 131.17, "C6H13NO2"))
    rt.add_residue(Residue("Lysine", "K", "LYS", 145.19, 145.19+2*35.45+2.02, "C6H14N2O2"))
    rt.add_residue(Residue("Methionine", "M", "MET", 149.21, 149.21, "C5H11NO2S"))
    rt.add_residue(Residue("Phenylalanine", "F", "PHE", 165.19, 165.19, "C9H11NO2"))
    rt.add_residue(Residue("Proline", "P", "PRO", 113.15, 113.15, "C5H9NO2", proline=True))
    rt.add_residue(Residue("Serine", "S", "SER", 105.09, 105.09, "C3H7NO3"))
    rt.add_residue(Residue("Threonine", "T", "THR", 119.12, 119.12, "C4H9NO3"))
    rt.add_residue(Residue("Tryptophan", "W", "TRP", 204.23, 204.23, "C11H12N2O2"))
    rt.add_residue(Residue("Tyrosine", "Y", "TYR", 181.19, 181.19, "C9H11NO3"))
    rt.add_residue(Residue("Valine", "V", "VAL", 117.15, 117.15, "C5H11NO2"))
    rt.add_residue(Residue("Cysteine", "C", "CYS", 121.16, 121.16, "C3H7NO2S"))
    return rt