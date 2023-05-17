import pubchempy as pcp


class PubChemHelper:
    def cas_to_pubchem(cas_number):
        try:
            compounds = pcp.get_compounds(cas_number, 'name')
            if compounds:
                return compounds[0].cid
            else:
                print("No PubChem ID found for this {0} CAS number", cas_number)
                return -1
        except Exception as e:
            print(str(e))
            return -2

    def cid_to_smiles(cid):
        try:
            compound = pcp.Compound.from_cid(cid)
            return compound.canonical_smiles
        except Exception as e:
            print(str(e))
            return "-1"
