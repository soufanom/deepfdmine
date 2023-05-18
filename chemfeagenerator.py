from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import deepchem as dc
import requests
import json
from base64 import b64decode

class ChemFeaGenerator:
    def get_chem_descriptors(smiles):
        # Create RDKit molecule object
        molecule = Chem.MolFromSmiles(smiles)

        # Get list of all descriptor functions in RDKit
        all_descriptors = [x[0] for x in Descriptors._descList]

        # Calculate all descriptors
        desc = []
        for descriptor in all_descriptors:
            descriptor_function = getattr(Descriptors, descriptor)
            desc.append(descriptor_function(molecule))
        return desc

    def get_morgan_fingerprint(smiles):
        # Create RDKit molecule object
        molecule = Chem.MolFromSmiles(smiles)

        # Generate Morgan fingerprint
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, 2, nBits=1024)
        return fingerprint


    def get_pubchem_fingerprint(cid):
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Fingerprint2D/JSON"
        response = requests.get(pubchem_url)
        data = json.loads(response.text)
        pcfp_base64 = data['PropertyTable']['Properties'][0]['Fingerprint2D']
        # Convert hexadecimal to binary
        pcfp_bitstring = "".join(["{:08b}".format(x) for x in b64decode(pcfp_base64)])[32:913]
        pcfp_bitstring = [char for char in pcfp_bitstring]

        return pcfp_bitstring

    def generate_mol2vec_features(smiles):
        featurizer = dc.feat.Mol2VecFingerprint()
        features = featurizer.featurize(smiles)
        return features

    def generate_features(smiles, cid):
        features = []
        desc = ChemFeaGenerator.get_chem_descriptors(smiles)
        morgans = ChemFeaGenerator.get_morgan_fingerprint(smiles)
        pubchemfp = ChemFeaGenerator.get_pubchem_fingerprint(cid)
        features.append(desc)
        features.append(morgans)
        features.append(pubchemfp)
        return features