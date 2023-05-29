from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import deepchem as dc
import requests
import json
from base64 import b64decode


class ChemFeaGenerator:
    def __init__(self, smiles, cid):
        self.smiles = smiles
        self.cid = cid

    def get_chem_descriptors(self):
        # Create RDKit molecule object
        molecule = Chem.MolFromSmiles(self.smiles)

        # Get list of all descriptor functions in RDKit
        all_descriptors = [x[0] for x in Descriptors._descList]

        # Calculate all descriptors
        desc = []
        for descriptor in all_descriptors:
            descriptor_function = getattr(Descriptors, descriptor)
            desc.append(descriptor_function(molecule))
        return desc

    def get_morgan_fingerprint(self):
        # Create RDKit molecule object
        molecule = Chem.MolFromSmiles(self.smiles)

        # Generate Morgan fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(molecule, 2, nBits=1024)
        return [int(x) for x in fp.ToBitString()]

    def get_pubchem_fingerprint(self):
        cid = self.cid
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Fingerprint2D/JSON"
        response = requests.get(pubchem_url)
        data = json.loads(response.text)
        pcfp_base64 = data['PropertyTable']['Properties'][0]['Fingerprint2D']
        # Convert hexadecimal to binary
        pcfp_bitstring = "".join(["{:08b}".format(x) for x in b64decode(pcfp_base64)])[32:913]
        pcfp_bitstring = [char for char in pcfp_bitstring]

        return pcfp_bitstring

    def generate_mol2vec_features(self):
        featurizer = dc.feat.Mol2VecFingerprint()
        features = featurizer.featurize(self.smiles)
        features = features.tolist()
        features = features[0]
        return features

    def append_list(self, list1, list2):
        for i in list2:
            list1.append(str(i))
        return list1

    def get_all_features(self):
        cid = self.cid
        smiles = self.smiles
        features = []
        desc = self.get_chem_descriptors()
        morgans = self.get_morgan_fingerprint()
        pubchemfp = self.get_pubchem_fingerprint()
        mol2vecfea = self.generate_mol2vec_features()
        self.append_list(features, desc)
        self.append_list(features, morgans)
        self.append_list(features, pubchemfp)
        self.append_list(features, mol2vecfea)
        return features

    def get_gcn_features(self):
        # Convert the SMILES string to a RDKit molecule
        mol = dc.feat.graph_features.ConvMolFeaturizer().featurize([self.smiles])
        features = mol[0].atom_features
        return features
