import json

from chemfeagenerator import ChemFeaGenerator
from pubchemhelper import PubChemHelper


class DataParser:
    def __init__(self, content_file_path, compounds_file_path):
        """

        :rtype: object
        """
        self.compounds_file_path = compounds_file_path
        self.content_file_path = content_file_path

    def parse_food_db_compounds_json_to_dict(self) -> object:
        """
        :rtype: object
        """
        compounds = {}
        with open(self.compounds_file_path) as file:
            for line in file:
                line = json.loads(line)
                key = line["id"]
                if key not in compounds:
                    compounds[key] = line
                else:
                    print(
                        "{0} seems a duplicate compound id in food db json file or it is missing from the file".format(
                            str(key)))

        return compounds

    def parse_food_db_content_json_to_dict(self, compounds) -> object:
        """
        :rtype: object
        """
        data = []
        with open(self.content_file_path) as file:
            for line in file:
                line = json.loads(line)
                source_type = line["source_type"]
                if (source_type == "Compound"):
                    food_id = "FOOD" + ("0" * (5 - len(str(line["food_id"])))) + str(line["food_id"])
                    orig_food_id = line["orig_food_id"]
                    compound_id = line["source_id"]
                    orig_content = line["orig_content"]
                    orig_unit = line["orig_unit"]

                    # check if the compound id is in the compounds json file
                    if compound_id in compounds:
                        compound_details = compounds[compound_id]
                        cas_number = compound_details["cas_number"]
                        cid = PubChemHelper.cas_to_pubchem(cas_number)
                        print(food_id)
                        print(orig_content)
                        print(compound_details)
                        print(compound_details["moldb_smiles"])
                        if cid != -1 and cid != -2:
                            print(cid)
                            smiles = PubChemHelper.cid_to_smiles(cid)
                            desc = ChemFeaGenerator.get_chem_descriptors(smiles)
                            pubchemfp = ChemFeaGenerator.get_pubchem_fingerprint(cid)
                            morgans = ChemFeaGenerator.get_morgan_fingerprint(smiles)
                            features = ChemFeaGenerator.generate_mol2vec_features(smiles)
                            print(features)
                            print(len(features))
                        exit(1)

        return data
