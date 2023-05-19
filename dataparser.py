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

    # write a method to write array elements to a file with comman separated values
    def write_array_to_file(self, array):
        outfile = open("array.txt", "w")
        for element in array:
            outfile.write(str(element) + ",")
        outfile.close()

    def write_compounds_fea_to_file(self, compounds):
        outfile = open("compounds_fea_food_db.txt", "w")
        for compound_id in compounds:
            compound_details = compounds[compound_id]
            cas_number = compound_details["cas_number"]
            public_id = compound_details["public_id"]
            smiles_moldb = compound_details["moldb_smiles"]
            print(cas_number)
            print(public_id)
            print(smiles_moldb)
            cid = PubChemHelper.cas_to_pubchem(cas_number, smiles_moldb)
            if cid is not None:
                outfile.write(str(cid) + "\t" + str(public_id) + "\t")
                smiles = PubChemHelper.cid_to_smiles(cid)
                outfile.write(str(smiles) + "\t")
                compound_obj = ChemFeaGenerator(smiles, cid)
                features = compound_obj.get_all_features()
                features = ",".join(features)
                outfile.write(features + "\n")
        outfile.close()

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
                        print(food_id)
                        print(orig_content)
                        print(compound_details)
                        print(compound_details["moldb_smiles"])


        return data
