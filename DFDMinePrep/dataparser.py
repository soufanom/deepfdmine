import json

from chemfeagenerator import ChemFeaGenerator
from pubchemhelper import PubChemHelper


class DataParser:
    def __init__(self, content_file_path, compounds_file_path, output_file_path):
        """

        :rtype: object
        """
        self.compounds_file_path = compounds_file_path
        self.content_file_path = content_file_path
        self.output_file_path = output_file_path

    def parse_food_db_compounds_json_to_dict(self) -> object:
        """
            Map the compounds json file to a dictionary (key is compound id and value is the compound details)
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

    def write_compounds_fea_to_file(self, compounds):
        outfile = open(self.output_file_path, "a")
        for index, key in enumerate(compounds):
            compound_id = key
            compound_details = compounds[compound_id]
            # For the following line of code, we are assuming that the name of the compound is the CAS number
            #   We can also use the name when cas number is not available
            cas_number = compound_details["name"]
            public_id = compound_details["public_id"]
            smiles_moldb = compound_details["moldb_smiles"]
            try:
                cid = PubChemHelper.cas_to_pubchem(cas_number, smiles_moldb)
                if cid is not None:
                    outfile.write(str(public_id) + "\t" + str(cid) + "\t")
                    smiles = PubChemHelper.cid_to_smiles(cid)
                    outfile.write(str(smiles) + "\t")
                    compound_obj = ChemFeaGenerator(smiles, cid)
                    features = compound_obj.get_all_features()
                    features = ",".join(features)
                    outfile.write(features + "\n")
            except Exception as e:
                print("Exception occurred for compound id: " + str(public_id) + " " + str(e))
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
