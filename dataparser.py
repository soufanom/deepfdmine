# write a class to parse the data from json file to dictionary
# and write a class to parse the data from dictionary to json file
import json


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
        with open(self.content_file_path) as file:
            for line in file:
                line = json.loads(line)
                key = line["source_id"]
                if key not in compounds:
                    compounds[key] = line
                else:
                    print("check why there are duplicate compounds ids in food db json file")

        return compounds

    def parse_food_db_content_json_to_dict(self, compounds) -> object:
        """
        :rtype: object
        """
        # read compounds json file
        compounds_data = json.load(open(self.compounds_file_path))
        data = []
        with open(self.content_file_path) as file:
            for line in file:
                line = json.loads(line)
                source_type = line["source_type"]
                if(source_type == "Compound"):
                    food_id = "FOOD"+("0"*(5-len(str(line["food_id"]))))+str(line["food_id"])
                    orig_food_id = line["orig_food_id"]
                    compound_id = line["source_id"]
                    orig_content = line["orig_content"]
                    orig_unit = line["orig_unit"]
                    print(food_id)

                    # check if the compound id is in the compounds json file
                    if compound_id in compounds:
                        compound_details = compounds[compound_id]
                        print(food_id)
                        print(orig_content)
                        print(compound_details)
                        exit(1)

        return data
