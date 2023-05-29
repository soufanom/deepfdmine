# write main function to call the class
#
import sys

from DFDMinePrep.dataparser import DataParser


def prepare_food_db_compounds():
    var_1 = sys.argv[1]
    path = "/Users/soufanom/Documents/Research/Databases/FooDB/foodb_2020_04_07_json/"
    compounds_file = path + var_1
    output_file = path + "/output/" + var_1 + "_out.txt"
    # create an object of DataParser class
    content_file = path + '/Content_compound_sample' + '.json'
    d1 = DataParser(content_file, compounds_file, output_file)
    # call parse_json_to_dict method
    compounds = d1.parse_food_db_compounds_json_to_dict()
    d1.write_compounds_fea_to_file(compounds)


if __name__ == '__main__':
    # Uncomment the following code to reproduce generating chemical features for compounds in FooDB
    # prepare_food_db_data()
    # content = d1.parse_food_db_content_json_to_dict(compounds)
    x = 10
    # TO-DO: Check for healthy food compounds
