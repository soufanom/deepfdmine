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
    content_output_path = "data/"
    d1 = DataParser(content_file, compounds_file, output_file, content_output_path)
    # call parse_json_to_dict method
    compounds = d1.parse_food_db_compounds_json_to_dict()
    d1.write_compounds_fea_to_file(compounds)


if __name__ == '__main__':
    var_1 = "Compound.json"#sys.argv[1]
    path = "/Users/soufanom/Documents/Research/Databases/FooDB/foodb_2020_04_07_json/"
    compounds_file = path + var_1
    output_file = path + "/output2/" + var_1 + "_out.txt"
    # create an object of DataParser class
    content_file = path + '/Content_compound' + '.json'
    content_output_path = "data2/"
    d1 = DataParser(content_file, compounds_file, output_file, content_output_path)

    compounds = d1.parse_food_db_compounds_json_to_dict()
    compounds_pubchem_id = DataParser.parse_food_db_pubchem_file("data/Foodb_final.txt")
    d1.parse_food_db_content_json_to_dict(compounds, compounds_pubchem_id)
    # Uncomment the following code to reproduce generating chemical features for compounds in FooDB
    # prepare_food_db_compounds()
    # content = d1.parse_food_db_content_json_to_dict(compounds)
    #prepared_food_db_content()
    # TO-DO: Check for healthy food compounds
