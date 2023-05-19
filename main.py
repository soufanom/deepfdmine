# write main function to call the class
#
from dataparser import DataParser

if __name__ == '__main__':
    # create an object of DataParser class
    content_file = '/Users/soufanom/Documents/Research/Databases/FooDB/foodb_2020_04_07_json/Content_compound_sample' \
                  '.json'
    compounds_file = '/Users/soufanom/Documents/Research/Databases/FooDB/foodb_2020_04_07_json/Compound.json'
    d1 = DataParser(content_file, compounds_file)
    # call parse_json_to_dict method
    compounds = d1.parse_food_db_compounds_json_to_dict()
    d1.write_compounds_fea_to_file(compounds)
    #content = d1.parse_food_db_content_json_to_dict(compounds)




