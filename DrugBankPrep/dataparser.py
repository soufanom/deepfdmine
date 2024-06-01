from lxml import etree

class Drug:
    def __init__(self, id, name, drug_type, drug_group, description, smiles, food_interactions, drug_interactions, description_interactions, pubchem_id, kegg_id):
        self.id = id
        self.name = name
        self.drug_type = drug_type
        self.drug_group = drug_group
        self.description = description
        self.smiles = smiles
        self.food_interactions = food_interactions
        self.drug_interactions = drug_interactions
        self.description_interactions = description_interactions
        self.pubchem_id = pubchem_id
        self.kegg_id = kegg_id

def safe_str(obj):
    """Convert None to an empty string"""
    return '' if obj is None else str(obj)

def parse_cancer_drugbank_xml(file_path):
    """
    Parses the DrugBank XML file using lxml and extracts essential information about each drug.
    """
    tree = etree.parse(file_path)
    root = tree.getroot()

    # Namespace used in DrugBank XML
    ns = {'db': 'http://www.drugbank.ca'}

    drugs = []

    for drug_elem in root.xpath('//db:drug', namespaces=ns):
        drug_id_elem = drug_elem.find('db:drugbank-id', namespaces=ns)
        name_elem = drug_elem.find('db:name', namespaces=ns)
        description_elem = drug_elem.find('db:description', namespaces=ns)
        smiles = drug_elem.find('.//db:calculated-properties/db:property[db:kind="SMILES"]/db:value', ns)
        if smiles is not None:
            smiles = smiles.text

        drug_type = drug_elem.get('type')  # Attribute of the drug element
        drug_group = drug_elem.find('db:groups/db:group', ns).text if drug_elem.find('db:groups/db:group',
                                                                           ns) is not None else None

        # Extracting food interactions
        food_interactions = [interaction.text for interaction in
                             drug_elem.findall('.//db:food-interactions/db:food-interaction', ns)]

        # Extracting drug interactions
        drug_interactions = [interaction.text for interaction in
                             drug_elem.findall('.//db:drug-interactions/db:drug-interaction/db:drugbank-id', ns)]

        # Extracting drug interactions
        description_interactions = [interaction.text for interaction in
                             drug_elem.findall('.//db:drug-interactions/db:drug-interaction/db:description', ns)]

        drug_id = drug_id_elem.text if drug_id_elem is not None else None
        name = name_elem.text if name_elem is not None else None
        description = description_elem.text if description_elem is not None else None

        # Extracting PubChem Compound ID
        pubchem_id = drug_elem.find(
            './/db:external-identifiers/db:external-identifier[db:resource="PubChem Compound"]/db:identifier', ns)
        pubchem_id = pubchem_id.text if pubchem_id is not None else None

        # Extracting KEGG Compound ID
        kegg_id = drug_elem.find(
            './/db:external-identifiers/db:external-identifier[db:resource="KEGG Compound"]/db:identifier', ns)
        kegg_id = kegg_id.text if kegg_id is not None else None

        cancer_drug = False
        for category in drug_elem.findall(".//db:categories/db:category/db:category", ns):
            if category.text.find("Antineoplastic") != -1:
                cancer_drug = True

        if drug_id and name and cancer_drug:  # Ensure we have the minimum required data
            if drug_type == "small molecule":
                drug = Drug(drug_id, name, drug_type, drug_group, description, smiles, food_interactions, drug_interactions, description_interactions, pubchem_id, kegg_id)
                drugs.append(drug)

    return drugs


def parse_drugbank_xml(file_path):
    """
    Parses the DrugBank XML file using lxml and extracts essential information about each drug.
    """
    tree = etree.parse(file_path)
    root = tree.getroot()

    # Namespace used in DrugBank XML
    ns = {'db': 'http://www.drugbank.ca'}

    drugs = []

    for drug_elem in root.xpath('//db:drug', namespaces=ns):
        drug_id_elem = drug_elem.find('db:drugbank-id', namespaces=ns)
        name_elem = drug_elem.find('db:name', namespaces=ns)
        description_elem = drug_elem.find('db:description', namespaces=ns)
        smiles = drug_elem.find('.//db:calculated-properties/db:property[db:kind="SMILES"]/db:value', ns)
        if smiles is not None:
            smiles = smiles.text

        drug_type = drug_elem.get('type')  # Attribute of the drug element
        drug_group = drug_elem.find('db:groups/db:group', ns).text if drug_elem.find('db:groups/db:group',
                                                                           ns) is not None else None

        # Extracting food interactions
        food_interactions = [interaction.text for interaction in
                             drug_elem.findall('.//db:food-interactions/db:food-interaction', ns)]

        # Extracting drug interactions
        drug_interactions = [interaction.text for interaction in
                             drug_elem.findall('.//db:drug-interactions/db:drug-interaction/db:drugbank-id', ns)]

        # Extracting drug interactions
        description_interactions = [interaction.text for interaction in
                             drug_elem.findall('.//db:drug-interactions/db:drug-interaction/db:description', ns)]

        drug_id = drug_id_elem.text if drug_id_elem is not None else None
        name = name_elem.text if name_elem is not None else None
        description = description_elem.text if description_elem is not None else None

        # Extracting PubChem Compound ID
        pubchem_id = drug_elem.find(
            './/db:external-identifiers/db:external-identifier[db:resource="PubChem Compound"]/db:identifier', ns)
        pubchem_id = pubchem_id.text if pubchem_id is not None else None

        # Extracting KEGG Compound ID
        kegg_id = drug_elem.find(
            './/db:external-identifiers/db:external-identifier[db:resource="KEGG Compound"]/db:identifier', ns)
        kegg_id = kegg_id.text if kegg_id is not None else None

        if drug_id and name:  # Ensure we have the minimum required data
            if drug_type == "small molecule":
                drug = Drug(drug_id, name, drug_type, drug_group, description, smiles, food_interactions, drug_interactions, description_interactions, pubchem_id, kegg_id)
                drugs.append(drug)

    return drugs

def query_drugs(drugs, query_id):
    """
    Queries the list of drugs for a specific DrugBank ID.
    """
    for drug in drugs:
        if drug.id == query_id:
            return drug
    return None

# Example Usage
file_path = '/Users/soufanom/Documents/Research/Databases/DrugBank/full database.xml'  # Replace with your DrugBank XML file path
#drugs = parse_drugbank_xml(file_path)
drugs = parse_cancer_drugbank_xml(file_path)

output_file = open('/Users/soufanom/PycharmProjects/deepfdmine/data/DrugBank_cancer_data.txt', "w")
#
for d in drugs:
    output_file.write(safe_str(d.id).strip()+"\t"+safe_str(d.pubchem_id).strip()+"\t"+safe_str(d.kegg_id).strip()+"\t"+safe_str(d.drug_type).strip()+"\t"+safe_str(d.drug_group).strip()+"\t"+safe_str(d.name).strip()+"\t"+safe_str(d.smiles).strip()+"\t"+safe_str(d.description).replace("\n", "").replace("\r", " ").replace("^M^M", " ").replace("\t", " ")+"\t")
    output_file.write(";".join(d.food_interactions).replace("\n", "").replace("\r", " ").replace("^M^M", " ").replace("\t", " ")+"\t"+";".join(d.drug_interactions).replace("\n", "").replace("\r", " ").replace("^M^M", " ").replace("\t", " ")+"\t"+";".join(d.description_interactions).replace("\n", "").replace("\r", " ").replace("^M^M", " ").replace("\t", " "))
    output_file.write("\n")

output_file.close()
# # Query for a specific drug by its DrugBank ID
# query_id = 'DB01048'  # Example DrugBank ID
# queried_drug = query_drugs(drugs, query_id)
#
# if queried_drug:
#     print(f"Drug ID: {queried_drug.id}, Name: {queried_drug.name}, Description: {queried_drug.description}")
# else:
#     print("Drug not found.")
