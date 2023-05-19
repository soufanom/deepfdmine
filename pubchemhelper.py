import pubchempy as pcp


class PubChemHelper:
    @staticmethod
    def cas_to_pubchem(cas_number, smiles):
        if cas_number is not None:
            compounds = pcp.get_compounds(cas_number, 'name')
            if compounds:
                return compounds[0].cid
            elif smiles is not None:
                # print("No PubChem ID found for this {0} CAS number", cas_number, "... Searching for the SMILES "
                #                                                                 "instead.")
                # Search for compounds using the provided SMILES
                query_cid = PubChemHelper.get_most_similar_compound(smiles)
                return query_cid
            else:
                return None
        elif smiles is not None:
            # Search for compounds using the provided SMILES
            query_cid = PubChemHelper.get_most_similar_compound(smiles)
            return query_cid
        else:
            print("No CAS number or SMILES provided.")
            return None

    @staticmethod
    def get_most_similar_compound(smiles):
        # Search for compounds using the provided SMILES
        results = pcp.get_compounds(smiles, 'smiles')

        if not results:
            print("No compounds found for the given SMILES.")
            return None

        # Get the first compound from the search results
        query_compound = results[0]
        if query_compound:
            # Get the CID (Compound Identifier) of the query compound
            query_cid = query_compound.cid
            return query_cid
        else:
            return None

    @staticmethod
    def cid_to_smiles(cid):
        try:
            compound = pcp.Compound.from_cid(cid)
            return compound.canonical_smiles
        except Exception as e:
            print(str(e))
            return "-1"
