import requests
import json

def get_pdb_ids(uniprot_id):
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        }
    }
    
    response = requests.post(
        search_url, 
        data=json.dumps(query),
        headers={"Content-Type": "application/json"}
    )
    
    if response.status_code == 200:
        results = response.json()
        pdb_ids = [hit['identifier'] for hit in results.get('result_set', [])]
        return pdb_ids
    elif response.status_code == 204:
        return []
    else:
        response.raise_for_status()

if __name__ == "__main__":
    uniprot_id = "O75469"
    pdb_ids = get_pdb_ids(uniprot_id)
    print(f"PDB IDs for UniProt {uniprot_id}:")
    for pdb_id in pdb_ids:
        print(pdb_id)
