from Bio import Entrez
import re

Entrez.email = "rzhan186@gmail.com"

def expand_keywords_with_mesh(term):
    """
    Expands a keyword using MeSH via Entrez and returns synonyms/related terms.
    Handles XML parsing properly.
    """
    try:
        # Search for the MeSH term ID
        search_handle = Entrez.esearch(db="mesh", term=term)
        search_result = Entrez.read(search_handle)
        search_handle.close()

        if not search_result["IdList"]:
            print(f"⚠️ No MeSH term found for '{term}'")
            return [term]

        mesh_id = search_result["IdList"][0]

        # Fetch the MeSH XML data
        fetch_handle = Entrez.efetch(db="mesh", id=mesh_id, rettype="xml", retmode="text")
        
        # Read the response
        fetch_data = fetch_handle.read()  # Read the response data
        fetch_handle.close()

        # We are assuming this is a valid XML string we can work with
        if "<mesh-terms>" not in fetch_data:  # Simple check to avoid parsing bad data
            print(f"⚠️ Unexpected data format for term '{term}': {fetch_data[:100]}...")  # Print first 100 chars
            return [term]

        # Parse synonyms from the data
        start = fetch_data.find("<Concept>")
        end = fetch_data.find("</Concept>") + len("</Concept>")

        if start == -1 or end == -1:
            print(f"⚠️ Couldn't locate concept tags in data for '{term}'")
            return [term]

        concept_data = fetch_data[start:end]
        terms = re.findall(r'<Term>(.*?)</Term>', concept_data)
        if not terms:
            print(f"⚠️ No terms found in the concept data for '{term}'")
            return [term]

        return [t for t in terms]

    except Exception as e:
        print(f"❌ Error expanding MeSH term '{term}': {e}")
        return [term]

STOPWORDS = {"the", "and", "of", "to", "in", "on", "with", "a", "an", "by", "for", "from", "at"}

def generate_keywords_query(user_input):
    """
    Converts user input into a PubMed-compatible search query.
    Expands terms using MeSH and targets both title and abstract fields.
    """
    # Extract words and remove stopwords
    raw_words = re.findall(r'\b\w+\b', user_input.lower())
    base_keywords = [word for word in raw_words if word not in STOPWORDS]

    if not base_keywords:
        raise ValueError("No valid keywords detected after filtering.")

    # Expand each word with MeSH synonyms
    all_keywords = set()
    for kw in base_keywords:
        expanded = expand_keywords_with_mesh(kw)
        all_keywords.update(expanded)

    # Format keywords for PubMed
    query_parts = [f'"{kw}"[Title/Abstract]' for kw in all_keywords]
    final_query = " OR ".join(query_parts)

    print(f"🔍 Final expanded query:\n{final_query}")
    return final_query

# Prompt user input
topic = input("Enter a topic: ")
keywords_query = generate_keywords_query(topic)

print(keywords_query)

# Dummy commented line for future use
# articles = fetch_pubmed_articles_by_date(formatted_journal, start_date, end_date, keywords_query)