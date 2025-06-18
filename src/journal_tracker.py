# first creating a python virtual env for this task
# shift+command+p --> python environment --> venv

# activate the virtual environment by select the correct python interpretor
# .venv

# run pip install biopython
# run pip install pandas

import pandas as pd
from Bio import Entrez

# 设置Entrez邮箱（PubMed API要求提供邮箱）
Entrez.email = "rzhan186@gmail.com"

# 目标期刊列表
# journals = [
#     "Nature", "Science", "Cell", 
#     "Nature Microbiology", "ISME Journal", "Environmental Microbiology",
#     "PNAS", "mSystems", "Microbiome", "Genome Biology"
# ]

journals = [
    "Nature"
]

# 生命科学相关关键词
keywords = [""]

# Function to fetch research articles only
# journal: A string representing the name of the journal to query (e.g., "Nature").
# max_results=20: The maximum number of results to return (default is 20).

def fetch_pubmed_research_papers(journal, max_results=20):

    """Fetch only research articles from a specific journal in PubMed"""
    
    # Filter for journal research articles only, excluding editorials, news, etc.
    # Excludes non-research content:
        # "news"[Publication Type]
        # "comment"[Publication Type]
        # "editorial"[Publication Type]
        # "review"[Publication Type]
    query = f'"{journal}"[Journal] AND "journal article"[Publication Type] NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type] OR "review"[Publication Type])'

    # Fetch article IDs
    # Performs a PubMed search using Entrez and retrieves article IDs, sorted by publication date.
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="pub date")

    # Entrez.read returns a dictionary containing various fields
    record = Entrez.read(handle)
    handle.close()
    
    pmid_list = record["IdList"]
    papers = []
    

    # Loops through each PMID in the list.
    for pmid in pmid_list:
        # Uses the NCBI Entrez API to fetch detailed information about a PubMed article based on its PubMed ID (PMID). Let's break it down:
        # Entrez.efetch	This function retrieves full records from an NCBI database.
        # db="pubmed"	Specifies that we are fetching data from PubMed (other options include db="nucleotide", db="protein", etc.).
        # id=pmid	Specifies the PubMed ID (PMID) of the article to retrieve. This can be a single ID (e.g., "38274581") or multiple IDs (e.g., "38274581,38265879").
        # rettype="xml"	Specifies that the article data should be returned in XML format, which is structured and machine-readable.
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
        paper_info = Entrez.read(handle)
        handle.close()
        
        # Parse article information
        # Extracts the main article data from the fetched XML response.

        article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
        title = article_data["ArticleTitle"]
        abstract = article_data.get("Abstract", {}).get("AbstractText", ["No abstract available"])[0]
        
        # Extract DOI if available
        # why not doi = []? Answer: A DOI (Digital Object Identifier) for a research article is typically a single string (e.g., "10.1038/s41586-024-05789"). Since we only need one DOI per paper, using None makes it explicit that there is either one DOI or none.
        doi = None 
        for id_item in paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]:
            if id_item.attributes["IdType"] == "doi":
                doi = id_item
                break

        papers.append({
            "Journal": journal,
            "Title": title,
            "Abstract": abstract,
            "DOI": f"https://doi.org/{doi}" if doi else "No DOI available"
        })
    
    return papers

# Filter papers using keywords
def filter_papers(papers, keywords):
    """Filter papers based on specific keywords in the title"""
    return [p for p in papers if any(kw.lower() in p["Title"].lower() for kw in keywords)]


# run the search
all_papers = []
for journal in journals:
    print(f"Fetching papers from {journal}...")
    papers = fetch_pubmed_research_papers(journal, max_results=20)
    filtered_papers = filter_papers(papers, keywords)
    all_papers.extend(filtered_papers)


# save as a markdown file
output_file = "latest_life_science_papers.md"
with open(output_file, "w", encoding="utf-8") as f:
    f.write("# Latest research article\n\n")
    for paper in all_papers:
        f.write(f"## {paper['Title']}\n")
        f.write(f"**Journal**: {paper['Journal']}\n\n")
        f.write(f"**Abstract**: {paper['Abstract']}\n\n")
        f.write(f"🔗 [DOI-link]({paper['DOI']})\n\n")
        f.write("---\n\n")

print(f"The list of articles is saved as {output_file}")



