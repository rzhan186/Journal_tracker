
import pandas as pd
from Bio import Entrez

# Set Entrez email
Entrez.email = "rzhan186@gmail.com"

journal="Nature"

def fetch_pubmed_articles(journal, volume=None, issue=None, issue_range=None):
    """
    Fetch research articles from a specific journal issue, or automatically detect the latest issue.
    
    Parameters:
        journal (str): Name of the journal (e.g., "Nature").
        volume (str, optional): Specific volume number (e.g., "589"). If None, fetches the latest.
        issue (str, optional): Specific issue number (e.g., "7841"). If None, fetches the latest.
        issue_range (list of tuples, optional): A list of (volume, issue) pairs for a range of issues.
        
    Returns:
        List of dictionaries containing article details (title, abstract, DOI).
    """

    # If no volume/issue is specified, fetch the latest issue
    if volume is None or issue is None:
        print(f"No volume/issue number specified, fetching the latest issue of {journal}...")
        latest_volume, latest_issue = get_latest_journal_issue(journal)
        if latest_volume and latest_issue:
            volume, issue = latest_volume, latest_issue
        else:
            print(f"Could not determine the latest issue for {journal}.")
            return []

    # If issue_range is provided, loop over multiple issues
    if issue_range:
        all_articles = []
        for vol, iss in issue_range:
            print(f"Fetching articles from {journal}, Volume {vol}, Issue {iss}...")
            articles = fetch_articles_by_issue(journal, vol, iss)
            all_articles.extend(articles)
        return all_articles

    # Fetch articles for a single issue
    print(f"Fetching articles from {journal}, Volume {volume}, Issue {issue}...")
    return fetch_articles_by_issue(journal, volume, issue)


def get_latest_journal_issue(journal):
    """
    Find the latest volume and issue number for a given journal using PubMed.
    
    Parameters:
        journal (str): Name of the journal.
        
    Returns:
        (latest_volume, latest_issue) tuple or (None, None) if not found.
    """
    
    # This searches for any article published in the specified journal.
    query = f'"{journal}"[Journal] AND "journal article"[Publication Type] NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type] OR "review"[Publication Type])'
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1, sort="pub date")
    record = Entrez.read(handle)
    handle.close()
    
    if not record["IdList"]:
        return None, None

    # Fetch details of the latest article
    pmid = record["IdList"][0]
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
    paper_info = Entrez.read(handle)
    handle.close()

    # Extract volume and issue
    article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
    volume = article_data.get("Journal", {}).get("JournalIssue", {}).get("Volume")
    issue = article_data.get("Journal", {}).get("JournalIssue", {}).get("Issue")

    return volume, issue





def fetch_articles_by_issue(journal, volume, issue):
    """
    Fetch all research articles from a specific journal issue in PubMed.
    
    Parameters:
        journal (str): Name of the journal.
        volume (str): Volume number.
        issue (str): Issue number.
        
    Returns:
        List of dictionaries containing article details (title, abstract, DOI).
    """

    # Construct PubMed query to filter by journal, volume, and issue
    query = (f'"{journal}"[Journal] AND "{volume}"[Volume] AND "{issue}"[Issue] '
             f'AND "journal article"[Publication Type] '
             f'NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type] OR "review"[Publication Type])')

    # Fetch article IDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="pub date")
    record = Entrez.read(handle)
    handle.close()
    
    pmid_list = record["IdList"]
    papers = []

    if not pmid_list:
        print(f"No articles found for {journal}, Volume {volume}, Issue {issue}.")
        return papers

    # Fetch details for each article
    for pmid in pmid_list:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
        paper_info = Entrez.read(handle)
        handle.close()

        # Extract article details
        article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
        title = article_data["ArticleTitle"]
        abstract = article_data.get("Abstract", {}).get("AbstractText", ["No abstract available"])[0]

        # Extract DOI if available
        doi = next((id_item for id_item in paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
                    if id_item.attributes["IdType"] == "doi"), None)

        papers.append({
            "Journal": journal,
            "Volume": volume,
            "Issue": issue,
            "Title": title,
            "Abstract": abstract,
            "DOI": f"https://doi.org/{doi}" if doi else "No DOI available"
        })

    return papers

# Example Usage

# Fetch the latest issue automatically
articles_latest = fetch_pubmed_articles("Nature")

# Fetch a specific issue manually
articles_specific = fetch_pubmed_articles("Nature", volume="589", issue="7841")

# Fetch a range of issues
issue_range = [("589", "7841"), ("589", "7842"), ("590", "7843")]
articles_range = fetch_pubmed_articles("Nature", issue_range=issue_range)

# Save results to CSV
df = pd.DataFrame(articles_latest)
df.to_csv("journal_articles.csv", index=False)

print(f"Fetched {len(df)} articles and saved to journal_articles.csv")
