# import modules
from datetime import datetime
from Bio import Entrez
from calendar import monthrange

import pandas as pd
import re
import requests
import json

# Set Entrez email
Entrez.email = "rzhan186@gmail.com"

######################################################################
# separate function for validating date inputs and checking the order.
def validate_date_input(date_str):
    """Validates user input for YYYY-MM or YYYY-MM-DD format."""
    pattern = r"^\d{4}-(0[1-9]|1[0-2])(-([0-2][0-9]|3[01]))?$"  # Matches YYYY-MM or YYYY-MM-DD
    return bool(re.match(pattern, date_str))

def validate_dates(start_date, end_date, today):

    today = datetime.today()

    """Validates the start and end date formats and order."""
    if not validate_date_input(start_date) or not validate_date_input(end_date):
        raise ValueError("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15).")

    if start_date > end_date:
        raise ValueError("❌ Start date cannot be later than end date.")

    if start_date > today.strftime('%Y-%m-%d'):
        raise ValueError("❌ Start date cannot be later than the current date")

    if end_date > today.strftime("%Y-%m-%d"):
        raise ValueError("❌ End date cannot be later than the current date")

def get_last_day_of_month(year, month,):
    """Returns the last valid day of a given month."""
    return monthrange(int(year), int(month))[1]  # Correct last day

######################################################################
# Format journal names to allow full or abbreviations

JOURNAL_LIST_URL = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt" 
    # modify this to host this urc on github
JOURNAL_LIST_FILE = "/Users/ruizhang/Desktop/pubmed_journal_abbreviations.json"  
    # modify this to host this on github

# function to get the Pubmed journal title and its abbreviations
def get_pubmed_journal_abbreviations():
    """
    Fetches and parses the latest journal abbreviations from PubMed's official source.
    Saves the data locally in a JSON file for future use.
    Returns a dictionary {abbreviation: full_name}.
    """
    print("🔍 Fetching the latest journal abbreviations from PubMed...")

    response = requests.get(JOURNAL_LIST_URL)
    response.raise_for_status()  # Ensure request was successful
    text_data = response.text.split("\n")  # Split response into lines

    journal_dict = {}
    current_abbr = None
    current_full = None

    for line in text_data:
        line = line.strip()  # Remove accidental spaces or newlines

        if not line or line.startswith("-"):  # Ignore empty lines and separators
            continue

        if line.startswith("MedAbbr: "):
            current_abbr = line.replace("MedAbbr: ", "").strip()
        elif line.startswith("JournalTitle: "):
            current_full = line.replace("JournalTitle: ", "").strip()

        # Store only when both abbreviation and full name are available
        if current_abbr and current_full:
            journal_dict[current_abbr] = current_full
            current_abbr, current_full = None, None  # Reset for the next journal entry

    # Save dictionary locally
    with open(JOURNAL_LIST_FILE, "w", encoding="utf-8") as file:
        json.dump(journal_dict, file, ensure_ascii=False, indent=4)
    
    print(f"✅ Journal abbreviations saved to {JOURNAL_LIST_FILE}")
    
    return journal_dict

get_pubmed_journal_abbreviations()

JOURNAL_LIST_FILE = "/Users/ruizhang/Desktop/pubmed_journal_abbreviations.json"  # Local file path

def load_pubmed_journal_abbreviations():
    """
    Loads the journal abbreviations dictionary from a local JSON file.
    If the file does not exist, fetches and saves it first.
    """
    try:
        with open(JOURNAL_LIST_FILE, "r", encoding="utf-8") as file:
            journal_dict = json.load(file)
        print(f"📂 Loaded journal abbreviations from {JOURNAL_LIST_FILE}")
        return journal_dict
    except FileNotFoundError:
        print("⚠️ No local file found. Fetching data from PubMed...")
        return get_pubmed_journal_abbreviations()  # Fetch and save if file does not exist


def format_journal_abbreviation(journal, journal_dict):
    """
    Optimized function to validate and format a journal name or abbreviation for PubMed.
    - Uses sets and dictionaries for O(1) lookups to handle large datasets efficiently.
    - If a full name is given, retrieves and formats its abbreviation.
    - If an abbreviation is given, ensures proper formatting.
    - If neither is found, halts immediately with an exception.
    """
    
    # Preprocess the journal_dict into two sets for quick lookup
    full_names_set = set(journal_dict.values())
    abbreviations_set = set(journal_dict.keys())
    
    # Check if input is a known abbreviation
    if journal in abbreviations_set:
        formatted_abbr = journal + "." if " " in journal and not journal.endswith(".") else journal
        print(f"✅ Found abbreviation. Using formatted abbreviation: {formatted_abbr}")
        return formatted_abbr

    # Check if input is a full journal name
    elif journal in full_names_set:
        abbr = [abbr for abbr, full_name in journal_dict.items() if full_name == journal][0]
        formatted_abbr = abbr + "." if " " in abbr and not abbr.endswith(".") else abbr
        print(f"✅ Found full name. Using abbreviation: {formatted_abbr}")
        return formatted_abbr

    # Raise an exception if the journal is not found
    raise ValueError(f"❌ Error: '{journal}' not found in PubMed journal list.")





######################################################################
# create a separate function to fetch articles using the constructed query.
def fetch_article_ids_from_pubmed(query):
    """Fetch article IDs from PubMed based on the provided query."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="pub date")
    record = Entrez.read(handle)
    handle.close()

    return record["IdList"], len(record["IdList"])

######################################################################
# Create a function dedicated to building the query.
def build_pubmed_query(journal, start_date, end_date):
    """Construct the PubMed query string."""
    return (
        f'"{journal}"[Journal] AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication]) '
        f'AND ("journal article"[Publication Type] OR "review"[Publication Type]) '
        f'NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type])'
    )

######################################################################
# Construct PubMed Query Function
def build_pubmed_query(journal, start_date, end_date):
    """Construct the PubMed query string."""
    return (
        f'"{journal}"[Journal] AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication]) '
        f'AND ("journal article"[Publication Type] OR "review"[Publication Type]) '
        f'NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type])'
    )

######################################################################
# parse the fetched articles into a more useful format.
def parse_pubmed_article(paper_info):
    """Parse details from a retrieved PubMed article."""
    article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
    title = article_data["ArticleTitle"]
    abstract = article_data.get("Abstract", {}).get("AbstractText", ["No abstract available"])[0]
    journal_info = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]
    
    pub_date = journal_info.get("PubDate", {})
    pub_year = pub_date.get("Year", "Unknown Year")
    pub_month = pub_date.get("Month", "Unknown Month")
    pub_day = pub_date.get("Day", "")
    
    doi = next((id_item for id_item in paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
                if id_item.attributes["IdType"] == "doi"), None)

    return {
        "Publication Date": f"{pub_year}-{pub_month}-{pub_day}".strip("-"),
        "Title": title,
        "Abstract": abstract,
        "DOI": f"https://doi.org/{doi}" if doi else "No DOI available"
    }

######################################################################
def fetch_pubmed_articles_by_date(journal, start_date=None, end_date=None):

    # Load journal abbreviations
    journal_dict = load_pubmed_journal_abbreviations()
    formatted_journal = format_journal_abbreviation(journal, journal_dict)

    today = datetime.today()
    current_month = today.strftime("%Y-%m")

    if start_date is None:
        start_date = current_month
    if end_date is None:
        end_date = current_month

        print(f"No date provided, fetching articles from {journal} of the current month {today.strftime('%Y-%m')}.")

    # Validate the user provided dates
    validate_dates(start_date, end_date, today)

    # Construct dates for the API request only if start_date is in YYYY-MM format
    if len(start_date) == 7:  # YYYY-MM
        start_date = f"{start_date}/01"  # Construct start date as YYYY-MM/01
    elif len(start_date) == 10:  # YYYY-MM-DD
        start_date = start_date  # Use as is

    if len(end_date) == 7:
        end_date = f"{end_date}/{get_last_day_of_month(end_date[:4], end_date[5:])}"
    elif len(end_date) == 10:  # YYYY-MM-DD
        end_date = end_date  # Use as is

    # start_date = f"{start_date}/01"
    # end_date = (today.strftime("%Y/%m/%d")
    #             if end_date == current_date
    #             else f"{end_date}/{get_last_day_of_month(end_date[:4], end_date[5:])}")

    print(f"Date range used for query is from {start_date} to {end_date}")

    # Build the PubMed query
    query = build_pubmed_query(formatted_journal, start_date, end_date)
    print(f"Fetching articles from {formatted_journal} published between {start_date} and {end_date}...")

    # Fetch article IDs
    pmid_list, count = fetch_article_ids_from_pubmed(query)
    print("Fetched article IDs, proceeding to fetch article details...")

    papers = []
    if not pmid_list:
        print(f"No articles found for {formatted_journal} from {start_date} to {end_date}.")
        return papers

    print(f"✅ {count} papers found.")

    # Fetch details for each article
    total_papers = len(pmid_list)

    for index, pmid in enumerate(pmid_list):
        print(f"Fetching details for article {index + 1} of {total_papers} (PMID: {pmid})...")
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
        paper_info = Entrez.read(handle)
        handle.close()

        # Extract article details
        article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
        title = article_data["ArticleTitle"]
        abstract = article_data.get("Abstract", {}).get("AbstractText", ["No abstract available"])[0]
        journal_info = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]

        pub_date = journal_info.get("PubDate", {})
        pub_year = pub_date.get("Year", "Unknown Year")
        pub_month = pub_date.get("Month", "Unknown Month")
        pub_day = pub_date.get("Day", "")

        # Extract DOI if available
        doi = next((id_item for id_item in paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
                    if id_item.attributes["IdType"] == "doi"), None)

        papers.append({
            "Journal": journal,
            "Publication Date": f"{pub_year}-{pub_month}-{pub_day}".strip("-"),
            "Title": title,
            "Abstract": abstract,
            "DOI": f"https://doi.org/{doi}" if doi else "No DOI available"
        })

        # Print progress after fetching each article
        print(f"✅ Fetched details for article {index + 1}: {title}")

    print("😊 All articles have been fetched successfully.")
    return papers

######################################################################
# function to export the article list to csv file

def export_fetched_articles_as_csv(articles,journal,start_date,end_date):

    df = pd.DataFrame(articles)
    df.to_csv(f"JournalTracker_{journal}_{start_date}_to_{end_date}.csv", index=False)

    print(f"Fetched {len(df)} articles and saved to JournalTracker_{journal}_{start_date}_to_{end_date}.csv")



# testing
# articles = fetch_pubmed_articles_by_date("Environ Sci Technol","2025-02-10","2025-02-12")


# update:
# host the pubmed urc on github
# host the abbreviation to full name list on github



