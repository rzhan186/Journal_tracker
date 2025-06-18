
import re
import pandas as pd
from Bio import Entrez
from datetime import datetime
from calendar import monthrange
import requests
import json

# Set Entrez email
Entrez.email = "rzhan186@gmail.com"

JOURNAL_LIST_URL = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt"
JOURNAL_LIST_FILE = "/Users/ruizhang/Desktop/pubmed_journal_abbreviations.json"  # Local file path

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


def get_last_day_of_month(year, month):

    """Returns the last valid day of a given month."""
    return monthrange(int(year), int(month))[1]  # Correct last day


def fetch_pubmed_articles_by_date(journal, start_month=None, end_month=None):

    """
    Fetch research articles from a specific journal published within a given date range.

    - Defaults to the current month if no input is provided.
    - If fetching the current month, sets `end_date` to today's date.
    - Ensures valid date format and order.

    Parameters:
        journal (str): Name of the journal (e.g., "Nature").
        start_month (str, optional): Start month in YYYY-MM format (defaults to current month).
        end_month (str, optional): End month in YYYY-MM format (defaults to current month).

    Returns:
        List of dictionaries containing article details (title, abstract, DOI).
    """
    # check if the entered journal name is valid
    # Load the stored journal dictionary
    journal_dict = load_pubmed_journal_abbreviations()

    # Validate and format journal abbreviation
    formatted_journal = format_journal_abbreviation(journal, journal_dict)
    
    # Get current month dynamically (YYYY-MM), if no values provided, use the default date
    
    today = datetime.today()
    current_month = today.strftime("%Y-%m")

    if start_month is None:
        start_month = today.strftime("%Y-%m")

    if end_month is None:
        end_month = today.strftime("%Y-%m")

        print(f"No date provided, fetching articles from {journal} of the current month {today.strftime("%Y-%m")}." )

    def validate_date_input(date_str):
        """Validates user input for YYYY-MM format."""
        pattern = r"^\d{4}-(0[1-9]|1[0-2])$"  # Ensures format YYYY-MM with months 01-12
        return bool(re.match(pattern, date_str))


    # Validate provided dates
    if not validate_date_input(start_month) or not validate_date_input(end_month):
        raise ValueError("❌ Invalid format! Dates must be in YYYY-MM format (e.g., 2024-01).")

    if start_month > end_month:
        raise ValueError("❌ Start month cannot be later than end month.")

    if start_month > today.strftime("%Y-%m"):
        raise ValueError("❌ Start month cannot be later than the current month")
    
    if end_month > today.strftime("%Y-%m"):
        raise ValueError("❌ End month cannot be later than the current month")

    # Convert YYYY-MM to YYYY/MM/DD for PubMed query
    # start_date = f"{start_month}/01"
    # end_date = f"{end_month}/31"  # Using 31 is safe; PubMed adjusts for month length

    # If no date is provided, ensure the end date is today's date

    if end_month == current_month:
        # If using the current month, set the end date to today's date
        end_date = today.strftime("%Y/%m/%d") 
    else:
        end_date = f"{end_month}/{get_last_day_of_month(start_month[:4], start_month[5:])}"

    # Always start from the first of the start_month
    start_date = f"{start_month}/01"

    # Construct PubMed query with date filter. Includes both research and review articles.
    query = (f'"{journal}"[Journal] AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication]) '
            f'AND ("journal article"[Publication Type] OR "review"[Publication Type]) '
            f'NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type])')

    print(f"Fetching articles from {journal} published between {start_date} and {end_date}...")

    # Fetch article IDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="pub date")
    record = Entrez.read(handle)
    handle.close()

    pmid_list = record["IdList"]
    papers = []

    if not pmid_list:
        print(f"No articles found for {journal} from {start_date} to {end_date}.")
        return papers

    print(f"✅ {len(pmid_list)} papers found.")
    
    # Fetch details for each article
    for pmid in pmid_list:
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

    return papers


# include an additional functionality to allow for keyword search

# Example Usage
journal = "Nat Rev Microbiol."
articles = fetch_pubmed_articles_by_date({journal}, "2024-10", "2025-01")

# Fetch articles from Nature published in this month
articles = fetch_pubmed_articles_by_date("Nature")

# Fetch articles from Nature published between January and March 2024
articles = fetch_pubmed_articles_by_date("Nature", "2024-02", "2024-03")

# Save results to CSV
df = pd.DataFrame(articles)
df.to_csv(f"{journal}_by_date.csv", index=False)

print(f"Fetched {len(df)} articles and saved to Journal_tracker_{journal}_by_date.csv")


