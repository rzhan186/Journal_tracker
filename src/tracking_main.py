# import modules
from datetime import datetime
from Bio import Entrez
from calendar import monthrange

import pandas as pd
import re
import requests
import json
import difflib


# Set Entrez email
Entrez.email = "rzhan186@gmail.com"

######################################################################
# separate function for validating date inputs and checking the order.
def validate_date_input(date_str):
    """Validates user input for YYYY-MM or YYYY-MM-DD format."""
    pattern = r"^\d{4}-(0[1-9]|1[0-2])(-([0-2][0-9]|3[01]))?$"  # Matches YYYY-MM or YYYY-MM-DD
    return bool(re.match(pattern, date_str))


def validate_dates(start_date, end_date):

    today = datetime.today()

    """Validates the start and end date formats and order."""

    if not start_date:
        return
    
    if not validate_date_input(start_date):
        raise ValueError("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15).")

    if not validate_date_input(end_date):
        raise ValueError("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15).")
    
    if start_date > end_date:
        raise ValueError("❌ Start date cannot be later than end date.")

    if start_date > today.strftime('%Y-%m-%d'):
        raise ValueError("❌ Start date cannot be later than the current date")

    if end_date > today.strftime("%Y-%m-%d"):
        raise ValueError("❌ End date cannot be later than the current date")

# i want to modify the code, when a start date is not provided, dont do anything and skip the program
# if a start date is provide then validate start and end date and all else.


def get_last_day_of_month(year, month,):
    """Returns the last valid day of a given month."""
    return monthrange(int(year), int(month))[1]  # Correct last day

######################################################################
# Format journal names to allow full or abbreviations

JOURNAL_LIST_URL = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt" 
FALLBACK_JOURNAL_LIST_URL="https://raw.githubusercontent.com/rzhan186/Journal_tracker/refs/heads/main/data/pubmed_journal_abbreviations.json?token=GHSAT0AAAAAADF4U3Q7CM3EE27ZVXD7VDQM2CTP7OQ"

JOURNAL_LIST_FILE = "pubmed_journal_abbreviations.json"  

def get_pubmed_journal_abbreviations():
    """
    Fetches and parses the latest journal abbreviations from PubMed's official source.
    Saves the data locally in a JSON file for future use.
    Returns a dictionary {abbreviation: full_name}.
    """
    
    print("🔍 Fetching the latest journal abbreviations from PubMed...")

    try: 
        response = requests.get(JOURNAL_LIST_URL)
        response.raise_for_status()  # Ensure request was successful
        print("✅ Successfully fetched from primary URL.")

    except (requests.HTTPError, requests.ConnectionError) as e:
        print(f"⚠️ Primary URL failed: {e}. Trying fallback URL...")
        try:
            response = requests.get(FALLBACK_JOURNAL_LIST_URL)
            response.raise_for_status()  # Ensure request was successful
            print("✅ Successfully fetched from fallback URL.")
    
        except (requests.HTTPError, requests.ConnectionError) as e:
            print(f"⚠️ Fallback URL also failed: {e}. Unable to fetch journal abbreviations.")
            print(f"⚠️ You need to download the pubmed_journal_abbreviations.json file manually from https://github.com/rzhan186/Journal_tracker/blob/main/data/pubmed_journal_abbreviations.json and save as a local json file")
            return {}  # Return an empty dictionary if both requests fail
    
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
    
    print(f"✅ Journal abbreviations saved to the directory {JOURNAL_LIST_FILE}")
    
    return journal_dict

def load_pubmed_journal_abbreviations():
    """
    Loads the journal abbreviations dictionary from a local JSON file.
    If the file does not exist, fetches and saves it first from Github .
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
    - Provides suggestions for typos using fuzzy matching.
    """
    
    # Preprocess the journal_dict into two sets for quick lookup
    full_names_set = set(journal_dict.values())
    abbreviations_set = set(journal_dict.keys())
    
    # Check if input is a known abbreviation
    if journal in abbreviations_set:
        formatted_abbr = journal  # Use the abbreviation directly, no need to add a dot.
        print(f"✅ Found abbreviation. Using formatted abbreviation: {formatted_abbr}")
        return formatted_abbr

    # Check if input is a full journal name
    elif journal in full_names_set:
        abbr = [abbr for abbr, full_name in journal_dict.items() if full_name == journal][0]
        formatted_abbr = abbr  # Get the abbreviation directly without adding a dot
        print(f"✅ Found full name. Using abbreviation: {formatted_abbr}")
        return formatted_abbr

    # No exact match — try suggesting similar names
    possible_matches = difflib.get_close_matches(journal, list(full_names_set.union(abbreviations_set)), n=5, cutoff=0.6)
    
    suggestion_msg = ""
    if possible_matches:
        suggestion_msg = "\n🛈 Did you mean: " + ", ".join(f"'{s}'" for s in possible_matches)

    raise ValueError(f"❌ Error: '{journal}' not found in PubMed journal list.{suggestion_msg}")


######################################################################
# create a separate function to fetch articles using the constructed query.
def fetch_article_ids_from_pubmed(query):
    """Fetch article IDs from PubMed based on the provided query."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="pub date")
    record = Entrez.read(handle)
    handle.close()

    return record["IdList"], len(record["IdList"])


######################################################################
# create a function to allow web of science style boolean search
def format_boolean_keywords_for_pubmed(raw_query):
    """
    Formats a Boolean keyword query for PubMed Title/Abstract search.
    Preserves wildcards (*, ?), handles phrases, and attaches [Title/Abstract].
    """
    # Tokenize input: phrases in quotes, operators, parentheses, or single words/wildcards
    tokens = re.findall(r'"[^"]+"|\(|\)|\bAND\b|\bOR\b|\bNOT\b|\*|\w[\w*?\-]*', raw_query, flags=re.IGNORECASE)
    formatted_tokens = []

    for token in tokens:
        upper_token = token.upper()

        if upper_token in {"AND", "OR", "NOT", "(", ")"}:
            formatted_tokens.append(upper_token)
        elif token.startswith('"') and token.endswith('"'):
            # Phrase — preserve as is and tag
            phrase = token.strip('"')
            formatted_tokens.append(f'"{phrase}"[Title/Abstract]')
        else:
            # Single word or wildcarded term
            formatted_tokens.append(f'{token}[Title/Abstract]')

    return " ".join(formatted_tokens)

# example input
# (climat* OR "global warming") AND (mercury OR pollution)

# example output
# climat*[Title/Abstract] OR "global warming"[Title/Abstract] AND mercury[Title/Abstract] OR pollution[Title/Abstract]

######################################################################
# Create a function dedicated to building the query.
# handles formatted Boolean queries
def build_pubmed_query(journal, start_date, end_date, keywords=None):
    """Construct the PubMed query string, supporting Boolean keyword search."""
    base_query = (
        f'"{journal}"[Journal] AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication]) '
        f'AND ("journal article"[Publication Type] OR "review"[Publication Type]) '
        f'NOT ("news"[Publication Type] OR "comment"[Publication Type] OR "editorial"[Publication Type])'
    )

    if keywords:
        keyword_clause = f' AND ({keywords})'
        return base_query + keyword_clause
    else:
        return base_query


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
def fetch_pubmed_articles_by_date(journal, start_date=None, end_date=None,keywords=None):

    # Load journal abbreviations
    journal_dict = load_pubmed_journal_abbreviations()
    formatted_journal = format_journal_abbreviation(journal, journal_dict)

    today = datetime.today()
    current_month = today.strftime("%Y-%m")

    if not start_date:
        start_date = current_month
        end_date = current_month
        # Fixing bug if no start date or end date is provided, the program won't work
        print(f"No date provided, fetching articles from {journal} of the current month {today.strftime('%Y-%m')}.")
    else:
        validate_dates(start_date, end_date or current_month)

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
    query = build_pubmed_query(formatted_journal, start_date, end_date, keywords)

    # ✅ Print the final query for verification
    print(f"\n🔍 Here is the final PubMed Query:\n{query}\n in case you are curious")

    print(f"Fetching articles from {formatted_journal} published between {start_date} and {end_date} using keywords {keywords}")

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
# Create placeholder CSV in case no search was run
def generate_placeholder_csv():
    placeholder_data = {
        "Title": ["No search performed yet"],
        "Authors": ["N/A"],
        "Journal": ["N/A"],
        "Publication Date": [str(datetime.today().date())],
        "Abstract": ["This is a placeholder. Run a search to get real results."]
    }
    df_placeholder = pd.DataFrame(placeholder_data)
    return df_placeholder.to_csv(index=False).encode("utf-8")


######################################################################
# function to export the article list to csv file

def export_fetched_articles_as_csv(articles, journal, start_date, end_date, timestamp=None):
    df = pd.DataFrame(articles)

    # Generate filename with optional timestamp
    filename = f"JournalTracker_{journal}_{start_date}_to_{end_date}"
    if timestamp:
        filename += f"_{timestamp}"
    filename += ".csv"

    df.to_csv(filename, index=False)

    print(f"📁 Fetched {len(df)} articles and saved to {filename}")


# More update update:
    # prompty the user to enter an email address ➡️ implement this at the end
    # prompt the user to indicate an export directory, export to current directory by defacult
    # output publication date not aligned with the journal publication date. 

# implement a search by keyword function