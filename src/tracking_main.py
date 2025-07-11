# tracking_main.py

# import modules
import os
from datetime import datetime, timedelta
from Bio import Entrez
from calendar import monthrange
import pandas as pd
import re
import requests
import json
import difflib
import time
from dotenv import load_dotenv
from concurrent.futures import ThreadPoolExecutor
import threading
import random

from RateLimit import PubMedRateLimit

# Load secrets
load_dotenv()

# Configure Entrez with API key and email
Entrez.email = os.getenv("PUBMED_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

# Verify configuration
if not Entrez.email:
    print("⚠️ Warning: PUBMED_EMAIL not found in environment variables")
if not Entrez.api_key:
    print("⚠️ Warning: NCBI_API_KEY not found in environment variables")
else:
    print("✅ NCBI API key configured successfully")


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

def format_journal_abbreviation(journal, journal_dict, rate_limiter=None):
    """
    Enhanced function to validate and format a journal name or abbreviation for PubMed.
    - First checks manual mappings for common journals
    - Then checks official PubMed list
    - Finally tests the journal name directly with PubMed
    - Provides suggestions for typos using fuzzy matching
    """
    
    # Manual mappings for commonly searched journals
    manual_mappings = {
        # Lancet family
        "Lancet": "Lancet",
        "The Lancet": "Lancet",
        "BMJ": "BMJ"
        
        # Add more as needed...
    }
    
    # Step 1: Check manual mappings first
    if journal in manual_mappings:
        formatted_abbr = manual_mappings[journal]
        print(f"✅ Found in manual mappings. Using: {formatted_abbr}")
        return formatted_abbr
    
    # Step 2: Check official PubMed list
    abbreviations_set = set(journal_dict.keys())
    full_names_set = set(journal_dict.values())
    
    # Check if input is a known abbreviation
    if journal in abbreviations_set:
        formatted_abbr = journal
        print(f"✅ Found abbreviation in official list. Using: {formatted_abbr}")
        return formatted_abbr

    # Check if input is a full journal name
    if journal in full_names_set:
        abbr = [abbr for abbr, full_name in journal_dict.items() if full_name == journal][0]
        formatted_abbr = abbr
        print(f"✅ Found full name in official list. Using abbreviation: {formatted_abbr}")
        return formatted_abbr
    
    # Step 3: Test directly with PubMed (for journals like "Lancet")
    print(f"🔍 Testing '{journal}' directly with PubMed...")
        
    if test_journal_name_on_pubmed(journal, rate_limiter):
            print(f"✅ '{journal}' works directly with PubMed!")
            return journal

    if test_journal_name_on_pubmed(journal):
        print(f"✅ '{journal}' works directly with PubMed!")
        return journal
    
    # Step 4: If nothing works, provide suggestions
    all_possible_names = list(
        full_names_set.union(abbreviations_set).union(manual_mappings.keys())
    )
    possible_matches = difflib.get_close_matches(journal, all_possible_names, n=5, cutoff=0.6)
    
    suggestion_msg = ""
    if possible_matches:
        suggestion_msg = "\n🛈 Did you mean: " + ", ".join(f"'{s}'" for s in possible_matches)

    raise ValueError(f"❌ Error: '{journal}' not found in any source.{suggestion_msg}")

def test_journal_name_on_pubmed(journal_name, rate_limiter=None):
    """
    Test if a journal name works directly with PubMed by doing a small search
    """
    try:
        # Test with a recent date range and limit to 1 result
        test_query = f'"{journal_name}"[Journal] AND ("2023/01/01"[Date - Publication] : "2024/12/31"[Date - Publication])'
        
        if rate_limiter:
            # Use the safe search method
            record = rate_limiter.safe_pubmed_search(test_query, max_results=1)
        else:
            # Fallback to direct Entrez call
            handle = Entrez.esearch(db="pubmed", term=test_query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
        
        # If we get any results, the journal name works
        return record and len(record.get("IdList", [])) > 0
        
    except Exception as e:
        # If there's an API error, assume the journal name doesn't work
        return False

def fetch_article_ids_from_pubmed(query, rate_limiter=None):
    """Fetch article IDs from PubMed based on the provided query."""
    
    if rate_limiter:
        # Use the safe search method
        record = rate_limiter.safe_pubmed_search(query, max_results=1000)
        if record:
            return record["IdList"], len(record["IdList"])
        else:
            return [], 0
    else:
        # Fallback to direct Entrez call
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

######################################################################
# Create a function dedicated to building the query.
# handles formatted Boolean queries
def build_pubmed_query(journal, start_date, end_date, keywords=None):
    """Construct the PubMed query string, supporting Boolean keyword search."""

    # Primary inclusion criteria - what we WANT
    include_types = [
        '"journal article"[Publication Type]',
        '"review"[Publication Type]',
        '"systematic review"[Publication Type]',
        '"meta-analysis"[Publication Type]',
        '"clinical trial"[Publication Type]',
        '"randomized controlled trial"[Publication Type]',
        '"comparative study"[Publication Type]',
        '"multicenter study"[Publication Type]'
    ]
    
    # Comprehensive exclusion criteria - what we DON'T want
    exclude_types = [
        '"news"[Publication Type]',
        '"comment"[Publication Type]',
        '"editorial"[Publication Type]',
        '"letter"[Publication Type]'
    ]
    
    # Build the base query with abstract requirement
    base_query = f'"{journal}"[Journal] AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication]) AND hasabstract[text]'
    
    # Add inclusion criteria (OR together)
    inclusion_clause = f' AND ({" OR ".join(include_types)})'
    
    # Add exclusion criteria (NOT each one)
    exclusion_clause = f' NOT ({" OR ".join(exclude_types)})'

    # Combine everything
    full_query = base_query + inclusion_clause + exclusion_clause

    # Add keywords if provided
    if keywords:
        keyword_clause = f' AND ({keywords})'
        full_query += keyword_clause
    
    return full_query


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


def fetch_pubmed_articles_by_date(journal, start_date=None, end_date=None, keywords=None, rate_limiter=None):
    from datetime import datetime
    import time

    # Load journal abbreviations
    journal_dict = load_pubmed_journal_abbreviations()
    formatted_journal = format_journal_abbreviation(journal, journal_dict, rate_limiter)

    today = datetime.today()
    current_month = today.strftime("%Y-%m")

    if not start_date:
        start_date = current_month
        end_date = current_month
        print(f"No date provided, fetching articles from {journal} of the current month {today.strftime('%Y-%m')}.")
    else:
        validate_dates(start_date, end_date or current_month)

    # Format dates for query
    if len(start_date) == 7:  # YYYY-MM
        start_date = f"{start_date}/01"
    if len(end_date) == 7:
        end_date = f"{end_date}/{get_last_day_of_month(end_date[:4], end_date[5:])}"

    print(f"Date range used for query is from {start_date} to {end_date}")

    # Build query
    query = build_pubmed_query(formatted_journal, start_date, end_date, keywords)

    print(f"\n🔍 Final PubMed Query:\n{query}\n")
    print(f"Fetching articles from {formatted_journal} between {start_date} and {end_date} using keywords {keywords}")

    # Fetch article IDs
    pmid_list, count = fetch_article_ids_from_pubmed(query, rate_limiter)
    print("Fetched article IDs, proceeding to fetch article details...")

    if not pmid_list:
        print(f"No articles found for {formatted_journal} from {start_date} to {end_date}.")
        return []

    print(f"✅ {count} papers found.")
    papers = []

    for index, pmid in enumerate(pmid_list):
        print(f"Fetching details for article {index + 1} of {len(pmid_list)} (PMID: {pmid})...")

        try:
            if rate_limiter:
                paper_xml = rate_limiter.safe_pubmed_fetch([pmid])
                if not paper_xml:
                    continue
                time.sleep(1)  # Optional manual delay
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            paper_info = Entrez.read(handle)
            handle.close()

            article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
            title = article_data["ArticleTitle"]
            
            # Enhanced abstract handling
            abstract_sections = article_data.get("Abstract", {}).get("AbstractText", [])
            if abstract_sections:
                if isinstance(abstract_sections[0], dict):
                    abstract = " ".join([section.get("text", "") for section in abstract_sections])
                else:
                    abstract = str(abstract_sections[0])
            else:
                abstract = "No abstract available"
            
            # Enhanced author extraction
            authors = []
            author_list = article_data.get("AuthorList", [])
            for author in author_list:
                last_name = author.get("LastName", "")
                first_name = author.get("ForeName", "")
                if last_name:
                    authors.append(f"{last_name}, {first_name}")
            
            # Journal information
            journal_info = article_data["Journal"]
            journal_title = journal_info.get("Title", journal)
            journal_abbrev = journal_info.get("ISOAbbreviation", "")
            
            # Issue information
            journal_issue = journal_info.get("JournalIssue", {})
            volume = journal_issue.get("Volume", "")
            issue = journal_issue.get("Issue", "")
            
            # Enhanced publication date handling
            pub_date = journal_issue.get("PubDate", {})
            pub_year = pub_date.get("Year", "Unknown Year")
            pub_month = pub_date.get("Month", "Unknown Month")
            pub_day = pub_date.get("Day", "")
            
            # Pages
            pagination = article_data.get("Pagination", {})
            pages = pagination.get("MedlinePgn", "")
            
            # Enhanced DOI handling
            doi = None
            article_ids = paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
            for id_item in article_ids:
                if id_item.attributes["IdType"] == "doi":
                    doi = str(id_item)
                    break
            
            # PMID
            pmid_value = paper_info["PubmedArticle"][0]["MedlineCitation"]["PMID"]

            papers.append({
                "Journal": journal_title or journal_abbrev or journal,
                "Publication Date": f"{pub_year}-{pub_month}-{pub_day}".strip("-"),
                "Title": str(title),
                "Authors": " and ".join(authors) if authors else "Unknown Author",
                "Abstract": abstract,
                "Volume": volume,
                "Issue": issue,
                "Pages": pages,
                "Year": str(pub_year),
                "Month": str(pub_month),
                "Day": str(pub_day),
                "DOI": f"https://doi.org/{doi}" if doi else "No DOI available",
                "PMID": str(pmid_value)
            })

            print(f"✅ Fetched details for article {index + 1}: {title}")

        except Exception as e:
            print(f"⚠️ Failed to fetch article {pmid}: {e}")
            continue

    print("😊 All articles have been fetched successfully.")
    return papers


def parse_pubmed_article_enhanced(paper_info):
    """Enhanced parser that extracts more detailed publication information"""
    try:
        article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
        
        # Title
        title = article_data["ArticleTitle"]
        
        # Authors
        authors = []
        author_list = article_data.get("AuthorList", [])
        for author in author_list:
            last_name = author.get("LastName", "")
            first_name = author.get("ForeName", "")
            if last_name:
                authors.append(f"{last_name}, {first_name}")
        
        # Journal information
        journal_info = article_data["Journal"]
        journal_title = journal_info.get("Title", "")
        journal_abbrev = journal_info.get("ISOAbbreviation", "")
        
        # Issue information
        journal_issue = journal_info.get("JournalIssue", {})
        volume = journal_issue.get("Volume", "")
        issue = journal_issue.get("Issue", "")
        
        # Publication date
        pub_date = journal_issue.get("PubDate", {})
        pub_year = pub_date.get("Year", "")
        pub_month = pub_date.get("Month", "")
        pub_day = pub_date.get("Day", "")
        
        # Pages
        pagination = article_data.get("Pagination", {})
        pages = pagination.get("MedlinePgn", "")
        
        # Abstract
        abstract = article_data.get("Abstract", {}).get("AbstractText", [""])[0]
        
        # DOI
        doi = None
        article_ids = paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]
        for id_item in article_ids:
            if id_item.attributes["IdType"] == "doi":
                doi = str(id_item)
                break
        
        # PMID
        pmid = paper_info["PubmedArticle"][0]["MedlineCitation"]["PMID"]
        
        return {
            "Title": title,
            "Authors": " and ".join(authors) if authors else "Unknown Author",
            "Journal": journal_title or journal_abbrev,
            "Volume": volume,
            "Issue": issue,
            "Pages": pages,
            "Year": pub_year,
            "Month": pub_month,
            "Day": pub_day,
            "Abstract": abstract,
            "DOI": doi,
            "PMID": str(pmid),
            "Publication Date": f"{pub_year}-{pub_month}-{pub_day}".strip("-")
        }
        
    except Exception as e:
        print(f"Error parsing article: {e}")
        return None



######################################################################
# function to allow fetch preprints to use the same keyword logic as pubmed

def compile_keyword_filter(raw_query, debug=False):
    """
    Compiles a Boolean keyword filter with debugging support.
    """
    
    def tokenize(query):
        tokens = re.findall(r'\(|\)|\bAND\b|\bOR\b|\bNOT\b|"[^"]+"|\w[\w*?]*', query, flags=re.IGNORECASE)
        if debug:
            print(f"DEBUG - Tokens: {tokens}")
        return tokens

    def to_python_expr(tokens):
        expr_parts = []
        for i, token in enumerate(tokens):
            upper = token.upper()
            if upper in {"AND", "OR", "NOT"}:
                # Handle implicit AND before NOT
                if upper == "NOT" and i > 0 and tokens[i-1].upper() not in {"AND", "OR", "("}:
                    expr_parts.append(" and ")
                expr_parts.append(f" {upper.lower()} ")
            elif token in "()":
                expr_parts.append(token)
            else:
                # Add implicit AND if needed
                if (i > 0 and tokens[i-1].upper() not in {"AND", "OR", "NOT", "("} and 
                    tokens[i-1] != "("):
                    expr_parts.append(" and ")
                
                clean_token = token.strip('"')
                escaped_token = re.escape(clean_token).replace(r'\*', '.*').replace(r'\?', '.')
                expr_parts.append(f'bool(re.search(r"{escaped_token}", text, re.IGNORECASE))')
        
        expr = ''.join(expr_parts)
        if debug:
            print(f"DEBUG - Expression: {expr}")
        return expr

    def matcher(text):
        try:
            tokens = tokenize(raw_query)
            expr = to_python_expr(tokens)
            
            if not expr.strip():
                return True
            
            result = eval(expr, {"re": re, "text": text, "bool": bool})
            
            if debug:
                print(f"DEBUG - Query: '{raw_query}'")
                print(f"DEBUG - Text preview: '{text[:100]}...'")
                print(f"DEBUG - Result: {result}")
                
                # Additional debugging for your specific case
                if 'cancer' in raw_query.lower() and 'polyamines' in raw_query.lower():
                    has_cancer = bool(re.search(r"cancer", text, re.IGNORECASE))
                    has_polyamines = bool(re.search(r"polyamines", text, re.IGNORECASE))
                    print(f"DEBUG - Text has 'cancer': {has_cancer}")
                    print(f"DEBUG - Text has 'polyamines': {has_polyamines}")
                    print(f"DEBUG - Expected result: {has_cancer and not has_polyamines}")
            
            return result
        except Exception as e:
            if debug:
                print(f"DEBUG - Error: {e}")
            return False

    return matcher

######################################################################
# function to fetch articles from preprints
import xml.etree.ElementTree as ET

# def fetch_preprints(server="biorxiv", start_date=None, end_date=None, keywords=None, max_results=50):
#     """
#     Fetches preprints from bioRxiv, medRxiv, or arXiv in a standardized format.
    
#     Args:
#         server (str): 'biorxiv', 'medrxiv', or 'arxiv'.
#         start_date (str): Start date in YYYY-MM-DD (bioRxiv/medRxiv only).
#         end_date (str): End date in YYYY-MM-DD (bioRxiv/medRxiv only).
#         keywords (str): Boolean keyword query (AND, OR, NOT, *, ?, parentheses).
#         max_results (int): Limit for arXiv results.

#     Returns:
#         List[Dict]: Each with keys Journal, Title, Abstract, Publication Date, DOI.
#     """
#     if server not in ["biorxiv", "medrxiv", "arxiv"]:
#         print("❌ Invalid server. Choose from 'biorxiv', 'medrxiv', or 'arxiv'.")
#         return []

#     matcher = compile_keyword_filter(keywords) if keywords else lambda x: True
#     results = []

#     if server == "arxiv":
#         if not keywords:
#             print("⚠️ arXiv requires keywords.")
#             return []
#         base_url = "http://export.arxiv.org/api/query"
#         query = f"all:{keywords}"
#         params = {
#             "search_query": query,
#             "start": 0,
#             "max_results": max_results,
#             "sortBy": "submittedDate",
#             "sortOrder": "descending"
#         }

#         try:
#             response = requests.get(base_url, params=params, timeout=10)
#             response.raise_for_status()

#             root = ET.fromstring(response.text)
#             ns = {'atom': 'http://www.w3.org/2005/Atom'}

#             for entry in root.findall('atom:entry', ns):
#                 title = entry.find('atom:title', ns).text.strip()
#                 summary = entry.find('atom:summary', ns).text.strip()
#                 published = entry.find('atom:published', ns).text.strip()
#                 doi = None
#                 for link in entry.findall('atom:link', ns):
#                     if 'doi.org' in link.get('href', ''):
#                         doi = link.get('href')
#                         break

#                 combined_text = f"{title} {summary}"
#                 if matcher(combined_text):
#                     results.append({
#                         "Journal": "arXiv",
#                         "Publication Date": published[:10],
#                         "Title": title,
#                         "Abstract": summary,
#                         "DOI": doi or "N/A"
#                     })
#         except Exception as e:
#             print(f"⚠️ Error fetching from arXiv: {e}")
#             return []

#     else:
#         if not (start_date and end_date):
#             print("⚠️ Start and end dates are required for bioRxiv/medRxiv.")
#             return []

#         url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}"
#         try:
#             response = requests.get(url, timeout=10)
#             response.raise_for_status()
#             data = response.json()

#             for item in data.get("collection", []):
#                 title = item.get("title", "")
#                 abstract = item.get("abstract", "")
#                 combined_text = f"{title} {abstract}"
#                 if matcher(combined_text):
#                     results.append({
#                         "Journal": server,
#                         "Publication Date": item["date"],
#                         "Title": title,
#                         "Abstract": abstract,
#                         "DOI": item.get("doi", "N/A")
#                     })
#         except Exception as e:
#             print(f"⚠️ Error fetching from {server}: {e}")
#             return []

#     print(f"✅ Fetched {len(results)} preprints from {server}")
#     return results

def fetch_preprints(server="biorxiv", start_date=None, end_date=None, keywords=None, max_results=None):
    """
    Fetches preprints from bioRxiv or medRxiv with enhanced functionality.
    
    Args:
        server (str): 'biorxiv' or 'medrxiv'
        start_date (str): Start date in YYYY-MM-DD format.
        end_date (str): End date in YYYY-MM-DD format.
        keywords (str): Boolean keyword query (AND, OR, NOT, *, ?, parentheses).
        max_results (int): Maximum results to fetch (None = unlimited).

    Returns:
        List[Dict]: Each with keys Journal, Title, Abstract, Publication Date, DOI.
    """
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    # Create keyword matcher if keywords provided
    matcher = compile_keyword_filter(keywords) if keywords else lambda x: True
    results = []

    if not (start_date and end_date):
        print("⚠️ Start and end dates are required for bioRxiv/medRxiv.")
        return []

    # Ensure dates are in YYYY-MM-DD format
    try:
        # Convert dates to proper format if needed
        if isinstance(start_date, str):
            if len(start_date) == 7:  # YYYY-MM format
                start_date = f"{start_date}-01"
            # Validate date format
            datetime.strptime(start_date, '%Y-%m-%d')
        
        if isinstance(end_date, str):
            if len(end_date) == 7:  # YYYY-MM format
                # Get last day of month
                year, month = end_date.split('-')
                last_day = monthrange(int(year), int(month))[1]
                end_date = f"{end_date}-{last_day:02d}"
            # Validate date format
            datetime.strptime(end_date, '%Y-%m-%d')
            
    except ValueError as e:
        print(f"❌ Invalid date format: {e}")
        return []

    # Pagination to get all results
    print(f"🔍 Fetching from {server}: {start_date} to {end_date}")
    
    try:
        initial_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/0"
        response = requests.get(initial_url, timeout=30)
        response.raise_for_status()
        
        if not response.headers.get('content-type', '').startswith('application/json'):
            print(f"❌ Non-JSON response from {server}")
            return []

        initial_data = response.json()
        
        if "collection" not in initial_data:
            print(f"⚠️ No 'collection' field in {server} response")
            return []

        messages = initial_data.get("messages", [])
        total_papers = None
        if messages and "total" in messages[0]:
            total_papers = int(messages[0]["total"])
            print(f"📊 {server}: Found {total_papers} preprints matching your search.")
        else:
            total_papers = len(initial_data.get("collection", []))
            print(f"📊 {server}: Found {total_papers} preprints (estimated; could be many more).")
                
        cursor = 0
        page_size = 100
        all_preprints = []
        
        while True:
            paginated_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/{cursor}"
            try:
                response = requests.get(paginated_url, timeout=30)
                response.raise_for_status()
                
                if not response.headers.get('content-type', '').startswith('application/json'):
                    print(f"❌ Non-JSON response from {server}")
                    break

                data = response.json()
                collection = data.get("collection", [])
                
                if not collection:
                    print(f"📭 No more results from {server} at cursor {cursor}")
                    break
                
                page_matches = 0
                for item in collection:
                    try:
                        title = item.get("title", "").strip()
                        abstract = item.get("abstract", "").strip()
                        
                        if not title and not abstract:
                            continue
                        
                        combined_text = f"{title} {abstract}"
                        if matcher(combined_text):
                            doi = item.get("doi", "")
                            if doi and not doi.startswith("http"):
                                doi = f"https://doi.org/{doi}"
                            
                            all_preprints.append({
                                "Journal": server,
                                "Publication Date": item.get("date", ""),
                                "Title": title,
                                "Abstract": abstract,
                                "DOI": doi if doi else "N/A"
                            })
                            page_matches += 1
                            
                    except Exception as e:
                        print(f"⚠️ Error processing preprint item: {e}")
                        continue
                
                print(f"📄 {server}: Page {cursor//page_size + 1} - {len(collection)} papers, {page_matches} matched keywords")
                
                cursor += len(collection)
                
                if len(collection) < page_size:
                    print(f"✅ {server}: Reached end of results")
                    break
                
                time.sleep(0.5)
                
                if max_results and len(all_preprints) >= max_results:
                    print(f"🛑 {server}: Reached max_results limit ({max_results})")
                    all_preprints = all_preprints[:max_results]
                    break
                    
            except requests.exceptions.Timeout:
                print(f"⏱️ Timeout fetching page from {server}")
                break
            except requests.exceptions.RequestException as e:
                print(f"⚠️ Request error fetching page from {server}: {e}")
                break
            except Exception as e:
                print(f"⚠️ Unexpected error fetching page from {server}: {e}")
                break
        
        results.extend(all_preprints)

    except requests.exceptions.Timeout:
        print(f"⏱️ Initial timeout fetching from {server}")
        return []
    except requests.exceptions.RequestException as e:
        print(f"⚠️ Initial request error fetching from {server}: {e}")
        return []
    except Exception as e:
        print(f"⚠️ Unexpected initial error fetching from {server}: {e}")
        return []

    print(f"✅ Successfully fetched {len(results)} preprints from {server}")
    return results

def fetch_preprints_parallel_pages(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """
    Fetch preprint pages in parallel to speed up data retrieval
    """
    import concurrent.futures
    from concurrent.futures import ThreadPoolExecutor
    
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    # Setup keyword matcher
    if keywords and keywords.strip():
        try:
            matcher = compile_keyword_filter(keywords.strip())
            has_keyword_filter = True
        except Exception as e:
            print(f"⚠️ Error compiling keyword filter: {e}")
            matcher = lambda x: True
            has_keyword_filter = False
    else:
        matcher = lambda x: True
        has_keyword_filter = False

    # Get total count first
    try:
        initial_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/0"
        response = requests.get(initial_url, timeout=30)
        response.raise_for_status()
        initial_data = response.json()
        
        messages = initial_data.get("messages", [])
        total_papers = int(messages[0]["total"]) if messages and "total" in messages[0] else len(initial_data.get("collection", []))
        
        print(f"📊 {server}: Found {total_papers} total preprints")
        
        # Calculate page ranges for parallel processing
        page_size = 100
        total_pages = (total_papers + page_size - 1) // page_size
        
        # Limit concurrent requests to avoid overwhelming the server
        max_workers = min(5, total_pages)
        
        def fetch_single_page(cursor):
            """Fetch a single page of results"""
            try:
                url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/{cursor}"
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                data = response.json()
                collection = data.get("collection", [])
                
                page_results = []
                for item in collection:
                    try:
                        title = item.get("title", "").strip()
                        abstract = item.get("abstract", "").strip()
                        
                        if not title and not abstract:
                            continue
                        
                        # Quick keyword filtering
                        combined_text = f"{title} {abstract}"
                        if matcher(combined_text):
                            doi = item.get("doi", "")
                            if doi and not doi.startswith("http"):
                                doi = f"https://doi.org/{doi}"
                            
                            page_results.append({
                                "Journal": server,
                                "Publication Date": item.get("date", ""),
                                "Title": title,
                                "Abstract": abstract,
                                "DOI": doi if doi else "N/A"
                            })
                    except Exception as e:
                        continue
                
                return page_results, len(collection)
                
            except Exception as e:
                print(f"⚠️ Error fetching page at cursor {cursor}: {e}")
                return [], 0

        # Process pages in parallel
        all_results = []
        total_processed = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all page fetch tasks
            cursor_values = list(range(0, total_papers, page_size))
            future_to_cursor = {
                executor.submit(fetch_single_page, cursor): cursor 
                for cursor in cursor_values
            }
            
            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_cursor):
                cursor = future_to_cursor[future]
                try:
                    page_results, page_count = future.result()
                    all_results.extend(page_results)
                    total_processed += page_count
                    
                    # Update progress
                    if progress_callback:
                        progress_callback(len(all_results))
                    
                    print(f"📄 {server}: Page {cursor//page_size + 1} - {page_count} processed, {len(page_results)} matched")
                    
                    # Early termination if we have enough results
                    if max_results and len(all_results) >= max_results:
                        print(f"🛑 {server}: Reached max_results limit ({max_results})")
                        # Cancel remaining futures
                        for remaining_future in future_to_cursor:
                            if remaining_future != future:
                                remaining_future.cancel()
                        break
                        
                except Exception as e:
                    print(f"⚠️ Error processing page: {e}")
                    continue
        
        # Limit results if max_results specified
        if max_results and len(all_results) > max_results:
            all_results = all_results[:max_results]
        
        print(f"✅ {server}: Found {len(all_results)} matching preprints out of {total_processed} processed")
        return all_results
        
    except Exception as e:
        print(f"⚠️ Error in parallel fetch from {server}: {e}")
        return []

def fetch_preprints_smart_filter(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """
    Optimized preprint search with smart filtering and early termination
    """
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    # Pre-compile keyword patterns for faster matching
    keyword_patterns = []
    if keywords and keywords.strip():
        try:
            # Extract individual terms for quick pre-filtering
            terms = re.findall(r'\b\w+\b', keywords.lower())
            keyword_patterns = [re.compile(term, re.IGNORECASE) for term in terms if len(term) > 2]
            matcher = compile_keyword_filter(keywords.strip())
            has_filters = True
        except Exception as e:
            print(f"⚠️ Error compiling filters: {e}")
            matcher = lambda x: True
            has_filters = False
    else:
        matcher = lambda x: True
        has_filters = False

    def quick_filter(text):
        """Quick pre-filter to eliminate obviously non-matching articles"""
        if not has_filters:
            return True
        
        # Quick check: if none of the key terms appear, skip expensive matching
        text_lower = text.lower()
        if keyword_patterns and not any(pattern.search(text_lower) for pattern in keyword_patterns):
            return False
        
        return True

    results = []
    
    try:
        # Get initial data
        initial_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/0"
        response = requests.get(initial_url, timeout=30)
        response.raise_for_status()
        initial_data = response.json()
        
        messages = initial_data.get("messages", [])
        total_papers = int(messages[0]["total"]) if messages and "total" in messages[0] else len(initial_data.get("collection", []))
        
        print(f"📊 {server}: Found {total_papers} total preprints")
        
        cursor = 0
        page_size = 100
        total_processed = 0
        quick_filtered = 0
        
        while cursor < total_papers:
            try:
                url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/{cursor}"
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                data = response.json()
                collection = data.get("collection", [])
                
                if not collection:
                    break
                
                page_matches = 0
                for item in collection:
                    try:
                        title = item.get("title", "").strip()
                        abstract = item.get("abstract", "").strip()
                        
                        if not title and not abstract:
                            continue
                        
                        combined_text = f"{title} {abstract}"
                        
                        # Quick pre-filter
                        if not quick_filter(combined_text):
                            continue
                        
                        quick_filtered += 1
                        
                        # Full keyword matching
                        if matcher(combined_text):
                            doi = item.get("doi", "")
                            if doi and not doi.startswith("http"):
                                doi = f"https://doi.org/{doi}"
                            
                            results.append({
                                "Journal": server,
                                "Publication Date": item.get("date", ""),
                                "Title": title,
                                "Abstract": abstract,
                                "DOI": doi if doi else "N/A"
                            })
                            page_matches += 1
                        
                        total_processed += 1
                        
                        # Update progress every 50 articles
                        if total_processed % 50 == 0 and progress_callback:
                            progress_callback(len(results))
                        
                    except Exception as e:
                        continue
                
                print(f"📄 {server}: Page {cursor//page_size + 1} - {len(collection)} processed, {page_matches} matched")
                
                cursor += len(collection)
                
                # Early termination if we have enough results
                if max_results and len(results) >= max_results:
                    print(f"🛑 {server}: Reached max_results limit ({max_results})")
                    results = results[:max_results]
                    break
                
                # Small delay to be nice to the API
                time.sleep(0.1)
                
            except Exception as e:
                print(f"⚠️ Error fetching page: {e}")
                break
        
        print(f"✅ {server}: {len(results)} matches found")
        print(f"📊 Filtering efficiency: {quick_filtered}/{total_processed} passed pre-filter")
        
        return results
        
    except Exception as e:
        print(f"⚠️ Error in smart filter from {server}: {e}")
        return [] 

def fetch_preprints_with_progress(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """
    Single optimized approach - parallel processing with smart filtering
    """
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    # Always use parallel processing (fastest approach)
    return fetch_preprints_parallel_pages(server, start_date, end_date, keywords, max_results, progress_callback)


######################################################################
# implement caching for repeated searches
import hashlib
import pickle
import os
from datetime import datetime, timedelta

def get_cache_key(server, start_date, end_date, keywords):
    """Generate cache key for search parameters"""
    cache_data = f"{server}_{start_date}_{end_date}_{keywords or ''}"
    return hashlib.md5(cache_data.encode()).hexdigest()
def fetch_preprints_with_cache(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """
    Cached version of preprint search - saves results for 1 hour
    """
    # Check cache first
    cache_key = get_cache_key(server, start_date, end_date, keywords)
    cache_file = f"cache/preprints_{cache_key}.pkl"
    
    if os.path.exists(cache_file):
        # Check if cache is less than 1 hour old
        if time.time() - os.path.getmtime(cache_file) < 3600:
            try:
                with open(cache_file, 'rb') as f:
                    cached_results = pickle.load(f)
                print(f"📋 Using cached results for {server} ({len(cached_results)} articles)")
                if progress_callback:
                    progress_callback(len(cached_results))
                return cached_results
            except Exception as e:
                print(f"⚠️ Cache read error: {e}")
    
    # Fetch fresh results using the optimized function
    results = fetch_preprints_smart_filter(server, start_date, end_date, keywords, max_results, progress_callback)
    
    # Cache the results
    try:
        os.makedirs("cache", exist_ok=True)
        with open(cache_file, 'wb') as f:
            pickle.dump(results, f)
        print(f"💾 Cached {len(results)} results for future use")
    except Exception as e:
        print(f"⚠️ Cache write error: {e}")
    
    return results

######################################################################
# performce optimized fetch preprint function
def fetch_preprints_with_progress(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """
    Ultra-optimized preprint search with multiple performance improvements
    """
    # Choose the best method based on expected dataset size
    try:
        # Quick size estimation
        initial_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/0"
        response = requests.get(initial_url, timeout=10)
        response.raise_for_status()
        initial_data = response.json()
        
        messages = initial_data.get("messages", [])
        estimated_total = int(messages[0]["total"]) if messages and "total" in messages[0] else 0
        
        print(f"📊 {server}: Estimated {estimated_total} articles to process")
        
        # Choose optimization strategy based on dataset size
        if estimated_total > 5000:
            print("🚀 Large dataset detected - using parallel processing")
            return fetch_preprints_parallel_pages(server, start_date, end_date, keywords, max_results, progress_callback)
        elif estimated_total > 1000:
            print("⚡ Medium dataset - using smart filtering")
            return fetch_preprints_smart_filter(server, start_date, end_date, keywords, max_results, progress_callback)
        else:
            print("📝 Small dataset - using cached approach")
            return fetch_preprints_with_cache(server, start_date, end_date, keywords, max_results, progress_callback)
            
    except Exception as e:
        print(f"⚠️ Falling back to original method due to error: {e}")
        # Fall back to your existing implementation if there are issues
        return fetch_preprints(server, start_date, end_date, keywords, max_results)


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
# add a function to merge pubmed and preprints output then highlight keywords
def merge_and_highlight_articles(articles, fields_to_highlight=["Title", "Abstract"], raw_keywords=""):
    """Highlights keywords in articles' titles and abstracts."""
    def highlight(text, terms):
        for term in sorted(terms, key=len, reverse=True):  # longest terms first
            pattern = re.compile(re.escape(term), re.IGNORECASE)
            text = pattern.sub(lambda m: f"**{m.group(0)}**", text)
        return text

    if not raw_keywords:
        return articles

    # Extract keywords while preserving operators
    terms = re.findall(r'"[^"]+"|\b\w+\*?\??\b', raw_keywords)
    terms = [term.strip('"') for term in terms if term.upper() not in ["AND", "OR", "NOT"]]

    for article in articles:
        for field in fields_to_highlight:
            if field in article and isinstance(article[field], str):
                article[field] = highlight(article[field], terms)

    return articles


######################################################################
# function to standardize date format for pubmed and preprint entries

def standardize_date_format(articles):
    """
    Standardize date formats across all articles to YYYY-MM-DD format.
    """
    from datetime import datetime
    
    for article in articles:
        if 'Publication Date' in article and article['Publication Date']:
            original_date = str(article['Publication Date'])
            
            try:
                # Handle PubMed format: 2025-Jun-25
                if '-' in original_date and any(month in original_date for month in 
                    ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
                    parsed_date = datetime.strptime(original_date, '%Y-%b-%d')
                    article['Publication Date'] = parsed_date.strftime('%Y-%m-%d')
                
                # Handle bioRxiv format: 2025/6/24 or 2025/06/24
                elif '/' in original_date:
                    # Split and pad with zeros if needed
                    parts = original_date.split('/')
                    if len(parts) == 3:
                        year, month, day = parts
                        # Pad month and day with leading zeros
                        month = month.zfill(2)
                        day = day.zfill(2)
                        standardized = f"{year}-{month}-{day}"
                        # Validate the date
                        datetime.strptime(standardized, '%Y-%m-%d')
                        article['Publication Date'] = standardized
                
                # Handle already standard format: 2025-01-01
                elif '-' in original_date and len(original_date) == 10:
                    # Validate it's in correct format
                    datetime.strptime(original_date, '%Y-%m-%d')
                    # Already correct, no change needed
                    
            except ValueError as e:
                # If any parsing fails, keep the original date
                print(f"Warning: Could not parse date '{original_date}': {e}")
                pass
        
        # Also check for 'Date' field (backup)
        elif 'Date' in article and article['Date']:
            original_date = str(article['Date'])
            
            try:
                # Same logic for 'Date' field
                if '-' in original_date and any(month in original_date for month in 
                    ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
                    parsed_date = datetime.strptime(original_date, '%Y-%b-%d')
                    article['Date'] = parsed_date.strftime('%Y-%m-%d')
                
                elif '/' in original_date:
                    parts = original_date.split('/')
                    if len(parts) == 3:
                        year, month, day = parts
                        month = month.zfill(2)
                        day = day.zfill(2)
                        standardized = f"{year}-{month}-{day}"
                        datetime.strptime(standardized, '%Y-%m-%d')
                        article['Date'] = standardized
                
                elif '-' in original_date and len(original_date) == 10:
                    datetime.strptime(original_date, '%Y-%m-%d')
                    
            except ValueError as e:
                print(f"Warning: Could not parse date '{original_date}': {e}")
                pass
    
    return articles

######################################################################
# function to standardize DOI format for pubmed and preprint entries

def standardize_doi_format(articles):
    """
    Standardize DOI formats across all articles to just the DOI identifier (without URL prefix).
    """
    
    for article in articles:
        # Check for DOI field
        if 'DOI' in article and article['DOI']:
            original_doi = str(article['DOI']).strip()
            
            # Remove common DOI URL prefixes
            prefixes_to_remove = [
                'https://doi.org/',
                'http://doi.org/',
                'https://dx.doi.org/',
                'http://dx.doi.org/',
                'doi:',
                'DOI:'
            ]
            
            standardized_doi = original_doi
            for prefix in prefixes_to_remove:
                if standardized_doi.startswith(prefix):
                    standardized_doi = standardized_doi[len(prefix):]
                    break
            
            # Update the article with standardized DOI
            article['DOI'] = standardized_doi
    
    return articles


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


######################################################################
# Implement Concurrent Processing
# Add to tracking_main.py
import concurrent.futures
import asyncio
from concurrent.futures import ThreadPoolExecutor
import threading
import time, random, logging
from datetime import timedelta


# class OptimizedPubMedFetcher:
#     def __init__(self, rate_limiter):
#         self.rate_limiter = rate_limiter
#         self.fetch_semaphore = threading.Semaphore(3)  # Limit concurrent requests
        
#     def fetch_articles_batch(self, pmid_batch, progress_callback=None):
#         """Fetch multiple articles in a single request"""
#         try:
#             with self.fetch_semaphore:
#                 # Convert pmid_batch to comma-separated string if it's a list
#                 if isinstance(pmid_batch, list):
#                     pmid_str = ",".join(pmid_batch)
#                 else:
#                     pmid_str = pmid_batch
                
#                 # Add delay for rate limiting
#                 time.sleep(random.uniform(0.5, 1.0))
                
#                 # Fetch articles using Entrez
#                 handle = Entrez.efetch(db="pubmed", id=pmid_str, rettype="xml")
#                 paper_info_list = Entrez.read(handle)
#                 handle.close()
                
#                 papers = []
                
#                 for paper_info in paper_info_list.get("PubmedArticle", []):
#                     try:
#                         article_data = paper_info["MedlineCitation"]["Article"]
#                         title = article_data["ArticleTitle"]
                        
#                         # Handle abstract
#                         abstract_sections = article_data.get("Abstract", {}).get("AbstractText", [])
#                         if abstract_sections:
#                             if isinstance(abstract_sections[0], dict):
#                                 abstract = " ".join([section.get("text", "") for section in abstract_sections])
#                             else:
#                                 abstract = str(abstract_sections[0])
#                         else:
#                             abstract = "No abstract available"
                        
#                         # FIXED: Handle publication date properly
#                         pub_date_str = "Unknown Date"
                        
#                         # Try multiple sources for publication date
#                         # 1. Try ArticleDate first (electronic publication date)
#                         article_date_list = article_data.get("ArticleDate", [])
#                         if article_date_list:
#                             try:
#                                 article_date = article_date_list[0]  # Get the first date
#                                 year = article_date.get("Year", "")
#                                 month = article_date.get("Month", "").zfill(2)
#                                 day = article_date.get("Day", "").zfill(2)
#                                 if year and month and day:
#                                     pub_date_str = f"{year}-{month}-{day}"
#                             except:
#                                 pass
                        
#                         # 2. If ArticleDate not available, try Journal PubDate
#                         if pub_date_str == "Unknown Date":
#                             try:
#                                 journal_info = article_data.get("Journal", {})
#                                 journal_issue = journal_info.get("JournalIssue", {})
#                                 pub_date = journal_issue.get("PubDate", {})
                                
#                                 year = pub_date.get("Year", "")
#                                 month = pub_date.get("Month", "")
#                                 day = pub_date.get("Day", "")
                                
#                                 # Convert month name to number if needed
#                                 if month and not month.isdigit():
#                                     month_map = {
#                                         'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
#                                         'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
#                                         'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12',
#                                         'January': '01', 'February': '02', 'March': '03', 'April': '04',
#                                         'May': '05', 'June': '06', 'July': '07', 'August': '08',
#                                         'September': '09', 'October': '10', 'November': '11', 'December': '12'
#                                     }
#                                     month = month_map.get(month, month)
                                
#                                 if year:
#                                     if month and day:
#                                         pub_date_str = f"{year}-{month.zfill(2)}-{day.zfill(2)}"
#                                     elif month:
#                                         pub_date_str = f"{year}-{month.zfill(2)}-01"
#                                     else:
#                                         pub_date_str = f"{year}-01-01"
#                             except:
#                                 pass
                        
#                         # 3. If still no date, try PubMedPubDate
#                         if pub_date_str == "Unknown Date":
#                             try:
#                                 pubmed_data = paper_info.get("PubmedData", {})
#                                 history = pubmed_data.get("History", {})
#                                 pub_med_pub_date = history.get("PubMedPubDate", [])
                                
#                                 for date_entry in pub_med_pub_date:
#                                     if date_entry.get("PubStatus") in ["pubmed", "entrez"]:
#                                         year = date_entry.get("Year", "")
#                                         month = date_entry.get("Month", "").zfill(2)
#                                         day = date_entry.get("Day", "").zfill(2)
#                                         if year and month and day:
#                                             pub_date_str = f"{year}-{month}-{day}"
#                                             break
#                             except:
#                                 pass
                        
#                         # Handle DOI
#                         doi = None
#                         article_ids = paper_info.get("PubmedData", {}).get("ArticleIdList", [])
#                         for id_item in article_ids:
#                             if hasattr(id_item, 'attributes') and id_item.attributes.get("IdType") == "doi":
#                                 doi = str(id_item)
#                                 break
                        
#                         papers.append({
#                             "Publication Date": pub_date_str,
#                             "Title": str(title),
#                             "Abstract": abstract,
#                             "DOI": f"https://doi.org/{doi}" if doi else "No DOI available"
#                         })
                        
#                         if progress_callback:
#                             progress_callback(1)  # Increment progress
                            
#                     except Exception as e:
#                         print(f"⚠️ Error parsing individual article: {e}")
#                         continue
                
#                 return papers
                
#         except Exception as e:
#             print(f"⚠️ Error fetching batch: {e}")
#             return []
    
#     def fetch_pubmed_articles_optimized(self, journal, start_date, end_date, keywords=None, progress_callback=None):
#         """Optimized version with batch processing and threading"""
#         try:
#             # Use your existing setup code
#             journal_dict = load_pubmed_journal_abbreviations()
#             formatted_journal = format_journal_abbreviation(journal, journal_dict, self.rate_limiter)
            
#             # Build query and get IDs using your existing function
#             query = build_pubmed_query(formatted_journal, start_date, end_date, keywords)
#             pmid_list, count = fetch_article_ids_from_pubmed(query, self.rate_limiter)
            
#             if not pmid_list:
#                 return []
            
#             # Split into batches of 50 (PubMed's recommended batch size)
#             batch_size = 50
#             pmid_batches = [pmid_list[i:i + batch_size] for i in range(0, len(pmid_list), batch_size)]
            
#             all_papers = []
            
#             # Process batches with limited concurrency
#             with ThreadPoolExecutor(max_workers=2) as executor:
#                 future_to_batch = {
#                     executor.submit(self.fetch_articles_batch, batch, progress_callback): batch 
#                     for batch in pmid_batches
#                 }
                
#                 for future in concurrent.futures.as_completed(future_to_batch):
#                     try:
#                         papers = future.result()
#                         all_papers.extend(papers)
#                     except Exception as e:
#                         print(f"⚠️ Batch processing error: {e}")
            
#             return all_papers
            
#         except Exception as e:
#             print(f"⚠️ Error in optimized fetch: {e}")
#             return []

class OptimizedPubMedFetcher:
    def __init__(self, rate_limiter):
        self.rate_limiter = rate_limiter  # Your PubMedRateLimit instance
        # Use rate limiter's configuration for max workers
        max_workers = 5 if rate_limiter.has_api_key else 2
        self.fetch_semaphore = threading.Semaphore(max_workers)
        self.max_workers = max_workers
        
    def fetch_articles_batch(self, pmid_batch, progress_callback=None):
        """Fetch multiple articles in a single request using the rate limiter"""
        try:
            with self.fetch_semaphore:
                # Convert pmid_batch to list if it's a string
                if isinstance(pmid_batch, str):
                    pmid_list = pmid_batch.split(',')
                else:
                    pmid_list = pmid_batch
                
                # Use rate limiter's thread-safe fetching
                self.rate_limiter._wait_for_rate_limit()
                
                # Optimized delay based on API key availability
                delay = random.uniform(0.05, 0.1) if self.rate_limiter.has_api_key else random.uniform(0.5, 1.0)
                time.sleep(delay)
                
                # Convert back to comma-separated string for Entrez
                pmid_str = ",".join(pmid_list)
                
                # Fetch articles using Entrez (will automatically use API key if configured)
                handle = Entrez.efetch(db="pubmed", id=pmid_str, rettype="xml")
                paper_info_list = Entrez.read(handle)
                handle.close()
                
                papers = []
                
                for paper_info in paper_info_list.get("PubmedArticle", []):
                    try:
                        article_data = paper_info["MedlineCitation"]["Article"]
                        title = article_data["ArticleTitle"]
                        
                        # Enhanced abstract handling
                        abstract_sections = article_data.get("Abstract", {}).get("AbstractText", [])
                        if abstract_sections:
                            if isinstance(abstract_sections[0], dict):
                                abstract = " ".join([section.get("text", "") for section in abstract_sections])
                            else:
                                abstract = str(abstract_sections[0])
                        else:
                            abstract = "No abstract available"
                        
                        # Enhanced author extraction
                        authors = []
                        author_list = article_data.get("AuthorList", [])
                        for author in author_list:
                            last_name = author.get("LastName", "")
                            first_name = author.get("ForeName", "")
                            if last_name:
                                authors.append(f"{last_name}, {first_name}")
                        
                        # Journal information
                        journal_info = article_data["Journal"]
                        journal_title = journal_info.get("Title", "")
                        journal_abbrev = journal_info.get("ISOAbbreviation", "")
                        
                        # Issue information
                        journal_issue = journal_info.get("JournalIssue", {})
                        volume = journal_issue.get("Volume", "")
                        issue = journal_issue.get("Issue", "")
                        
                        # Pages
                        pagination = article_data.get("Pagination", {})
                        pages = pagination.get("MedlinePgn", "")
                        
                        # Enhanced publication date handling (keep your existing logic)
                        pub_date_str = "Unknown Date"
                        
                        # Try multiple sources for publication date
                        # 1. Try ArticleDate first (electronic publication date)
                        article_date_list = article_data.get("ArticleDate", [])
                        if article_date_list:
                            try:
                                article_date = article_date_list[0]  # Get the first date
                                year = article_date.get("Year", "")
                                month = article_date.get("Month", "").zfill(2)
                                day = article_date.get("Day", "").zfill(2)
                                if year and month and day:
                                    pub_date_str = f"{year}-{month}-{day}"
                            except:
                                pass
                        
                        # 2. If ArticleDate not available, try Journal PubDate
                        if pub_date_str == "Unknown Date":
                            try:
                                journal_info = article_data.get("Journal", {})
                                journal_issue = journal_info.get("JournalIssue", {})
                                pub_date = journal_issue.get("PubDate", {})
                                
                                year = pub_date.get("Year", "")
                                month = pub_date.get("Month", "")
                                day = pub_date.get("Day", "")
                                
                                # Convert month name to number if needed
                                if month and not month.isdigit():
                                    month_map = {
                                        'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
                                        'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
                                        'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12',
                                        'January': '01', 'February': '02', 'March': '03', 'April': '04',
                                        'May': '05', 'June': '06', 'July': '07', 'August': '08',
                                        'September': '09', 'October': '10', 'November': '11', 'December': '12'
                                    }
                                    month = month_map.get(month, month)
                                
                                if year:
                                    if month and day:
                                        pub_date_str = f"{year}-{month.zfill(2)}-{day.zfill(2)}"
                                    elif month:
                                        pub_date_str = f"{year}-{month.zfill(2)}-01"
                                    else:
                                        pub_date_str = f"{year}-01-01"
                            except:
                                pass
                        
                        # 3. If still no date, try PubMedPubDate
                        if pub_date_str == "Unknown Date":
                            try:
                                pubmed_data = paper_info.get("PubmedData", {})
                                history = pubmed_data.get("History", {})
                                pub_med_pub_date = history.get("PubMedPubDate", [])
                                
                                for date_entry in pub_med_pub_date:
                                    if date_entry.get("PubStatus") in ["pubmed", "entrez"]:
                                        year = date_entry.get("Year", "")
                                        month = date_entry.get("Month", "").zfill(2)
                                        day = date_entry.get("Day", "").zfill(2)
                                        if year and month and day:
                                            pub_date_str = f"{year}-{month}-{day}"
                                            break
                            except:
                                pass
                        
                        # Handle DOI
                        # Enhanced DOI handling
                        doi = None
                        article_ids = paper_info.get("PubmedData", {}).get("ArticleIdList", [])
                        for id_item in article_ids:
                            if hasattr(id_item, 'attributes') and id_item.attributes.get("IdType") == "doi":
                                doi = str(id_item)
                                break
                        
                        # PMID
                        pmid_value = paper_info["MedlineCitation"]["PMID"]
                        
                        # Extract year, month, day from your existing pub_date_str logic
                        year = month = day = ""
                        if pub_date_str != "Unknown Date":
                            try:
                                parts = pub_date_str.split('-')
                                year = parts[0] if len(parts) > 0 else ""
                                month = parts[1] if len(parts) > 1 else ""
                                day = parts[2] if len(parts) > 2 else ""
                            except:
                                pass
                        
                        papers.append({
                            "Publication Date": pub_date_str,
                            "Title": str(title),
                            "Authors": " and ".join(authors) if authors else "Unknown Author",
                            "Abstract": abstract,
                            "Journal": journal_title or journal_abbrev,
                            "Volume": volume,
                            "Issue": issue,
                            "Pages": pages,
                            "Year": year,
                            "Month": month,
                            "Day": day,
                            "DOI": f"https://doi.org/{doi}" if doi else "No DOI available",
                            "PMID": str(pmid_value)
                        })
                        
                        if progress_callback:
                            progress_callback(1)  # Increment progress
                            
                    except Exception as e:
                        print(f"⚠️ Error parsing individual article: {e}")
                        continue
                
                return papers
                
        except Exception as e:
            print(f"⚠️ Error fetching batch: {e}")
            return []
    
    def fetch_pubmed_articles_optimized(self, journal, start_date, end_date, keywords=None, progress_callback=None):
        """Optimized version with batch processing and threading - works with PubMedRateLimit"""
        try:
            # Use your existing setup code
            journal_dict = load_pubmed_journal_abbreviations()
            formatted_journal = format_journal_abbreviation(journal, journal_dict, self.rate_limiter)
            
            # Build query and get IDs using your existing function
            query = build_pubmed_query(formatted_journal, start_date, end_date, keywords)
            pmid_list, count = fetch_article_ids_from_pubmed(query, self.rate_limiter)
            
            if not pmid_list:
                return []
            
            # Optimize batch size based on API key availability
            batch_size = 200 if self.rate_limiter.has_api_key else 50
            pmid_batches = [pmid_list[i:i + batch_size] for i in range(0, len(pmid_list), batch_size)]
            
            all_papers = []
            
            # Process batches with optimized concurrency
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_batch = {
                    executor.submit(self.fetch_articles_batch, batch, progress_callback): batch 
                    for batch in pmid_batches
                }
                
                for future in concurrent.futures.as_completed(future_to_batch):
                    try:
                        papers = future.result()
                        all_papers.extend(papers)
                    except Exception as e:
                        print(f"⚠️ Batch processing error: {e}")
            
            return all_papers
            
        except Exception as e:
            print(f"⚠️ Error in optimized fetch: {e}")
            return []
    
    def fetch_articles_with_progress(self, journal_list, start_date, end_date, keywords=None, progress_callback=None):
        """
        Fetch articles from multiple journals with detailed progress tracking
        """
        all_articles = []
        
        for i, journal in enumerate(journal_list):
            try:
                print(f"🔍 Processing journal {i+1}/{len(journal_list)}: {journal}")
                
                # Create journal-specific progress callback
                def journal_progress_callback(count):
                    if progress_callback:
                        progress_callback(count, journal)
                
                # Fetch articles for this journal
                articles = self.fetch_pubmed_articles_optimized(
                    journal, start_date, end_date, keywords, journal_progress_callback
                )
                
                # Add journal information to each article
                for article in articles:
                    article["Journal"] = journal
                    article["Source"] = "PubMed"
                
                all_articles.extend(articles)
                
                print(f"✅ {journal}: {len(articles)} articles found")
                
            except Exception as e:
                print(f"❌ Error processing {journal}: {e}")
                continue
        
        return all_articles
    
    def get_performance_stats(self):
        """Get performance statistics"""
        rate_info = self.rate_limiter.get_rate_limit_info()
        return {
            "api_key_enabled": rate_info["has_api_key"],
            "max_concurrent_requests": self.max_workers,
            "batch_size": 200 if rate_info["has_api_key"] else 50,
            "total_api_requests": rate_info["total_requests"],
            "rate_limit_per_second": rate_info["max_requests_per_second"]
        }

def execute_subscription_search(journals, keywords, include_preprints, frequency):
    """Execute search based on subscription parameters"""
    from datetime import datetime, timedelta
    
    # Calculate date range based on frequency
    today = datetime.today().date()
    
    if frequency == "weekly":
        days_back = 7
    elif frequency == "monthly":
        days_back = 30
    elif frequency.startswith("every"):
        # Extract number from "every X days"
        try:
            days_back = int(frequency.split()[1])
        except:
            days_back = 7  # fallback
    else:
        days_back = 7  # Default fallback
    
    start_date = str(today - timedelta(days=days_back))
    end_date = str(today)
    
    print(f"🔍 Searching from {start_date} to {end_date}")
    
    all_articles = []
    
    try:
        # Search PubMed journals
        if journals:
            for journal in journals:
                try:
                    articles = fetch_pubmed_articles_by_date(
                        journal, start_date, end_date, 
                        format_boolean_keywords_for_pubmed(keywords) if keywords else None
                    )
                    for article in articles:
                        article["Journal"] = journal
                        article["Source"] = "PubMed"
                    all_articles.extend(articles)
                    print(f"✅ PubMed - {journal}: {len(articles)} articles")
                except Exception as e:
                    print(f"❌ Error searching PubMed {journal}: {e}")
                    continue
        
        # Search preprints
        if include_preprints:
            for server in ["biorxiv", "medrxiv"]:
                try:
                    preprints = fetch_preprints(
                        server=server,
                        start_date=start_date,
                        end_date=end_date,
                        keywords=keywords  # Pass raw keywords, not formatted
                    )
                    for article in preprints:
                        article["Journal"] = server
                        article["Source"] = "Preprint"
                    all_articles.extend(preprints)
                    print(f"✅ Preprints - {server}: {len(preprints)} articles")
                except Exception as e:
                    print(f"❌ Error searching preprints {server}: {e}")
                    continue
        
        if all_articles:
            # Process results
            all_articles = standardize_date_format(all_articles)
            all_articles = standardize_doi_format(all_articles)
            merged = merge_and_highlight_articles(all_articles, [], keywords)
            return pd.DataFrame(merged)
        else:
            return pd.DataFrame()
            
    except Exception as e:
        print(f"❌ Error in subscription search: {e}")
        return pd.DataFrame()

def test_api_configuration_with_rate_limiter():
    """Test if API key and email are properly configured using PubMedRateLimit"""
    try:
        print("🔍 Testing NCBI API configuration with rate limiter...")
        
        # Initialize your rate limiter (this configures Entrez)
        rate_limiter = PubMedRateLimit()
        
        # Test with a simple search using the rate limiter
        search_result = rate_limiter.safe_pubmed_search("cancer[Title]", max_results=5)
        
        if search_result and search_result.get("IdList"):
            print(f"✅ API test successful! Found {len(search_result['IdList'])} articles")
            print(f"📧 Email: {Entrez.email}")
            print(f"🔑 API Key: {'*' * (len(Entrez.api_key) - 4) + Entrez.api_key[-4:] if Entrez.api_key else 'Not configured'}")
            
            # Show rate limiter configuration
            config = rate_limiter.get_rate_limit_info()
            print(f"🚀 Rate Limits: {config['max_requests_per_second']} req/sec, {config['max_requests_per_minute']} req/min")
            print(f"⚡ Performance Mode: {'High (API Key)' if config['has_api_key'] else 'Standard (No API Key)'}")
            return True
        else:
            print("❌ API test failed - no results returned")
            return False
            
    except Exception as e:
        print(f"❌ API test failed: {e}")
        return False


#####################################
# generate BibTex format

# def generate_bibtex_from_dataframe(df):
#     """Generate BibTeX format from search results DataFrame"""
#     bibtex_entries = []
    
#     for index, row in df.iterrows():
#         # Clean up the title (remove markdown formatting)
#         title = str(row.get('Title', 'Unknown Title')).replace('**', '')
        
#         # Clean up author field if it exists
#         authors = row.get('Authors', 'Unknown Author')
#         if pd.isna(authors) or authors == 'Unknown Author':
#             authors = 'Unknown Author'
        
#         # Get journal name
#         journal = str(row.get('Journal', 'Unknown Journal'))
        
#         # Get year from publication date
#         pub_date = str(row.get('Publication Date', ''))
#         year = 'Unknown Year'
#         if pub_date and pub_date != 'nan':
#             try:
#                 year = pub_date.split('-')[0] if '-' in pub_date else pub_date[:4]
#             except:
#                 year = 'Unknown Year'
        
#         # Get DOI
#         doi = str(row.get('DOI', ''))
#         doi_clean = doi.replace('https://doi.org/', '') if doi and doi != 'nan' and doi != 'N/A' else ''
        
#         # Create entry key (first author + year + first word of title)
#         first_word = title.split()[0] if title else 'Unknown'
#         entry_key = f"{authors.split(',')[0].replace(' ', '')}_{year}_{first_word}".replace(' ', '').replace('.', '')
        
#         # Create BibTeX entry
#         bibtex_entry = f"""@article{{{entry_key},
#     title = {{{title}}},
#     author = {{{authors}}},
#     journal = {{{journal}}},
#     year = {{{year}}}"""
        
#         if doi_clean:
#             bibtex_entry += f",\n    doi = {{{doi_clean}}}"
        
#         if doi:
#             bibtex_entry += f",\n    url = {{{doi}}}"
            
#         bibtex_entry += "\n}\n"
#         bibtex_entries.append(bibtex_entry)
    
#     return "\n".join(bibtex_entries)


def generate_bibtex_from_dataframe(df):
    """Generate comprehensive BibTeX format from search results DataFrame"""
    bibtex_entries = []
    
    for index, row in df.iterrows():
        # Clean up the title (remove markdown formatting)
        title = str(row.get('Title', 'Unknown Title')).replace('**', '')
        
        # Authors
        authors = row.get('Authors', 'Unknown Author')
        if pd.isna(authors) or authors == 'Unknown Author':
            authors = 'Unknown Author'
        
        # Journal information
        journal = str(row.get('Journal', 'Unknown Journal'))
        
        # Publication details
        year = str(row.get('Year', ''))
        if not year or year == 'nan':
            pub_date = str(row.get('Publication Date', ''))
            if pub_date and pub_date != 'nan':
                try:
                    year = pub_date.split('-')[0] if '-' in pub_date else pub_date[:4]
                except:
                    year = 'Unknown Year'
        
        volume = str(row.get('Volume', ''))
        issue = str(row.get('Issue', ''))
        pages = str(row.get('Pages', ''))
        month = str(row.get('Month', ''))
        
        # DOI and PMID
        doi = str(row.get('DOI', ''))
        pmid = str(row.get('PMID', ''))
        
        # Abstract
        abstract = str(row.get('Abstract', ''))
        
        # Create entry key
        first_author = authors.split(' and ')[0].split(',')[0].replace(' ', '') if authors != 'Unknown Author' else 'UnknownAuthor'
        first_word = title.split()[0] if title else 'Unknown'
        entry_key = f"{first_author}_{year}_{first_word}".replace(' ', '').replace('.', '')
        
        # Create BibTeX entry
        bibtex_entry = f"""@article{{{entry_key},
    title = {{{title}}},
    author = {{{authors}}},
    journal = {{{journal}}},
    year = {{{year}}}"""
        
        # Add optional fields if available
        if volume and volume != 'nan':
            bibtex_entry += f",\n    volume = {{{volume}}}"
        
        if issue and issue != 'nan':
            bibtex_entry += f",\n    number = {{{issue}}}"
            
        if pages and pages != 'nan':
            bibtex_entry += f",\n    pages = {{{pages}}}"
            
        if month and month != 'nan':
            bibtex_entry += f",\n    month = {{{month}}}"
        
        if doi and doi != 'nan' and doi != 'N/A':
            doi_clean = doi.replace('https://doi.org/', '')
            bibtex_entry += f",\n    doi = {{{doi_clean}}}"
        
        if pmid and pmid != 'nan':
            bibtex_entry += f",\n    pmid = {{{pmid}}}"
        
        if doi and doi != 'nan' and doi != 'N/A':
            bibtex_entry += f",\n    url = {{{doi}}}"
        
        # Add abstract if available and not too long
        if abstract and abstract != 'nan' and len(abstract) > 10:
            # Truncate very long abstracts
            if len(abstract) > 500:
                abstract = abstract[:500] + "..."
            bibtex_entry += f",\n    abstract = {{{abstract}}}"
            
        bibtex_entry += "\n}\n"
        bibtex_entries.append(bibtex_entry)
    
    return "\n".join(bibtex_entries)