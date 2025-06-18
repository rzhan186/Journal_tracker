import os
import pandas as pd
import re
import requests
import json

from Bio import Entrez
from datetime import datetime
from calendar import monthrange

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



# Now I want that when I run

# def fetch_pubmed_articles_by_date(journal, start_month=None, end_month=None):

# The user can enter either the full name or the abbreviation of the journal stored in journal_dict.
# The function will halt if neither full name or the abbreviation is found in the journal_dict. 

# If full name is found, the function will automatically select its abbreviation and format it to be used in the query.
# If abbreviation is used, the function will format it and be used in the query. 



def format_journal_abbreviation(user_input, journal_dict):
    """
    Ensures proper abbreviation formatting for PubMed searches:
    - If abbreviation has multiple words, ensure it ends with a dot.
    - If abbreviation is a single word, leave it as is.
    """
    # Normalize user input (remove accidental extra spaces)
    user_input = user_input.strip()

    # Check if input is a known abbreviation
    if user_input in journal_dict:
        words = user_input.split()
        # If multiple words, ensure a dot at the end
        formatted_abbr = user_input + "." if len(words) > 1 and not user_input.endswith(".") else user_input
        print(f"✅ Formatted journal abbreviation: {formatted_abbr}")
        return formatted_abbr

    #print(f"⚠️'{user_input}' not found in PubMed journal list. Searching as entered.")
    return user_input  # Use raw input if not found






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

journal_dict = load_pubmed_journal_abbreviations()


def format_journal_abbreviation(journal, journal_dict):
    """
    Validates and formats a journal input for a PubMed query.
    - If the user enters the full name, it retrieves its abbreviation.
    - If the user enters an abbreviation, it ensures proper formatting.
    - If neither is found, the function halts execution.
    """
    # Check if input is a full journal name (value in the dictionary)
    if journal in journal_dict.values():
        for abbr, full_name in journal_dict.items():
            if full_name == journal:
                formatted_abbr = abbr + "." if " " in abbr and not abbr.endswith(".") else abbr
                print(f"✅ Found full name. Using abbreviation: {formatted_abbr}")
                return formatted_abbr

    # Check if input is an abbreviation (key in the dictionary)
    elif journal in journal_dict:
        formatted_abbr = journal + "." if " " in journal and not journal.endswith(".") else journal
        print(f"✅ Found abbreviation. Using formatted abbreviation: {formatted_abbr}")
        return formatted_abbr

    # Halt execution if the journal is not found
    print(f"❌ Error: '{journal}' not found in PubMed journal list.")
    exit(1)  # Terminate execution



def format_journal_abbreviation(journal, journal_dict):
    """
    Optimized function to validate and format a journal name or abbreviation for PubMed.
    - Uses dictionary lookups for O(1) efficiency.
    - If a full name is given, retrieves and formats its abbreviation.
    - If an abbreviation is given, ensures proper formatting.
    - If neither is found, halts immediately.
    """
    # Create a reverse lookup dictionary for full names → abbreviations
    full_to_abbr = {full_name: abbr for abbr, full_name in journal_dict.items()}

    # ✅ Fast check: If input is a known abbreviation
    if journal in journal_dict:
        formatted_abbr = journal + "." if " " in journal and not journal.endswith(".") else journal
        print(f"✅ Found abbreviation. Using formatted abbreviation: {formatted_abbr}")
        return formatted_abbr

    # ✅ Fast check: If input is a full journal name
    elif journal in full_to_abbr:
        abbr = full_to_abbr[journal]
        formatted_abbr = abbr + "." if " " in abbr and not abbr.endswith(".") else abbr
        print(f"✅ Found full name. Using abbreviation: {formatted_abbr}")
        return formatted_abbr

    # ❌ Journal not found → Exit immediately
    print(f"❌ Error: '{journal}' not found in PubMed journal list.")
    exit(1)



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


journal="Limnology and"
journal="Nature"

formatted_journal = format_journal_abbreviation(journal, journal_dict)
