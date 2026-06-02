import os
import re
import json
import time
import random
import difflib
import concurrent.futures
import threading
from calendar import monthrange
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

from Bio import Entrez
import pandas as pd
import requests
from dotenv import load_dotenv

from RateLimit import PubMedRateLimit

load_dotenv()

Entrez.email = os.getenv("PUBMED_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

if not Entrez.email:
    print("⚠️ Warning: PUBMED_EMAIL not found in environment variables")
if not Entrez.api_key:
    print("⚠️ Warning: NCBI_API_KEY not found in environment variables")
else:
    print("✅ NCBI API key configured successfully")


def validate_date_input(date_str):
    """Validates user input for YYYY-MM or YYYY-MM-DD format."""
    pattern = r"^\d{4}-(0[1-9]|1[0-2])(-([0-2][0-9]|3[01]))?$"
    return bool(re.match(pattern, date_str))


def validate_dates(start_date, end_date):
    """Validates the start and end date formats and order."""
    today = datetime.today()

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


def get_last_day_of_month(year, month):
    """Returns the last valid day of a given month."""
    return monthrange(int(year), int(month))[1]


JOURNAL_LIST_URL = "https://ftp.ncbi.nih.gov/pubmed/J_Medline.txt"
FALLBACK_JOURNAL_LIST_URL = "https://raw.githubusercontent.com/rzhan186/Journal_tracker/refs/heads/main/src/pubmed_journal_abbreviations.json"
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
        response.raise_for_status()
        print("✅ Successfully fetched from primary URL.")

    except (requests.HTTPError, requests.ConnectionError) as e:
        print(f"⚠️ Primary URL failed: {e}. Trying fallback URL...")
        try:
            response = requests.get(FALLBACK_JOURNAL_LIST_URL)
            response.raise_for_status()
            print("✅ Successfully fetched from fallback URL.")

        except (requests.HTTPError, requests.ConnectionError) as e:
            print(f"⚠️ Fallback URL also failed: {e}. Unable to fetch journal abbreviations.")
            return {}

    text_data = response.text.split("\n")

    journal_dict = {}
    current_abbr = None
    current_full = None

    for line in text_data:
        line = line.strip()

        if not line or line.startswith("-"):
            continue

        if line.startswith("MedAbbr: "):
            current_abbr = line.replace("MedAbbr: ", "").strip()
        elif line.startswith("JournalTitle: "):
            current_full = line.replace("JournalTitle: ", "").strip()

        if current_abbr and current_full:
            journal_dict[current_abbr] = current_full
            current_abbr, current_full = None, None

    with open(JOURNAL_LIST_FILE, "w", encoding="utf-8") as file:
        json.dump(journal_dict, file, ensure_ascii=False, indent=4)

    print(f"✅ Journal abbreviations saved to {JOURNAL_LIST_FILE}")
    return journal_dict


def load_pubmed_journal_abbreviations():
    """
    Loads journal abbreviations from a local JSON file.
    If the file does not exist, fetches and saves it first.
    """
    try:
        with open(JOURNAL_LIST_FILE, "r", encoding="utf-8") as file:
            journal_dict = json.load(file)
        print(f"📂 Loaded journal abbreviations from {JOURNAL_LIST_FILE}")
        return journal_dict
    except FileNotFoundError:
        print("⚠️ No local file found. Fetching data from PubMed...")
        return get_pubmed_journal_abbreviations()


def format_journal_abbreviation(journal, journal_dict, rate_limiter=None):
    """
    Validates and formats a journal name or abbreviation for PubMed.
    Checks manual mappings, the official PubMed list, and tests directly with PubMed.
    Provides fuzzy suggestions for unrecognized names.
    """
    manual_mappings = {
        "Lancet": "Lancet",
        "The Lancet": "Lancet",
        "BMJ": "BMJ"
    }

    if journal in manual_mappings:
        formatted_abbr = manual_mappings[journal]
        print(f"✅ Found in manual mappings. Using: {formatted_abbr}")
        return formatted_abbr

    abbreviations_set = set(journal_dict.keys())
    full_names_set = set(journal_dict.values())

    if journal in abbreviations_set:
        print(f"✅ Found abbreviation in official list. Using: {journal}")
        return journal

    if journal in full_names_set:
        abbr = [k for k, v in journal_dict.items() if v == journal][0]
        print(f"✅ Found full name in official list. Using abbreviation: {abbr}")
        return abbr

    print(f"🔍 Testing '{journal}' directly with PubMed...")

    if test_journal_name_on_pubmed(journal, rate_limiter):
        print(f"✅ '{journal}' works directly with PubMed!")
        return journal

    if test_journal_name_on_pubmed(journal):
        print(f"✅ '{journal}' works directly with PubMed!")
        return journal

    all_possible_names = list(full_names_set.union(abbreviations_set).union(manual_mappings.keys()))
    possible_matches = difflib.get_close_matches(journal, all_possible_names, n=5, cutoff=0.6)

    suggestion_msg = ""
    if possible_matches:
        suggestion_msg = "\n🛈 Did you mean: " + ", ".join(f"'{s}'" for s in possible_matches)

    raise ValueError(f"❌ Error: '{journal}' not found in any source.{suggestion_msg}")


def test_journal_name_on_pubmed(journal_name, rate_limiter=None):
    """Test if a journal name works directly with PubMed by doing a small search."""
    try:
        test_query = f'"{journal_name}"[Journal] AND ("2023/01/01"[Date - Publication] : "2024/12/31"[Date - Publication])'

        if rate_limiter:
            record = rate_limiter.safe_pubmed_search(test_query, max_results=1)
        else:
            handle = Entrez.esearch(db="pubmed", term=test_query, retmax=1)
            record = Entrez.read(handle)
            handle.close()

        return record and len(record.get("IdList", [])) > 0

    except Exception:
        return False


def fetch_article_ids_from_pubmed(query, rate_limiter=None):
    """Fetch article IDs from PubMed based on the provided query."""
    if rate_limiter:
        record = rate_limiter.safe_pubmed_search(query, max_results=1000)
        if record:
            return record["IdList"], len(record["IdList"])
        else:
            return [], 0
    else:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="pub date")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"], len(record["IdList"])


def format_boolean_keywords_for_pubmed(raw_query):
    """
    Formats a Boolean keyword query for PubMed Title/Abstract search.
    Preserves wildcards (*, ?), handles phrases, and attaches [Title/Abstract].
    """
    tokens = re.findall(r'"[^"]+"|\(|\)|\bAND\b|\bOR\b|\bNOT\b|\*|\w[\w*?\-]*', raw_query, flags=re.IGNORECASE)
    formatted_tokens = []

    for token in tokens:
        upper_token = token.upper()

        if upper_token in {"AND", "OR", "NOT", "(", ")"}:
            formatted_tokens.append(upper_token)
        elif token.startswith('"') and token.endswith('"'):
            phrase = token.strip('"')
            formatted_tokens.append(f'"{phrase}"[Title/Abstract]')
        else:
            formatted_tokens.append(f'{token}[Title/Abstract]')

    return " ".join(formatted_tokens)


def build_pubmed_query(journal, start_date, end_date, keywords=None):
    """Construct the PubMed query string, supporting Boolean keyword search."""
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

    exclude_types = [
        '"news"[Publication Type]',
        '"comment"[Publication Type]',
        '"editorial"[Publication Type]',
        '"letter"[Publication Type]'
    ]

    base_query = f'"{journal}"[Journal] AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication]) AND hasabstract[text]'
    inclusion_clause = f' AND ({" OR ".join(include_types)})'
    exclusion_clause = f' NOT ({" OR ".join(exclude_types)})'
    full_query = base_query + inclusion_clause + exclusion_clause

    if keywords:
        full_query += f' AND ({keywords})'

    return full_query


def fetch_pubmed_articles_by_date(journal, start_date=None, end_date=None, keywords=None, rate_limiter=None):
    journal_dict = load_pubmed_journal_abbreviations()
    formatted_journal = format_journal_abbreviation(journal, journal_dict, rate_limiter)

    today = datetime.today()
    current_month = today.strftime("%Y-%m")

    if not start_date:
        start_date = current_month
        end_date = current_month
        print(f"No date provided, fetching articles from {journal} for the current month {today.strftime('%Y-%m')}.")
    else:
        validate_dates(start_date, end_date or current_month)

    if len(start_date) == 7:
        start_date = f"{start_date}/01"
    if len(end_date) == 7:
        end_date = f"{end_date}/{get_last_day_of_month(end_date[:4], end_date[5:])}"

    print(f"Date range: {start_date} to {end_date}")

    query = build_pubmed_query(formatted_journal, start_date, end_date, keywords)
    print(f"\n🔍 Final PubMed Query:\n{query}\n")

    pmid_list, count = fetch_article_ids_from_pubmed(query, rate_limiter)

    if not pmid_list:
        print(f"No articles found for {formatted_journal} from {start_date} to {end_date}.")
        return []

    print(f"✅ {count} papers found.")
    papers = []

    for index, pmid in enumerate(pmid_list):
        print(f"Fetching article {index + 1}/{len(pmid_list)} (PMID: {pmid})...")

        try:
            if rate_limiter:
                paper_xml = rate_limiter.safe_pubmed_fetch([pmid])
                if not paper_xml:
                    continue
                time.sleep(1)
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            paper_info = Entrez.read(handle)
            handle.close()

            article_data = paper_info["PubmedArticle"][0]["MedlineCitation"]["Article"]
            title = article_data["ArticleTitle"]

            abstract_sections = article_data.get("Abstract", {}).get("AbstractText", [])
            if abstract_sections:
                if isinstance(abstract_sections[0], dict):
                    abstract = " ".join([s.get("text", "") for s in abstract_sections])
                else:
                    abstract = str(abstract_sections[0])
            else:
                abstract = "No abstract available"

            authors = []
            for author in article_data.get("AuthorList", []):
                last_name = author.get("LastName", "")
                first_name = author.get("ForeName", "")
                if last_name:
                    authors.append(f"{last_name}, {first_name}")

            journal_info = article_data["Journal"]
            journal_title = journal_info.get("Title", journal)
            journal_abbrev = journal_info.get("ISOAbbreviation", "")

            journal_issue = journal_info.get("JournalIssue", {})
            volume = journal_issue.get("Volume", "")
            issue = journal_issue.get("Issue", "")

            pub_date = journal_issue.get("PubDate", {})
            pub_year = pub_date.get("Year", "Unknown Year")
            pub_month = pub_date.get("Month", "Unknown Month")
            pub_day = pub_date.get("Day", "")

            pagination = article_data.get("Pagination", {})
            pages = pagination.get("MedlinePgn", "")

            doi = None
            for id_item in paper_info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]:
                if id_item.attributes["IdType"] == "doi":
                    doi = str(id_item)
                    break

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

            print(f"✅ Fetched: {title}")

        except Exception as e:
            print(f"⚠️ Failed to fetch article {pmid}: {e}")
            continue

    print("😊 All articles fetched successfully.")
    return papers


def compile_keyword_filter(raw_query):
    """
    Compiles a Boolean keyword filter for matching against article text.
    Supports AND, OR, NOT operators, quoted phrases, and wildcards (*, ?).
    """
    def tokenize(query):
        return re.findall(r'\(|\)|\bAND\b|\bOR\b|\bNOT\b|"[^"]+"|\w[\w*?]*', query, flags=re.IGNORECASE)

    def to_python_expr(tokens):
        expr_parts = []
        for i, token in enumerate(tokens):
            upper = token.upper()
            if upper in {"AND", "OR", "NOT"}:
                if upper == "NOT" and i > 0 and tokens[i-1].upper() not in {"AND", "OR", "("}:
                    expr_parts.append(" and ")
                expr_parts.append(f" {upper.lower()} ")
            elif token in "()":
                expr_parts.append(token)
            else:
                if (i > 0 and tokens[i-1].upper() not in {"AND", "OR", "NOT", "("}
                        and tokens[i-1] != "("):
                    expr_parts.append(" and ")
                clean_token = token.strip('"')
                escaped_token = re.escape(clean_token).replace(r'\*', '.*').replace(r'\?', '.')
                expr_parts.append(f'bool(re.search(r"{escaped_token}", text, re.IGNORECASE))')
        return ''.join(expr_parts)

    def matcher(text):
        try:
            tokens = tokenize(raw_query)
            expr = to_python_expr(tokens)
            if not expr.strip():
                return True
            return eval(expr, {"re": re, "text": text, "bool": bool})
        except Exception:
            return False

    return matcher


def fetch_preprints(server="biorxiv", start_date=None, end_date=None, keywords=None, max_results=None):
    """
    Fetches preprints from bioRxiv or medRxiv.

    Args:
        server: 'biorxiv' or 'medrxiv'
        start_date: Start date in YYYY-MM-DD format.
        end_date: End date in YYYY-MM-DD format.
        keywords: Boolean keyword query (AND, OR, NOT, *, ?, parentheses).
        max_results: Maximum results to fetch (None = unlimited).

    Returns:
        List[Dict] with keys Journal, Title, Abstract, Publication Date, DOI.
    """
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    matcher = compile_keyword_filter(keywords) if keywords else lambda x: True

    if not (start_date and end_date):
        print("⚠️ Start and end dates are required for bioRxiv/medRxiv.")
        return []

    try:
        if isinstance(start_date, str) and len(start_date) == 7:
            start_date = f"{start_date}-01"
        if isinstance(start_date, str):
            datetime.strptime(start_date, '%Y-%m-%d')

        if isinstance(end_date, str) and len(end_date) == 7:
            year, month = end_date.split('-')
            last_day = monthrange(int(year), int(month))[1]
            end_date = f"{end_date}-{last_day:02d}"
        if isinstance(end_date, str):
            datetime.strptime(end_date, '%Y-%m-%d')

    except ValueError as e:
        print(f"❌ Invalid date format: {e}")
        return []

    print(f"🔍 Fetching from {server}: {start_date} to {end_date}")

    results = []
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
        total_papers = int(messages[0]["total"]) if messages and "total" in messages[0] else len(initial_data.get("collection", []))
        print(f"📊 {server}: Found {total_papers} preprints.")

        cursor = 0
        page_size = 100
        all_preprints = []

        while True:
            paginated_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/{cursor}"
            try:
                response = requests.get(paginated_url, timeout=30)
                response.raise_for_status()

                if not response.headers.get('content-type', '').startswith('application/json'):
                    break

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
                        print(f"⚠️ Error processing item: {e}")
                        continue

                print(f"📄 {server}: Page {cursor//page_size + 1} — {len(collection)} papers, {page_matches} matched")

                cursor += len(collection)

                if len(collection) < page_size:
                    break

                time.sleep(0.5)

                if max_results and len(all_preprints) >= max_results:
                    all_preprints = all_preprints[:max_results]
                    break

            except requests.exceptions.Timeout:
                print(f"⏱️ Timeout fetching page from {server}")
                break
            except requests.exceptions.RequestException as e:
                print(f"⚠️ Request error: {e}")
                break
            except Exception as e:
                print(f"⚠️ Unexpected error: {e}")
                break

        results.extend(all_preprints)

    except requests.exceptions.Timeout:
        print(f"⏱️ Timeout fetching from {server}")
        return []
    except requests.exceptions.RequestException as e:
        print(f"⚠️ Request error: {e}")
        return []
    except Exception as e:
        print(f"⚠️ Error: {e}")
        return []

    print(f"✅ Fetched {len(results)} preprints from {server}")
    return results


def fetch_preprints_parallel_pages(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """Fetch preprint pages in parallel to speed up data retrieval."""
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    if keywords and keywords.strip():
        try:
            matcher = compile_keyword_filter(keywords.strip())
        except Exception as e:
            print(f"⚠️ Error compiling keyword filter: {e}")
            matcher = lambda x: True
    else:
        matcher = lambda x: True

    try:
        initial_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/0"
        response = requests.get(initial_url, timeout=30)
        response.raise_for_status()
        initial_data = response.json()

        messages = initial_data.get("messages", [])
        total_papers = int(messages[0]["total"]) if messages and "total" in messages[0] else len(initial_data.get("collection", []))

        print(f"📊 {server}: Found {total_papers} total preprints")

        page_size = 100
        total_pages = (total_papers + page_size - 1) // page_size
        max_workers = min(5, total_pages)

        def fetch_single_page(cursor):
            try:
                url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/{cursor}"
                r = requests.get(url, timeout=30)
                r.raise_for_status()

                collection = r.json().get("collection", [])
                page_results = []

                for item in collection:
                    try:
                        title = item.get("title", "").strip()
                        abstract = item.get("abstract", "").strip()

                        if not title and not abstract:
                            continue

                        if matcher(f"{title} {abstract}"):
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
                    except Exception:
                        continue

                return page_results, len(collection)

            except Exception as e:
                print(f"⚠️ Error fetching page at cursor {cursor}: {e}")
                return [], 0

        all_results = []
        total_processed = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            cursor_values = list(range(0, total_papers, page_size))
            future_to_cursor = {executor.submit(fetch_single_page, c): c for c in cursor_values}

            for future in concurrent.futures.as_completed(future_to_cursor):
                cursor = future_to_cursor[future]
                try:
                    page_results, page_count = future.result()
                    all_results.extend(page_results)
                    total_processed += page_count

                    if progress_callback:
                        progress_callback(len(all_results))

                    print(f"📄 {server}: Page {cursor//page_size + 1} — {page_count} processed, {len(page_results)} matched")

                    if max_results and len(all_results) >= max_results:
                        for f in future_to_cursor:
                            if f != future:
                                f.cancel()
                        break

                except Exception as e:
                    print(f"⚠️ Error processing page: {e}")
                    continue

        if max_results and len(all_results) > max_results:
            all_results = all_results[:max_results]

        print(f"✅ {server}: {len(all_results)} matching preprints out of {total_processed} processed")
        return all_results

    except Exception as e:
        print(f"⚠️ Error in parallel fetch from {server}: {e}")
        return []


def fetch_preprints_smart_filter(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """Preprint search with a quick pre-filter pass to skip non-matching articles early."""
    if server not in ["biorxiv", "medrxiv"]:
        print("❌ Invalid server. Choose 'biorxiv' or 'medrxiv'.")
        return []

    keyword_patterns = []
    has_filters = False
    if keywords and keywords.strip():
        try:
            terms = re.findall(r'\b\w+\b', keywords.lower())
            keyword_patterns = [re.compile(term, re.IGNORECASE) for term in terms if len(term) > 2]
            matcher = compile_keyword_filter(keywords.strip())
            has_filters = True
        except Exception as e:
            print(f"⚠️ Error compiling filters: {e}")
            matcher = lambda x: True
    else:
        matcher = lambda x: True

    def quick_filter(text):
        if not has_filters:
            return True
        text_lower = text.lower()
        return not keyword_patterns or any(p.search(text_lower) for p in keyword_patterns)

    results = []
    try:
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

        while cursor < total_papers:
            try:
                url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/{cursor}"
                response = requests.get(url, timeout=30)
                response.raise_for_status()

                collection = response.json().get("collection", [])
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
                        if quick_filter(combined_text) and matcher(combined_text):
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

                        if total_processed % 50 == 0 and progress_callback:
                            progress_callback(len(results))

                    except Exception:
                        continue

                print(f"📄 {server}: Page {cursor//page_size + 1} — {len(collection)} processed, {page_matches} matched")

                cursor += len(collection)

                if max_results and len(results) >= max_results:
                    results = results[:max_results]
                    break

                time.sleep(0.1)

            except Exception as e:
                print(f"⚠️ Error fetching page: {e}")
                break

        print(f"✅ {server}: {len(results)} matches found")
        return results

    except Exception as e:
        print(f"⚠️ Error in smart filter from {server}: {e}")
        return []


def fetch_preprints_with_progress(server, start_date, end_date, keywords=None, max_results=None, progress_callback=None):
    """Fetches preprints using the optimal strategy based on estimated dataset size."""
    try:
        initial_url = f"https://api.biorxiv.org/details/{server}/{start_date}/{end_date}/0"
        response = requests.get(initial_url, timeout=10)
        response.raise_for_status()
        initial_data = response.json()

        messages = initial_data.get("messages", [])
        estimated_total = int(messages[0]["total"]) if messages and "total" in messages[0] else 0

        print(f"📊 {server}: Estimated {estimated_total} articles to process")

        if estimated_total > 5000:
            print("🚀 Large dataset — using parallel processing")
            return fetch_preprints_parallel_pages(server, start_date, end_date, keywords, max_results, progress_callback)
        else:
            print("⚡ Using smart filtering")
            return fetch_preprints_smart_filter(server, start_date, end_date, keywords, max_results, progress_callback)

    except Exception as e:
        print(f"⚠️ Falling back to sequential fetch: {e}")
        return fetch_preprints(server, start_date, end_date, keywords, max_results)


def merge_and_highlight_articles(articles, fields_to_highlight=["Title", "Abstract"], raw_keywords=""):
    """Highlights keywords in articles' titles and abstracts."""
    def highlight(text, terms):
        for term in sorted(terms, key=len, reverse=True):
            pattern = re.compile(re.escape(term), re.IGNORECASE)
            text = pattern.sub(lambda m: f"**{m.group(0)}**", text)
        return text

    if not raw_keywords:
        return articles

    terms = re.findall(r'"[^"]+"|\b\w+\*?\??\b', raw_keywords)
    terms = [term.strip('"') for term in terms if term.upper() not in ["AND", "OR", "NOT"]]

    for article in articles:
        for field in fields_to_highlight:
            if field in article and isinstance(article[field], str):
                article[field] = highlight(article[field], terms)

    return articles


def standardize_date_format(articles):
    """Standardize date formats across all articles to YYYY-MM-DD."""
    for article in articles:
        for date_field in ('Publication Date', 'Date'):
            if date_field not in article or not article[date_field]:
                continue

            original_date = str(article[date_field])
            try:
                if '-' in original_date and any(m in original_date for m in
                        ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
                    article[date_field] = datetime.strptime(original_date, '%Y-%b-%d').strftime('%Y-%m-%d')

                elif '/' in original_date:
                    parts = original_date.split('/')
                    if len(parts) == 3:
                        year, month, day = parts
                        standardized = f"{year}-{month.zfill(2)}-{day.zfill(2)}"
                        datetime.strptime(standardized, '%Y-%m-%d')
                        article[date_field] = standardized

                elif '-' in original_date and len(original_date) == 10:
                    datetime.strptime(original_date, '%Y-%m-%d')

            except ValueError as e:
                print(f"Warning: Could not parse date '{original_date}': {e}")

    return articles


def standardize_doi_format(articles):
    """Strip URL prefixes from DOIs, keeping only the bare identifier."""
    prefixes = [
        'https://doi.org/', 'http://doi.org/',
        'https://dx.doi.org/', 'http://dx.doi.org/',
        'doi:', 'DOI:'
    ]

    for article in articles:
        if 'DOI' in article and article['DOI']:
            doi = str(article['DOI']).strip()
            for prefix in prefixes:
                if doi.startswith(prefix):
                    doi = doi[len(prefix):]
                    break
            article['DOI'] = doi

    return articles


def export_fetched_articles_as_csv(articles, journal, start_date, end_date, timestamp=None):
    df = pd.DataFrame(articles)
    filename = f"JournalTracker_{journal}_{start_date}_to_{end_date}"
    if timestamp:
        filename += f"_{timestamp}"
    filename += ".csv"
    df.to_csv(filename, index=False)
    print(f"📁 Saved {len(df)} articles to {filename}")


class OptimizedPubMedFetcher:
    def __init__(self, rate_limiter):
        self.rate_limiter = rate_limiter
        max_workers = 5 if rate_limiter.has_api_key else 2
        self.fetch_semaphore = threading.Semaphore(max_workers)
        self.max_workers = max_workers

    def fetch_articles_batch(self, pmid_batch, progress_callback=None):
        """Fetch multiple articles in a single batched Entrez request."""
        try:
            with self.fetch_semaphore:
                pmid_list = pmid_batch.split(',') if isinstance(pmid_batch, str) else pmid_batch

                self.rate_limiter._wait_for_rate_limit()

                delay = random.uniform(0.05, 0.1) if self.rate_limiter.has_api_key else random.uniform(0.5, 1.0)
                time.sleep(delay)

                handle = Entrez.efetch(db="pubmed", id=",".join(pmid_list), rettype="xml")
                paper_info_list = Entrez.read(handle)
                handle.close()

                papers = []

                for paper_info in paper_info_list.get("PubmedArticle", []):
                    try:
                        article_data = paper_info["MedlineCitation"]["Article"]
                        title = article_data["ArticleTitle"]

                        abstract_sections = article_data.get("Abstract", {}).get("AbstractText", [])
                        if abstract_sections:
                            if isinstance(abstract_sections[0], dict):
                                abstract = " ".join([s.get("text", "") for s in abstract_sections])
                            else:
                                abstract = str(abstract_sections[0])
                        else:
                            abstract = "No abstract available"

                        authors = []
                        for author in article_data.get("AuthorList", []):
                            last_name = author.get("LastName", "")
                            first_name = author.get("ForeName", "")
                            if last_name:
                                authors.append(f"{last_name}, {first_name}")

                        journal_info = article_data["Journal"]
                        journal_title = journal_info.get("Title", "")
                        journal_abbrev = journal_info.get("ISOAbbreviation", "")

                        journal_issue = journal_info.get("JournalIssue", {})
                        volume = journal_issue.get("Volume", "")
                        issue = journal_issue.get("Issue", "")
                        pages = article_data.get("Pagination", {}).get("MedlinePgn", "")

                        pub_date_str = "Unknown Date"

                        article_date_list = article_data.get("ArticleDate", [])
                        if article_date_list:
                            try:
                                d = article_date_list[0]
                                y, m, day = d.get("Year", ""), d.get("Month", "").zfill(2), d.get("Day", "").zfill(2)
                                if y and m and day:
                                    pub_date_str = f"{y}-{m}-{day}"
                            except Exception:
                                pass

                        if pub_date_str == "Unknown Date":
                            try:
                                pub_date = journal_info.get("JournalIssue", {}).get("PubDate", {})
                                y = pub_date.get("Year", "")
                                m = pub_date.get("Month", "")
                                day = pub_date.get("Day", "")

                                if m and not m.isdigit():
                                    month_map = {
                                        'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
                                        'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
                                        'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12',
                                        'January': '01', 'February': '02', 'March': '03', 'April': '04',
                                        'June': '06', 'July': '07', 'August': '08',
                                        'September': '09', 'October': '10', 'November': '11', 'December': '12'
                                    }
                                    m = month_map.get(m, m)

                                if y:
                                    if m and day:
                                        pub_date_str = f"{y}-{m.zfill(2)}-{day.zfill(2)}"
                                    elif m:
                                        pub_date_str = f"{y}-{m.zfill(2)}-01"
                                    else:
                                        pub_date_str = f"{y}-01-01"
                            except Exception:
                                pass

                        if pub_date_str == "Unknown Date":
                            try:
                                for date_entry in paper_info.get("PubmedData", {}).get("History", {}).get("PubMedPubDate", []):
                                    if date_entry.get("PubStatus") in ["pubmed", "entrez"]:
                                        y = date_entry.get("Year", "")
                                        m = date_entry.get("Month", "").zfill(2)
                                        day = date_entry.get("Day", "").zfill(2)
                                        if y and m and day:
                                            pub_date_str = f"{y}-{m}-{day}"
                                            break
                            except Exception:
                                pass

                        doi = None
                        for id_item in paper_info.get("PubmedData", {}).get("ArticleIdList", []):
                            if hasattr(id_item, 'attributes') and id_item.attributes.get("IdType") == "doi":
                                doi = str(id_item)
                                break

                        pmid_value = paper_info["MedlineCitation"]["PMID"]

                        year = month = day = ""
                        if pub_date_str != "Unknown Date":
                            try:
                                parts = pub_date_str.split('-')
                                year = parts[0] if len(parts) > 0 else ""
                                month = parts[1] if len(parts) > 1 else ""
                                day = parts[2] if len(parts) > 2 else ""
                            except Exception:
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
                            progress_callback(1)

                    except Exception as e:
                        print(f"⚠️ Error parsing article: {e}")
                        continue

                return papers

        except Exception as e:
            print(f"⚠️ Error fetching batch: {e}")
            return []

    def fetch_pubmed_articles_optimized(self, journal, start_date, end_date, keywords=None, progress_callback=None):
        """Optimized article fetching with batch processing and threading."""
        try:
            journal_dict = load_pubmed_journal_abbreviations()
            formatted_journal = format_journal_abbreviation(journal, journal_dict, self.rate_limiter)

            query = build_pubmed_query(formatted_journal, start_date, end_date, keywords)
            pmid_list, count = fetch_article_ids_from_pubmed(query, self.rate_limiter)

            if not pmid_list:
                return []

            batch_size = 200 if self.rate_limiter.has_api_key else 50
            pmid_batches = [pmid_list[i:i + batch_size] for i in range(0, len(pmid_list), batch_size)]

            all_papers = []

            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_batch = {
                    executor.submit(self.fetch_articles_batch, batch, progress_callback): batch
                    for batch in pmid_batches
                }

                for future in concurrent.futures.as_completed(future_to_batch):
                    try:
                        all_papers.extend(future.result())
                    except Exception as e:
                        print(f"⚠️ Batch processing error: {e}")

            return all_papers

        except Exception as e:
            print(f"⚠️ Error in optimized fetch: {e}")
            return []

    def fetch_articles_with_progress(self, journal_list, start_date, end_date, keywords=None, progress_callback=None):
        """Fetch articles from multiple journals with progress tracking."""
        all_articles = []

        for i, journal in enumerate(journal_list):
            try:
                print(f"🔍 Processing journal {i+1}/{len(journal_list)}: {journal}")

                def journal_progress_callback(count):
                    if progress_callback:
                        progress_callback(count, journal)

                articles = self.fetch_pubmed_articles_optimized(
                    journal, start_date, end_date, keywords, journal_progress_callback
                )

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
        """Return current rate limit and concurrency configuration."""
        rate_info = self.rate_limiter.get_rate_limit_info()
        return {
            "api_key_enabled": rate_info["has_api_key"],
            "max_concurrent_requests": self.max_workers,
            "batch_size": 200 if rate_info["has_api_key"] else 50,
            "total_api_requests": rate_info["total_requests"],
            "rate_limit_per_second": rate_info["max_requests_per_second"]
        }


def execute_subscription_search(journals, keywords, include_preprints, frequency):
    """Execute a search based on subscription parameters."""
    today = datetime.today().date()

    if frequency == "weekly":
        days_back = 7
    elif frequency == "monthly":
        days_back = 30
    elif frequency.startswith("every"):
        try:
            days_back = int(frequency.split()[1])
        except Exception:
            days_back = 7
    else:
        days_back = 7

    start_date = str(today - timedelta(days=days_back))
    end_date = str(today)

    print(f"🔍 Searching from {start_date} to {end_date}")

    all_articles = []

    try:
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

        if include_preprints:
            for server in ["biorxiv", "medrxiv"]:
                try:
                    preprints = fetch_preprints(
                        server=server,
                        start_date=start_date,
                        end_date=end_date,
                        keywords=keywords
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
            all_articles = standardize_date_format(all_articles)
            all_articles = standardize_doi_format(all_articles)
            merged = merge_and_highlight_articles(all_articles, [], keywords)
            return pd.DataFrame(merged)
        else:
            return pd.DataFrame()

    except Exception as e:
        print(f"❌ Error in subscription search: {e}")
        return pd.DataFrame()


def generate_bibtex_from_dataframe(df):
    """Generate BibTeX entries from a search results DataFrame."""
    bibtex_entries = []

    for index, row in df.iterrows():
        title = str(row.get('Title', 'Unknown Title')).replace('**', '')

        authors = row.get('Authors', 'Unknown Author')
        if pd.isna(authors) or authors == 'Unknown Author':
            authors = 'Unknown Author'

        journal = str(row.get('Journal', 'Unknown Journal'))

        year = str(row.get('Year', ''))
        if not year or year == 'nan':
            pub_date = str(row.get('Publication Date', ''))
            if pub_date and pub_date != 'nan':
                try:
                    year = pub_date.split('-')[0] if '-' in pub_date else pub_date[:4]
                except Exception:
                    year = 'Unknown Year'

        volume = str(row.get('Volume', ''))
        issue = str(row.get('Issue', ''))
        pages = str(row.get('Pages', ''))
        month = str(row.get('Month', ''))
        doi = str(row.get('DOI', ''))
        pmid = str(row.get('PMID', ''))
        abstract = str(row.get('Abstract', ''))

        first_author = authors.split(' and ')[0].split(',')[0].replace(' ', '') if authors != 'Unknown Author' else 'UnknownAuthor'
        first_word = title.split()[0] if title else 'Unknown'
        entry_key = f"{first_author}_{year}_{first_word}".replace(' ', '').replace('.', '')

        bibtex_entry = f"""@article{{{entry_key},
    title = {{{title}}},
    author = {{{authors}}},
    journal = {{{journal}}},
    year = {{{year}}}"""

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
        if abstract and abstract != 'nan' and len(abstract) > 10:
            if len(abstract) > 500:
                abstract = abstract[:500] + "..."
            bibtex_entry += f",\n    abstract = {{{abstract}}}"

        bibtex_entry += "\n}\n"
        bibtex_entries.append(bibtex_entry)

    return "\n".join(bibtex_entries)
