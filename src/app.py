# app.py - Streamlit Web Interface
import streamlit as st
from tracking_main import (
    fetch_pubmed_articles_by_date,
    load_pubmed_journal_abbreviations,
    format_journal_abbreviation,
    format_boolean_keywords_for_pubmed,
)
import pandas as pd

st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
st.title("📚 PubMed Journal Tracker")

st.markdown("""
Enter your search criteria below. 
- Use Boolean operators: AND, OR, NOT  
- Use quotes for phrases: "global warming"  
- ❗ Wildcards (*, ?) are not supported in PubMed API
""")

with st.form("search_form"):
    journal = st.text_input("Journal name or abbreviation (e.g., Environmental Science & Technology or  Environ Sci Technol):")
    start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
    end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
    raw_keywords = st.text_area("Keyword logic (optional)", height=100)
    submitted = st.form_submit_button("🔍 Search")

if submitted:
    try:
        journal_dict = load_pubmed_journal_abbreviations()
        formatted_journal = format_journal_abbreviation(journal.strip(), journal_dict)

        keywords = None
        if raw_keywords.strip():
            if raw_keywords.count('(') != raw_keywords.count(')'):
                st.warning("⚠️ Unbalanced parentheses in your keyword logic.")
                st.stop()
            keywords = format_boolean_keywords_for_pubmed(raw_keywords.strip())

        st.info("Querying PubMed. Please wait...")
        articles = fetch_pubmed_articles_by_date(formatted_journal, start_date, end_date, keywords)

        if articles:
            df = pd.DataFrame(articles)
            st.success(f"✅ Found {len(df)} articles!")

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download CSV", csv, file_name="PubMed_Results.csv", mime="text/csv")
        else:
            st.warning("❌ No articles found. Try different keywords or date ranges.")

    except Exception as e:
        st.error(f"❌ An error occurred: {e}")
