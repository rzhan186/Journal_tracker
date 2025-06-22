# app.py - Streamlit Web Interface
import streamlit as st
from tracking_main import (
    fetch_pubmed_articles_by_date,
    load_pubmed_journal_abbreviations,
    format_journal_abbreviation,
    format_boolean_keywords_for_pubmed,
)
import pandas as pd
import os
from store_subscription import store_user_subscription  # New: for Supabase integration

# Optional: Load .env locally
from dotenv import load_dotenv
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
st.title("📚 PubMed Journal Tracker")

st.markdown("""
Enter your search criteria below. 
- Use Boolean operators: AND, OR, NOT  
- Use quotes for phrases: "global warming"  
- ❗ Wildcards (*, ?) are not supported in PubMed API
""")

with st.form("search_form"):
    email = st.text_input("📧 Your email (optional, for update subscription):")
    journal = st.text_input("Journal name or abbreviation (e.g., Environmental Science & Technology or Environ Sci Technol):")
    start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
    end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
    raw_keywords = st.text_area("Keyword logic (optional)", height=100)
    subscribe = st.checkbox("📬 Subscribe to automatic updates")
    frequency = st.selectbox("Update frequency (if subscribed):", ["daily", "weekly", "monthly", "custom"], disabled=not subscribe)
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

            if subscribe and email:
                store_user_subscription(
                    email=email,
                    journals=[formatted_journal],
                    keywords=keywords,
                    start_date=start_date,
                    end_date=end_date,
                    frequency=frequency
                )
                st.success(f"📬 You will receive {frequency} updates at {email}.")

        else:
            st.warning("❌ No articles found. Try different keywords or date ranges.")

    except Exception as e:
        st.error(f"❌ An error occurred: {e}")
