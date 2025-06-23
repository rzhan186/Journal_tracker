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
from datetime import datetime, timedelta

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

journal_dict = load_pubmed_journal_abbreviations()
# Flip dictionary: abbreviation -> full name
journal_fullname_to_abbrev = {v: k for k, v in journal_dict.items()}

with st.form("search_form"):
    email = st.text_input("📧 Your email (optional, for update subscription):")
    journal_inputs = st.multiselect(
        "Select journals (start typing to search):",
        options=list(journal_fullname_to_abbrev.keys()),
        help="You can add multiple journals."
    )

    date_option = st.selectbox(
        "Select date range:",
        ["Custom", "Past Week", "Past Month", "Past Year"]
    )

    start_date = ""
    end_date = ""
    if date_option == "Custom":
        start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
        end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
    else:
        today = datetime.today().date()
        if date_option == "Past Week":
            start_date = str(today - timedelta(days=7))
        elif date_option == "Past Month":
            start_date = str(today - timedelta(days=30))
        elif date_option == "Past Year":
            start_date = str(today - timedelta(days=365))
        end_date = str(today)

    raw_keywords = st.text_area("Keyword logic (optional)", height=100)
    subscribe = st.checkbox("📬 Subscribe to automatic updates")

    frequency = None
    custom_days = None

    if subscribe:
        freq_choice = st.selectbox("Update frequency:", ["daily", "weekly", "monthly", "custom"])
        if freq_choice == "custom":
            custom_days = st.number_input("🔧 Enter custom interval in days:", min_value=1, step=1)
            frequency = f"every {custom_days} days"
        else:
            frequency = freq_choice

    submitted = st.form_submit_button("🔍 Search")

if submitted:
    try:
        formatted_journals = []
        for j in journal_inputs:
            try:
                abbrev = journal_fullname_to_abbrev[j.strip()]
                formatted_journals.append(abbrev)
            except KeyError:
                st.error(f"❌ Error: '{j}' not found in PubMed journal list.")
                st.stop()

        keywords = None
        if raw_keywords.strip():
            if raw_keywords.count('(') != raw_keywords.count(')'):
                st.warning("⚠️ Unbalanced parentheses in your keyword logic.")
                st.stop()
            keywords = format_boolean_keywords_for_pubmed(raw_keywords.strip())

        st.info("Querying PubMed. Please wait...")
        all_articles = []

        for journal in formatted_journals:
            articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, keywords)
            for article in articles:
                article["Journal"] = journal  # tag article with journal name
            all_articles.extend(articles)

        if all_articles:
            df = pd.DataFrame(all_articles)
            st.success(f"✅ Found {len(df)} articles!")

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download CSV", csv, file_name="PubMed_Results.csv", mime="text/csv")

            if subscribe and email:
                res = store_user_subscription(
                    email=email,
                    journals=formatted_journals,
                    keywords=keywords,
                    start_date=start_date,
                    end_date=end_date,
                    frequency=frequency
                )
                st.success(f"📬 You will receive {frequency} updates at {email}.")
                st.write("🛠️ Supabase insert result:", res)

        else:
            st.warning("❌ No articles found. Try different keywords or date ranges.")

    except Exception as e:
        st.error(f"❌ An error occurred: {e}")
