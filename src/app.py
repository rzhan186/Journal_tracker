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
- 🔍 Keyword logic is optional and supports advanced combinations like `(cadmium OR "cadmium exposure") AND rice`
""")

journal_dict = load_pubmed_journal_abbreviations()

with st.form("search_form"):
    st.markdown("**📧 Email Address** *(Optional – required by NCBI Entrez API, used only for contact)*")
    email = st.text_input("Enter your email (Optional):",
        help="Required by NCBI Entrez API. This is optional and only used for API compliance.")

    journal_options = list(journal_dict.keys())
    selected_journals = st.multiselect(
        "📘 Select journal(s):",
        options=journal_options,
        help="Start typing to search and select one or more journal names. Case-sensitive based on PubMed official names.")

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

    raw_keywords = st.text_area("Enter your search keyword (Optional) :", height=100,
                                help="Use AND, OR, NOT. Wrap phrases in quotes. E.g., (cadmium OR \"cadmium exposure\") AND rice")

    subscribe = st.checkbox("📬 Subscribe to automatic updates")

    frequency = None
    custom_days = None
    subscriber_email = None

    if subscribe:
        st.markdown("**🔁 Update Frequency**")
        freq_choice = st.selectbox("How often do you want to receive updates?", ["weekly", "monthly", "custom"])
        if freq_choice == "custom":
            custom_days = st.number_input("🔧 Enter custom interval in days:", min_value=1, step=1)
            frequency = f"every {custom_days} days"
        else:
            frequency = freq_choice

        subscriber_email = st.text_input("📧 Email to receive updates:",
                                         help="Where your automated journal updates will be sent")

    submitted = st.form_submit_button("🔍 Search")

if submitted:
    try:
        formatted_journals = []
        for name in selected_journals:
            if name in journal_dict:
                formatted_journals.append(journal_dict[name])
            else:
                st.error(f"❌ Error: '{name}' not found in PubMed journal list.")
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

            if subscribe and subscriber_email:
                res = store_user_subscription(
                    email=subscriber_email,
                    journals=formatted_journals,
                    keywords=keywords,
                    start_date=start_date,
                    end_date=end_date,
                    frequency=frequency
                )
                st.success(f"📬 You will receive {frequency} updates at {subscriber_email}.")
                st.write("🛠️ Supabase insert result:", res)

        else:
            st.warning("❌ No articles found. Try different keywords or date ranges.")

    except Exception as e:
        st.error(f"❌ An error occurred: {e}")
