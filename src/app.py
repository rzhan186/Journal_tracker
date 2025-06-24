
# app.py - Streamlit Web Interface
import streamlit as st
from tracking_main import (
    fetch_pubmed_articles_by_date,
    load_pubmed_journal_abbreviations,
    format_boolean_keywords_for_pubmed,
    build_pubmed_query
)

import pandas as pd
import os
from store_subscription import store_user_subscription
from datetime import datetime, timedelta
from dotenv import load_dotenv
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
st.title("📚 PubMed Journal Tracker")

st.markdown("""
Welcome to the **PubMed Journal Tracker**!  
Use this tool to search for articles by journal, date range, and keywords.
Subscribe to periodic updates from your journals and topics of interest.
""")

journal_dict = load_pubmed_journal_abbreviations()
email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance. Optional.")

# Journal selection
full_to_abbrev = {v: k for k, v in journal_dict.items()}
journal_options = list(full_to_abbrev.keys())

selected_journals = st.multiselect("📘 Select journal(s):", options=journal_options,
    help="Start typing to search journal names. Case-sensitive based on PubMed's names.")

# Date range
date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
today = datetime.today().date()
start_date = end_date = None

if date_option == "Custom":
    col1, col2 = st.columns(2)
    with col1:
        start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
    with col2:
        end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
else:
    if date_option == "Past Week":
        start_date = str(today - timedelta(days=7))
    elif date_option == "Past Month":
        start_date = str(today - timedelta(days=30))
    elif date_option == "Past Year":
        start_date = str(today - timedelta(days=365))
    end_date = str(today)

# Keyword input
raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
    help="""🔎 **Search Tips**  
- Use **AND**, **OR**, **NOT** for logic  
- Wrap phrases in **parentheses**: *(cadmium exposure)*  
- Wildcards: `metagenom*`, `wom?n` for flexible matching  
""")

# Track button visibility
if "show_search" not in st.session_state:
    st.session_state.show_search = True

if st.session_state.show_search:
    if st.button("🔍 Search"):
        try:
            if not selected_journals:
                st.error("❌ Please select at least one journal.")
                st.stop()

            formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)]
            keywords = None
            if raw_keywords.strip():
                if raw_keywords.count('(') != raw_keywords.count(')'):
                    st.warning("⚠️ Unbalanced parentheses in keyword logic.")
                    st.stop()
                keywords = format_boolean_keywords_for_pubmed(raw_keywords.strip())

            st.caption("🔍 Formatted PubMed keyword logic:")
            st.code(keywords if keywords else "(None)", language="none")

            if selected_journals:
                query_preview = build_pubmed_query(
                    journal=selected_journals[0],
                    start_date=start_date,
                    end_date=end_date,
                    keywords=keywords
                )
                st.caption("📄 Final PubMed query (1st journal shown):")
                st.code(query_preview, language="none")

            with st.status("🔍 Searching PubMed...", expanded=True) as status:
                all_articles = []
                for i, full_name in enumerate(selected_journals):
                    st.write(f"🔎 Searching: **{full_name}** ({i+1}/{len(selected_journals)})")
                    articles = fetch_pubmed_articles_by_date(full_name, start_date, end_date, keywords)
                    for article in articles:
                        article["Journal"] = full_name
                    all_articles.extend(articles)

                if all_articles:
                    status.update(label=f"✅ Found {len(all_articles)} article(s).", state="complete")
                else:
                    status.update(label="⚠️ No articles found. Try different keywords or dates.", state="error")

            if all_articles:
                df = pd.DataFrame(all_articles)
                st.success(f"✅ Found {len(df)} articles!")
                st.download_button("📥 Download CSV", df.to_csv(index=False).encode("utf-8"),
                    file_name="PubMed_Results.ƒSecsv", mime="text/csv")

        except Exception as e:
            st.error(f"❌ Error occurred: {e}")

# --- Subscription Section Below Search ---
subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

if subscribe:
    st.session_state.show_search = False

    col1, col2 = st.columns(2)
    with col1:
        freq_choice = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"])
    with col2:
        subscriber_email = st.text_input("📧 Email to receive updates")

    if freq_choice == "custom":
        custom_days = st.number_input("🔧 Custom interval (days):", min_value=1, step=1)
        frequency = f"every {custom_days} days"
    else:
        frequency = freq_choice

    st.markdown("### ✅ Confirm your subscription")
    st.info(f"""
**Email**: {subscriber_email or "Not provided"}  
**Journals**: {', '.join(selected_journals) if selected_journals else "None selected"}  
**Keywords**: {raw_keywords if raw_keywords else "None"}  
**Frequency**: {frequency}
""")

    if st.button("📩 Confirm and Subscribe"):
        if not subscriber_email:
            st.error("❌ Please provide an email address.")
        elif not selected_journals:
            st.error("❌ Please select at least one journal.")
        else:
            formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)]
            result = store_user_subscription(
                email=subscriber_email,
                journals=formatted_journals,
                keywords=raw_keywords,
                start_date=start_date,
                end_date=end_date,
                frequency=frequency
            )
            st.success(f"📬 Subscribed! You'll receive {frequency} updates at {subscriber_email}.")
            st.write("🛠️ Supabase insert result:", result)
