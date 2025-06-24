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
Welcome to the **PubMed Journal Tracker**!  
Use this tool to search for articles by journal, date range, and keywords.
Subscribe to periodic update from your journals and topic of interests.
""")

journal_dict = load_pubmed_journal_abbreviations()

email = st.text_input("📧 Enter your email (Optional):",
    help="Required by NCBI Entrez API. This is optional and only used for API compliance.")

# Journal selection
full_to_abbrev = {v: k for k, v in journal_dict.items()}
journal_options = list(full_to_abbrev.keys())

selected_journals = st.multiselect(
    "📘 Select journal(s):",
    options=journal_options,
    help="Start typing to search and select one or more journal names. Case-sensitive based on PubMed official names."
)

# Date range selection
date_option = st.selectbox(
    "📅 Select date range:",
    ["Past Week", "Past Month", "Past Year", "Custom"],
    index=0
)

today = datetime.today().date()
start_date = None
end_date = None

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
raw_keywords = st.text_area("❓ Enter your search keyword (Optional) :", height=100,
    help="""🔎 **Search Tips**  
        - Use **AND**, **OR**, **NOT** for logical combinations.  
        - Wrap phrases in **parentheses**: *(cadmium exposure)*  
        - Use **wildcards** for flexible matches:  
        • Asterisk `*` matches multiple characters → `metagenom*` matches *metagenome*, *metagenomics*  
        • Question mark `?` matches a single character → `wom?n` matches *woman*, *women*  
        """)

# Subscription options
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

# Manual search button
submitted = st.button("🔍 Search")

if submitted:
    try:
        if not selected_journals:
            st.error("❌ Please select at least one journal.")
            st.stop()

        formatted_journals = []
        for full_name in selected_journals:
            abbrev = full_to_abbrev.get(full_name)
            if abbrev:
                formatted_journals.append(abbrev)
            else:
                st.error(f"❌ '{full_name}' not found in the journal list.")
                st.stop()

        keywords = None

        if raw_keywords.strip():
            if raw_keywords.count('(') != raw_keywords.count(')'):
                st.warning("⚠️ Unbalanced parentheses in your keyword logic.")
                st.stop()
            keywords = format_boolean_keywords_for_pubmed(raw_keywords.strip())

        # ✅ Show formatted keyword query (whether or not keywords were entered)
        st.caption("🔍 Formatted PubMed keyword logic:")
        st.code(keywords if keywords else "(None)", language="none")

        # ✅ Optional full query preview
        if selected_journals:
            query_preview = build_pubmed_query(
                journal=selected_journals[0],
                start_date=start_date,
                end_date=end_date,
                keywords=keywords
            )
            st.caption("📄 Final PubMed query (1st journal shown):")
            st.code(query_preview, language="none")


        with st.status("🔍 Starting PubMed search...", expanded=True) as status:
            all_articles = []

            for i, full_name in enumerate(selected_journals):
                st.write(f"📖 Searching: **{full_name}** ({i+1}/{len(selected_journals)})")

                articles = fetch_pubmed_articles_by_date(full_name, start_date, end_date, keywords)
                for article in articles:
                    article["Journal"] = full_name
                all_articles.extend(articles)

            if all_articles:
                status.update(label=f"✅ Found {len(all_articles)} article(s).", state="complete")
            else:
                status.update(label="⚠️ No articles found. Try different keywords or date ranges.", state="error")

        if all_articles:
            df = pd.DataFrame(all_articles)
            st.success(f"✅ Found {len(df)} articles!")

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download CSV", csv, file_name="PubMed_Results.csv", mime="text/csv")

            # --- Subscription section shown only if articles are found ---
            if all_articles:
                st.divider()
                st.subheader("📬 Subscribe to Automatic Updates")

                with st.form("subscribe_form"):
                    st.markdown("Stay informed with regular updates from your selected journals.")

                    sub_email = st.text_input("📧 Email to receive updates:", value=subscriber_email or "",
                                            help="Where your automated journal updates will be sent")

                    freq = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"], index=0)
                    if freq == "custom":
                        interval_days = st.number_input("🔧 Custom interval (in days):", min_value=1, step=1)
                        freq_text = f"every {interval_days} days"
                    else:
                        freq_text = freq

                    confirm = st.form_submit_button("📨 Subscribe Now")

                    if confirm:
                        if not sub_email:
                            st.warning("❗ Please enter an email address to subscribe.")
                        else:
                            # Insert into Supabase (or your backend)
                            result = store_user_subscription(
                                email=sub_email,
                                journals=formatted_journals,
                                keywords=keywords,
                                start_date=start_date,
                                end_date=end_date,
                                frequency=freq_text
                            )

                            st.success(f"✅ You’re now subscribed to updates {freq_text} at **{sub_email}**.")
                            st.caption("You'll receive updates based on your current selection:")
                            st.write(f"**Journals:** {', '.join(selected_journals)}")
                            st.write(f"**Date range:** {start_date} to {end_date}")
                            st.write(f"**Keywords:** {raw_keywords if raw_keywords else '(none)'}")
                            st.write("🛠️ Supabase response:", result)

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