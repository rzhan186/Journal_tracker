# # app.py - Streamlit Web Interface

# import streamlit as st
# from tracking_main import (
#     fetch_pubmed_articles_by_date,
#     load_pubmed_journal_abbreviations,
#     format_boolean_keywords_for_pubmed,
#     build_pubmed_query,
#     generate_placeholder_csv,
# )
# import pandas as pd
# import os
# from store_subscription import store_user_subscription

# from datetime import datetime, timedelta
# from dotenv import load_dotenv
# from itsdangerous import URLSafeSerializer, BadSignature

# from email_dispatcher import send_email

# load_dotenv()

# EMAIL_ADDRESS = st.secrets["EMAIL_ADDRESS"]
# EMAIL_PASSWORD = st.secrets["EMAIL_PASSWORD"]

# SUPABASE_URL = os.getenv("SUPABASE_URL")
# SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
# serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

# BASE_URL = "https://journaltracker.streamlit.app"

# # Streamlit app configuration
# st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
# st.title("📚 PubMed Journal Tracker")

# # Check for unsubscribe token in the URL
# if 'token' in st.query_params:
#     token = st.query_params['token']  # Get the token from the query parameters
    
#     # Call the unsubscribe handling function
#     from unsubscribe import handle_unsubscribe  # Make sure this imports the unsubscribe logic
#     handle_unsubscribe(token)  # Call the function in unsubscribe.py
    
# else:
#     # If no token is found, display the main application interface

#     st.markdown("""
#     Use this tool to search PubMed by journal, date range, and keywords.  
#     You can also subscribe to automatic updates.
#     """)

#     journal_dict = load_pubmed_journal_abbreviations()
#     full_to_abbrev = {v: k for k, v in journal_dict.items()}
#     journal_options = list(full_to_abbrev.keys())

#     email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance.")
#     selected_journals = st.multiselect("📘 Select journal(s):", options=journal_options)

#     date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
#     today = datetime.today().date()

#     if date_option == "Custom":
#         col1, col2 = st.columns(2)
#         with col1:
#             start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
#         with col2:
#             end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
#     else:
#         days = {"Past Week": 7, "Past Month": 30, "Past Year": 365}[date_option]
#         start_date = str(today - timedelta(days=days))
#         end_date = str(today)

#     raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
#         help="""🔎 **Search Tips**  
#         - Use **AND**, **OR**, **NOT**  
#         - Wrap phrases in **parentheses**: *(cadmium exposure)*  
#         - Wildcards: `metagenom*`, `wom?n`
#         """)

#     # Manual Search Button
#     if st.button("🔍 Search"):
#         try:
#             if not selected_journals:
#                 st.error("❌ Please select at least one journal.")
#                 st.stop()

#             formatted_journals = [full_to_abbrev.get(j) for j in selected_journals if j in full_to_abbrev]
#             if not formatted_journals:
#                 st.error("❌ Invalid journal names.")
#                 st.stop()

#             # Validate keywords
#             keywords = None
#             if raw_keywords.strip():
#                 if raw_keywords.count("(") != raw_keywords.count(")"):
#                     st.warning("⚠️ Unbalanced parentheses.")
#                     st.stop()
#                 keywords = format_boolean_keywords_for_pubmed(raw_keywords.strip())

#             st.caption("🔍 Formatted PubMed keyword logic:")
#             st.code(keywords if keywords else "(None)", language="none")

#             query_preview = build_pubmed_query(
#                 journal=selected_journals[0], start_date=start_date, end_date=end_date, keywords=keywords
#             )
#             st.caption("📄 Final PubMed query (1st journal shown):")
#             st.code(query_preview, language="none")

#             with st.status("🔍 Searching PubMed...", expanded=True) as status:
#                 all_articles = []
#                 for i, journal in enumerate(selected_journals):
#                     st.write(f"🔎 Searching: **{journal}** ({i+1}/{len(selected_journals)})")
#                     articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, keywords)
#                     for article in articles:
#                         article["Journal"] = journal
#                     all_articles.extend(articles)

#                 if all_articles:
#                     status.update(label=f"✅ Found {len(all_articles)} article(s).", state="complete")
#                 else:
#                     status.update(label="⚠️ No articles found.", state="error")

#             if all_articles:
#                 df = pd.DataFrame(all_articles)
#                 st.download_button("📥 Download CSV", df.to_csv(index=False).encode("utf-8"),
#                                    file_name="PubMed_Results.csv", mime="text/csv")
#         except Exception as e:
#             st.error(f"❌ Error: {e}")

#     # --- Subscribe toggle ---
#     subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

#     # Subscription section only renders when checked
#     if subscribe:
#         col1, col2 = st.columns(2)
#         with col1:
#             freq_choice = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"])
#         with col2:
#             subscriber_email = st.text_input("📧 Email to receive updates")

#         if freq_choice == "custom":
#             custom_days = st.number_input("🔧 Custom interval (days):", min_value=1, step=1)
#             frequency = f"every {custom_days} days"
#         else:
#             frequency = freq_choice

#         st.markdown("✅ Confirm your subscription")
#         st.info(f"""
#                 **Email**: {subscriber_email or "Not provided"}  
#                 **Journals**: {', '.join(selected_journals) if selected_journals else "None selected"}  
#                 **Keywords**: {raw_keywords if raw_keywords else "None"}  
#                 **Frequency**: {frequency}
#                 """)

#         if st.button("📩 Confirm and Subscribe"):
#             if not subscriber_email:
#                 st.error("❌ Please provide an email address.")
#             elif not selected_journals:
#                 st.error("❌ Please select at least one journal.")
#             else:
#                 formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)]
#                 csv_bytes = df.to_csv(index=False).encode("utf-8") if "df" in locals() else generate_placeholder_csv()

#                 result = store_user_subscription(
#                     email=subscriber_email,
#                     journals=formatted_journals,
#                     keywords=raw_keywords,
#                     start_date=start_date,
#                     end_date=end_date,
#                     frequency=frequency,
#                 )
#                 st.success(f"📬 Subscribed! You'll receive {frequency} updates at {subscriber_email}.")
#                 st.write("🛠️ Supabase insert result:", result)

#                 if result["status"] == "success":
#                     unsubscribe_token = result["unsubscribe_token"]
#                     unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"

#                     email_body = f"""Hi {subscriber_email},

#                 You have successfully subscribed to automatic PubMed updates.

#                 📘 Journals: {', '.join(formatted_journals)}
#                 🔑 Keywords: {raw_keywords or 'None'}
#                 🔁 Frequency: {frequency}
#                 📅 Date Range: {start_date} to {end_date}

#                 If you wish to unsubscribe, click the link below:
#                 🔓 {unsubscribe_link}

#                 – PubMed Tracker Team
#                     """
                    
#                     try:
#                         send_email(
#                             to_email=subscriber_email,
#                             subject="📬 Journal Tracker: Subscription Confirmed",
#                             body=email_body
#                         )
#                         st.success("✅ A confirmation email has been sent.")
#                     except Exception as e:
#                         st.warning(f"⚠️ Subscription saved, but email failed: {e}")


# app.py - Streamlit Web Interface

# import streamlit as st
# from tracking_main import (
#     fetch_pubmed_articles_by_date,
#     fetch_preprints,
#     load_pubmed_journal_abbreviations,
#     format_boolean_keywords_for_pubmed,
#     build_pubmed_query,
#     generate_placeholder_csv,
#     merge_and_highlight_articles,
#     compile_keyword_filter)

# import pandas as pd
# import os
# from store_subscription import store_user_subscription

# from datetime import datetime, timedelta
# from dotenv import load_dotenv
# from itsdangerous import URLSafeSerializer, BadSignature

# from email_dispatcher import send_email

# load_dotenv()

# EMAIL_ADDRESS = st.secrets["EMAIL_ADDRESS"]
# EMAIL_PASSWORD = st.secrets["EMAIL_PASSWORD"]

# SUPABASE_URL = os.getenv("SUPABASE_URL")
# SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
# serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

# BASE_URL = "https://journaltracker.streamlit.app"

# # Streamlit app configuration
# st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
# st.title("📚 PubMed Journal Tracker")

# # Check for unsubscribe token in the URL
# if 'token' in st.query_params:
#     token = st.query_params['token']  # Get the token from the query parameters
    
#     # Call the unsubscribe handling function
#     from unsubscribe import handle_unsubscribe  # Make sure this imports the unsubscribe logic
#     handle_unsubscribe(token)  # Call the function in unsubscribe.py
    
# else:
#     # If no token is found, display the main application interface

#     st.markdown("""
#     Use this tool to search PubMed by journal, date range, and keywords.  
#     You can also subscribe to automatic updates.
#     """)

#     journal_dict = load_pubmed_journal_abbreviations()
#     full_to_abbrev = {v: k for k, v in journal_dict.items()}
#     journal_options = list(full_to_abbrev.keys())

#     email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance.")
#     col1, col2 = st.columns([2, 1])
#     with col1:
#         selected_journals = st.multiselect("📘 Select journal(s):", options=journal_options)
#     with col2:
#         include_preprints = st.checkbox("📑 Include preprints", help="Currently supports bioRxiv and medRxiv.")

#     date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
#     today = datetime.today().date()

#     if date_option == "Custom":
#         col1, col2 = st.columns(2)
#         with col1:
#             start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
#         with col2:
#             end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
#     else:
#         days = {"Past Week": 7, "Past Month": 30, "Past Year": 365}[date_option]
#         start_date = str(today - timedelta(days=days))
#         end_date = str(today)

#     raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
#         help="""🔎 **Search Tips**  
#         - Use **AND**, **OR**, **NOT**  
#         - Wrap phrases in **parentheses**: *(cadmium exposure)*  
#         - Wildcards: `metagenom*`, `wom?n`
#         """)

#     # Manual Search Button
#     if st.button("🔍 Search"):
#         try:
#             if not selected_journals:
#                 st.error("❌ Please select at least one journal.")
#                 st.stop()

#             formatted_journals = [full_to_abbrev.get(j) for j in selected_journals if j in full_to_abbrev]
#             if not formatted_journals:
#                 st.error("❌ Invalid journal names.")
#                 st.stop()

#             if raw_keywords.strip():
#                 if raw_keywords.count("(") != raw_keywords.count(")"):
#                     st.warning("⚠️ Unbalanced parentheses.")
#                     st.stop()
#                 user_keywords = raw_keywords.strip()
#                 compiled_filter = compile_keyword_filter(raw_keywords)
#                 pubmed_keywords = format_boolean_keywords_for_pubmed(raw_keywords)

#                 keywords = pubmed_keywords
#             else:
#                 keywords = None

#             if keywords:
#                 query_preview = build_pubmed_query(
#                     journal=selected_journals[0],
#                     start_date=start_date,
#                     end_date=end_date,
#                     keywords=keywords
#                 )
#                 st.caption("📄 Final PubMed query (1st journal shown):")
#                 st.code(query_preview, language="none")

#             all_articles = []
#             for journal in selected_journals:
#                 articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, pubmed_keywords)

#                 for article in articles:
#                     article["Journal"] = journal
#                 all_articles.extend(articles)

#             if include_preprints:
#                 for server in ["biorxiv", "medrxiv"]:
#                     preprints = fetch_preprints(
#                         server=server,
#                         start_date=start_date,
#                         end_date=end_date,
#                         keywords=raw_keywords
#                     )
#                     for article in preprints:
#                         article["Journal"] = server
#                         article["Source"] = "Preprint"
#                     all_articles.extend(preprints)

#             if all_articles:
#                 merged = merge_and_highlight_articles(all_articles, [], raw_keywords)
#                 for article in all_articles:
#                     if "source" not in article:
#                         article["Source"] = "PubMed"  # Add this field to match CSV header

#                 df = pd.DataFrame(merge_and_highlight_articles(all_articles, [], raw_keywords))

#                 st.success(f"✅ Found {len(df)} articles.")
#                 st.download_button("📥 Download CSV", df.to_csv(index=False).encode("utf-8"), file_name="Combined_Results.csv")
#             else:
#                 st.warning("⚠️ No articles found.")

#         except Exception as e:
#             st.error(f"❌ Error: {e}")

#     # --- Subscribe toggle ---
#     subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

#     # Subscription section only renders when checked
#     if subscribe:
#         col1, col2 = st.columns(2)
#         with col1:
#             freq_choice = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"])
#         with col2:
#             subscriber_email = st.text_input("📧 Email to receive updates")

#         if freq_choice == "custom":
#             custom_days = st.number_input("🔧 Custom interval (days):", min_value=1, step=1)
#             frequency = f"every {custom_days} days"
#         else:
#             frequency = freq_choice

#         st.markdown("✅ Confirm your subscription")
#         st.info(f"""
#                 **Email**: {subscriber_email or "Not provided"}  
#                 **Journals**: {', '.join(selected_journals) if selected_journals else "None selected"}  
#                 **Keywords**: {raw_keywords if raw_keywords else "None"}  
#                 **Frequency**: {frequency}
#                 """)

#         if st.button("📩 Confirm and Subscribe"):
#             if not subscriber_email:
#                 st.error("❌ Please provide an email address.")
#             elif not selected_journals:
#                 st.error("❌ Please select at least one journal.")
#             else:
#                 formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)]
#                 csv_bytes = df.to_csv(index=False).encode("utf-8") if "df" in locals() and not df.empty else generate_placeholder_csv()

#                 result = store_user_subscription(
#                     email=subscriber_email,
#                     journals=formatted_journals,
#                     keywords=raw_keywords,
#                     start_date=start_date,
#                     end_date=end_date,
#                     frequency=frequency,
#                 )
#                 st.success(f"📬 Subscribed! You'll receive {frequency} updates at {subscriber_email}.")
#                 st.write("🛠️ Supabase insert result:", result)

#                 if result["status"] == "success":
#                     unsubscribe_token = result["unsubscribe_token"]
#                     unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"

#                     email_body = f"""Hi {subscriber_email},

#                 You have successfully subscribed to automatic PubMed updates.

#                 📘 Journals: {', '.join(formatted_journals)}
#                 🔑 Keywords: {raw_keywords or 'None'}
#                 🔁 Frequency: {frequency}
#                 📅 Date Range: {start_date} to {end_date}

#                 If you wish to unsubscribe, click the link below:
#                 🔓 {unsubscribe_link}

#                 – PubMed Tracker Team
#                     """
                    
#                     try:
#                         send_email(
#                             to_email=subscriber_email,
#                             subject="📬 Journal Tracker: Subscription Confirmed",
#                             body=email_body
#                         )
#                         st.success("✅ A confirmation email has been sent.")
#                     except Exception as e:
#                         st.warning(f"⚠️ Subscription saved, but email failed: {e}")


# import streamlit as st
# from tracking_main import (
#     fetch_pubmed_articles_by_date,
#     fetch_preprints,
#     load_pubmed_journal_abbreviations,
#     format_boolean_keywords_for_pubmed,
#     build_pubmed_query,
#     generate_placeholder_csv,
#     merge_and_highlight_articles,
#     compile_keyword_filter)

# import pandas as pd
# import os
# from store_subscription import store_user_subscription

# from datetime import datetime, timedelta
# from dotenv import load_dotenv
# from itsdangerous import URLSafeSerializer, BadSignature

# from email_dispatcher import send_email

# load_dotenv()

# EMAIL_ADDRESS = st.secrets["EMAIL_ADDRESS"]
# EMAIL_PASSWORD = st.secrets["EMAIL_PASSWORD"]

# SUPABASE_URL = os.getenv("SUPABASE_URL")
# SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
# serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

# BASE_URL = "https://journaltracker.streamlit.app"

# # Streamlit app configuration
# st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
# st.title("📚 PubMed Journal Tracker")

# # Check for unsubscribe token in the URL
# if 'token' in st.query_params:
#     token = st.query_params['token']  # Get the token from the query parameters
    
#     # Call the unsubscribe handling function
#     from unsubscribe import handle_unsubscribe  # Make sure this imports the unsubscribe logic
#     handle_unsubscribe(token)  # Call the function in unsubscribe.py
    
# else:
#     # If no token is found, display the main application interface

#     st.markdown("""
#     Use this tool to search PubMed by journal, date range, and keywords.  
#     You can also subscribe to automatic updates.
#     """)

#     journal_dict = load_pubmed_journal_abbreviations()
#     full_to_abbrev = {v: k for k, v in journal_dict.items()}
#     journal_options = list(full_to_abbrev.keys())

#     email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance.")
#     col1, col2 = st.columns([2, 1])
#     with col1:
#         selected_journals = st.multiselect("📘 Select journal(s):", options=journal_options)
#     with col2:
#         include_preprints = st.checkbox("📑 Include preprints", help="Currently supports bioRxiv and medRxiv.")

#     date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
#     today = datetime.today().date()

#     if date_option == "Custom":
#         col1, col2 = st.columns(2)
#         with col1:
#             start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
#         with col2:
#             end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
#     else:
#         days = {"Past Week": 7, "Past Month": 30, "Past Year": 365}[date_option]
#         start_date = str(today - timedelta(days=days))
#         end_date = str(today)

#     raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
#         help="""🔎 **Search Tips**  
#         - Use **AND**, **OR**, **NOT**  
#         - Wrap phrases in **parentheses**: *(cadmium exposure)*  
#         - Wildcards: `metagenom*`, `wom?n`
#         """)

#     # Manual Search Button
#     if st.button("🔍 Search"):
#         try:
#             if not selected_journals:
#                 st.error("❌ Please select at least one journal.")
#                 st.stop()

#             formatted_journals = [full_to_abbrev.get(j) for j in selected_journals if j in full_to_abbrev]
#             if not formatted_journals:
#                 st.error("❌ Invalid journal names.")
#                 st.stop()

#             # Validate keywords early if provided
#             if raw_keywords.strip():
#                 if raw_keywords.count("(") != raw_keywords.count(")"):
#                     st.warning("⚠️ Unbalanced parentheses.")
#                     st.stop()
#                 user_keywords = raw_keywords.strip()
#                 compiled_filter = compile_keyword_filter(raw_keywords)
#                 pubmed_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
#                 keywords = pubmed_keywords
#             else:
#                 keywords = None

#             # Create placeholders for status updates
#             status_placeholder = st.empty()
#             progress_bar = st.progress(0)
#             details_container = st.container()
            
#             # Show search details in collapsible sections
#             with details_container:
#                 with st.expander("🔍 Search Details", expanded=False):
#                     st.write(f"**Selected Journals:** {', '.join(selected_journals)}")
#                     st.write(f"**Date Range:** {start_date} to {end_date}")
#                     st.write(f"**Keywords:** {raw_keywords if raw_keywords else 'None'}")
#                     st.write(f"**Include Preprints:** {'Yes' if include_preprints else 'No'}")
                    
#                     if keywords:
#                         query_preview = build_pubmed_query(
#                             journal=selected_journals[0],
#                             start_date=start_date,
#                             end_date=end_date,
#                             keywords=keywords
#                         )
#                         st.caption("📄 PubMed query example (first journal):")
#                         st.code(query_preview, language="none")

#             # Initialize search
#             status_placeholder.info("🚀 Starting search...")
#             progress_bar.progress(10)
            
#             all_articles = []
#             total_journals = len(selected_journals)
#             preprint_servers = ["biorxiv", "medrxiv"] if include_preprints else []
#             total_sources = total_journals + len(preprint_servers)
            
#             current_step = 0
            
#             # Search PubMed journals
#             for i, journal in enumerate(selected_journals):
#                 current_step += 1
#                 status_placeholder.info(f"🔍 Searching {journal}... ({current_step}/{total_sources})")
#                 progress_bar.progress(int((current_step / total_sources) * 80))
                
#                 try:
#                     articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, pubmed_keywords)
#                     for article in articles:
#                         article["Journal"] = journal
#                     all_articles.extend(articles)
                    
#                     # Show intermediate results
#                     if articles:
#                         status_placeholder.success(f"✅ Found {len(articles)} articles in {journal}")
#                     else:
#                         status_placeholder.info(f"📭 No articles found in {journal}")
                        
#                 except Exception as e:
#                     st.error(f"❌ Error searching {journal}: {str(e)}")
#                     continue

#             # Search preprints if requested
#             if include_preprints:
#                 for server in preprint_servers:
#                     current_step += 1
#                     status_placeholder.info(f"🔍 Searching {server} preprints... ({current_step}/{total_sources})")
#                     progress_bar.progress(int((current_step / total_sources) * 80))
                    
#                     try:
#                         preprints = fetch_preprints(
#                             server=server,
#                             start_date=start_date,
#                             end_date=end_date,
#                             keywords=raw_keywords
#                         )
#                         for article in preprints:
#                             article["Journal"] = server
#                             article["Source"] = "Preprint"
#                         all_articles.extend(preprints)
                        
#                         # Show intermediate results
#                         if preprints:
#                             status_placeholder.success(f"✅ Found {len(preprints)} preprints in {server}")
#                         else:
#                             status_placeholder.info(f"📭 No preprints found in {server}")
                            
#                     except Exception as e:
#                         st.error(f"❌ Error searching {server}: {str(e)}")
#                         continue

#             # Finalize results
#             progress_bar.progress(90)
#             status_placeholder.info("📊 Processing results...")

#             if all_articles:
#                 # Add Source field for PubMed articles
#                 for article in all_articles:
#                     if "Source" not in article:
#                         article["Source"] = "PubMed"

#                 merged = merge_and_highlight_articles(all_articles, [], raw_keywords)
#                 df = pd.DataFrame(merged)
                
#                 progress_bar.progress(100)
#                 status_placeholder.success(f"🎉 Search completed! Found {len(df)} articles total.")
                
#                 # Show results summary
#                 st.success(f"✅ Found {len(df)} articles across all sources.")
                
#                 # Results breakdown
#                 with st.expander("📊 Results Breakdown", expanded=True):
#                     if "Source" in df.columns:
#                         source_counts = df["Source"].value_counts()
#                         for source, count in source_counts.items():
#                             st.write(f"**{source}:** {count} articles")
                    
#                     if "Journal" in df.columns:
#                         journal_counts = df["Journal"].value_counts()
#                         st.write("\n**By Journal/Source:**")
#                         for journal, count in journal_counts.items():
#                             st.write(f"- {journal}: {count} articles")

#                 # Download button
#                 csv_data = df.to_csv(index=False).encode("utf-8")
#                 st.download_button(
#                     label="📥 Download Results as CSV",
#                     data=csv_data,
#                     file_name=f"search_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
#                     mime="text/csv"
#                 )
                
#                 # Show sample results
#                 with st.expander("👀 Preview Results", expanded=False):
#                     st.dataframe(df.head(10))
                
#             else:
#                 progress_bar.progress(100)
#                 status_placeholder.warning("📭 No articles found matching your criteria.")
                
#                 # Show suggestions
#                 st.info("""
#                 **No results found. Try:**
#                 - Expanding your date range
#                 - Using broader keywords
#                 - Checking different journals
#                 - Including preprints
#                 """)

#         except Exception as e:
#             st.error(f"❌ Unexpected error: {e}")
#             # Clear progress indicators on error
#             if 'status_placeholder' in locals():
#                 status_placeholder.empty()
#             if 'progress_bar' in locals():
#                 progress_bar.empty()

#     # --- Subscribe toggle ---
#     subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

#     # Subscription section only renders when checked
#     if subscribe:
#         col1, col2 = st.columns(2)
#         with col1:
#             freq_choice = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"])
#         with col2:
#             subscriber_email = st.text_input("📧 Email to receive updates")

#         if freq_choice == "custom":
#             custom_days = st.number_input("🔧 Custom interval (days):", min_value=1, step=1)
#             frequency = f"every {custom_days} days"
#         else:
#             frequency = freq_choice

#         st.markdown("✅ Confirm your subscription")
#         st.info(f"""
#                 **Email**: {subscriber_email or "Not provided"}  
#                 **Journals**: {', '.join(selected_journals) if selected_journals else "None selected"}  
#                 **Keywords**: {raw_keywords if raw_keywords else "None"}  
#                 **Frequency**: {frequency}
#                 """)

#         if st.button("📩 Confirm and Subscribe"):
#             if not subscriber_email:
#                 st.error("❌ Please provide an email address.")
#             elif not selected_journals:
#                 st.error("❌ Please select at least one journal.")
#             else:
#                 formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)]
#                 csv_bytes = df.to_csv(index=False).encode("utf-8") if "df" in locals() and not df.empty else generate_placeholder_csv()

#                 result = store_user_subscription(
#                     email=subscriber_email,
#                     journals=formatted_journals,
#                     keywords=raw_keywords,
#                     start_date=start_date,
#                     end_date=end_date,
#                     frequency=frequency,
#                 )
#                 st.success(f"📬 Subscribed! You'll receive {frequency} updates at {subscriber_email}.")
#                 st.write("🛠️ Supabase insert result:", result)

#                 if result["status"] == "success":
#                     unsubscribe_token = result["unsubscribe_token"]
#                     unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"

#                     email_body = f"""Hi {subscriber_email},

#                 You have successfully subscribed to automatic PubMed updates.

#                 📘 Journals: {', '.join(formatted_journals)}
#                 🔑 Keywords: {raw_keywords or 'None'}
#                 🔁 Frequency: {frequency}
#                 📅 Date Range: {start_date} to {end_date}

#                 If you wish to unsubscribe, click the link below:
#                 🔓 {unsubscribe_link}

#                 – PubMed Tracker Team
#                     """
                    
#                     try:
#                         send_email(
#                             to_email=subscriber_email,
#                             subject="📬 Journal Tracker: Subscription Confirmed",
#                             body=email_body
#                         )
#                         st.success("✅ A confirmation email has been sent.")
#                     except Exception as e:
#                         st.warning(f"⚠️ Subscription saved, but email failed: {e}")


# import streamlit as st
# from tracking_main import (
#     fetch_pubmed_articles_by_date,
#     fetch_preprints,
#     load_pubmed_journal_abbreviations,
#     format_boolean_keywords_for_pubmed,
#     build_pubmed_query,
#     generate_placeholder_csv,
#     merge_and_highlight_articles,
#     compile_keyword_filter)

# import pandas as pd
# import os
# from store_subscription import store_user_subscription

# from datetime import datetime, timedelta
# from dotenv import load_dotenv
# from itsdangerous import URLSafeSerializer, BadSignature

# from email_dispatcher import send_email

# load_dotenv()

# EMAIL_ADDRESS = st.secrets["EMAIL_ADDRESS"]
# EMAIL_PASSWORD = st.secrets["EMAIL_PASSWORD"]

# SUPABASE_URL = os.getenv("SUPABASE_URL")
# SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
# serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

# BASE_URL = "https://journaltracker.streamlit.app"

# # Streamlit app configuration
# st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
# st.title("📚 PubMed Journal Tracker")

# # Check for unsubscribe token in the URL
# if 'token' in st.query_params:
#     token = st.query_params['token']  # Get the token from the query parameters
    
#     # Call the unsubscribe handling function
#     from unsubscribe import handle_unsubscribe  # Make sure this imports the unsubscribe logic
#     handle_unsubscribe(token)  # Call the function in unsubscribe.py
    
# else:
#     # If no token is found, display the main application interface

#     st.markdown("""
#     Use this tool to search PubMed by journal, date range, and keywords.  
#     You can also subscribe to automatic updates.
#     """)

#     journal_dict = load_pubmed_journal_abbreviations()
#     full_to_abbrev = {v: k for k, v in journal_dict.items()}
#     journal_options = list(full_to_abbrev.keys())

#     email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance.")
    
#     # Journal selection
#     selected_journals = st.multiselect("📘 Select journal(s) (Optional):", options=journal_options)
    
#     # Include preprints checkbox moved here
#     include_preprints = st.checkbox("📑 Include preprints", help="Currently supports bioRxiv and medRxiv.")

#     date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
#     today = datetime.today().date()

#     if date_option == "Custom":
#         col1, col2 = st.columns(2)
#         with col1:
#             start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
#         with col2:
#             end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
#     else:
#         days = {"Past Week": 7, "Past Month": 30, "Past Year": 365}[date_option]
#         start_date = str(today - timedelta(days=days))
#         end_date = str(today)

#     raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
#         help="""🔎 **Search Tips**  
#         - Use **AND**, **OR**, **NOT**  
#         - Wrap phrases in **parentheses**: *(cadmium exposure)*  
#         - Wildcards: `metagenom*`, `wom?n`
#         """)

#     # Manual Search Button
#     if st.button("🔍 Search"):
#         try:
#             # Modified validation - allow search if either journals are selected OR preprints are included
#             if not selected_journals and not include_preprints:
#                 st.error("❌ Please select at least one journal or include preprints.")
#                 st.stop()

#             # Only process journal validation if journals are selected
#             formatted_journals = []
#             if selected_journals:
#                 formatted_journals = [full_to_abbrev.get(j) for j in selected_journals if j in full_to_abbrev]
#                 if not formatted_journals:
#                     st.error("❌ Invalid journal names.")
#                     st.stop()

#             # Validate keywords early if provided
#             if raw_keywords.strip():
#                 if raw_keywords.count("(") != raw_keywords.count(")"):
#                     st.warning("⚠️ Unbalanced parentheses.")
#                     st.stop()
#                 user_keywords = raw_keywords.strip()
#                 compiled_filter = compile_keyword_filter(raw_keywords)
#                 pubmed_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
#                 keywords = pubmed_keywords
#             else:
#                 keywords = None

#             # Create placeholders for status updates
#             status_placeholder = st.empty()
#             progress_bar = st.progress(0)
#             details_container = st.container()
            
#             # Show search details in collapsible sections
#             with details_container:
#                 with st.expander("🔍 Search Details", expanded=False):
#                     st.write(f"**Selected Journals:** {', '.join(selected_journals) if selected_journals else 'None'}")
#                     st.write(f"**Date Range:** {start_date} to {end_date}")
#                     st.write(f"**Keywords:** {raw_keywords if raw_keywords else 'None'}")
#                     st.write(f"**Include Preprints:** {'Yes' if include_preprints else 'No'}")
                    
#                     # Only show PubMed query if journals are selected
#                     if keywords and selected_journals:
#                         query_preview = build_pubmed_query(
#                             journal=selected_journals[0],
#                             start_date=start_date,
#                             end_date=end_date,
#                             keywords=keywords
#                         )
#                         st.caption("📄 PubMed query example (first journal):")
#                         st.code(query_preview, language="none")

#             # Initialize search
#             status_placeholder.info("🚀 Starting search...")
#             progress_bar.progress(10)
            
#             all_articles = []
#             total_journals = len(selected_journals) if selected_journals else 0
#             preprint_servers = ["biorxiv", "medrxiv"] if include_preprints else []
#             total_sources = total_journals + len(preprint_servers)
            
#             current_step = 0
            
#             # Search PubMed journals (only if journals are selected)
#             if selected_journals:
#                 for i, journal in enumerate(selected_journals):
#                     current_step += 1
#                     status_placeholder.info(f"🔍 Searching {journal}... ({current_step}/{total_sources})")
#                     progress_bar.progress(int((current_step / total_sources) * 80))
                    
#                     try:
#                         articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, pubmed_keywords)
#                         for article in articles:
#                             article["Journal"] = journal
#                         all_articles.extend(articles)
                        
#                         # Show intermediate results
#                         if articles:
#                             status_placeholder.success(f"✅ Found {len(articles)} articles in {journal}")
#                         else:
#                             status_placeholder.info(f"📭 No articles found in {journal}")
                            
#                     except Exception as e:
#                         st.error(f"❌ Error searching {journal}: {str(e)}")
#                         continue

#             # Search preprints if requested
#             if include_preprints:
#                 for server in preprint_servers:
#                     current_step += 1
#                     status_placeholder.info(f"🔍 Searching {server} preprints... ({current_step}/{total_sources})")
#                     progress_bar.progress(int((current_step / total_sources) * 80))
                    
#                     try:
#                         preprints = fetch_preprints(
#                             server=server,
#                             start_date=start_date,
#                             end_date=end_date,
#                             keywords=raw_keywords
#                         )
#                         for article in preprints:
#                             article["Journal"] = server
#                             article["Source"] = "Preprint"
#                         all_articles.extend(preprints)
                        
#                         # Show intermediate results
#                         if preprints:
#                             status_placeholder.success(f"✅ Found {len(preprints)} preprints in {server}")
#                         else:
#                             status_placeholder.info(f"📭 No preprints found in {server}")
                            
#                     except Exception as e:
#                         st.error(f"❌ Error searching {server}: {str(e)}")
#                         continue

#             # Finalize results
#             progress_bar.progress(90)
#             status_placeholder.info("📊 Processing results...")

#             if all_articles:
#                 # Add Source field for PubMed articles
#                 for article in all_articles:
#                     if "Source" not in article:
#                         article["Source"] = "PubMed"

#                 merged = merge_and_highlight_articles(all_articles, [], raw_keywords)
#                 df = pd.DataFrame(merged)
                
#                 progress_bar.progress(100)
#                 status_placeholder.success(f"🎉 Search completed! Found {len(df)} articles total.")
                
#                 # Show results summary
#                 st.success(f"✅ Found {len(df)} articles across all sources.")
                
#                 # Results breakdown
#                 with st.expander("📊 Results Breakdown", expanded=True):
#                     if "Source" in df.columns:
#                         source_counts = df["Source"].value_counts()
#                         for source, count in source_counts.items():
#                             st.write(f"**{source}:** {count} articles")
                    
#                     if "Journal" in df.columns:
#                         journal_counts = df["Journal"].value_counts()
#                         st.write("\n**By Journal/Source:**")
#                         for journal, count in journal_counts.items():
#                             st.write(f"- {journal}: {count} articles")

#                 # Download button
#                 csv_data = df.to_csv(index=False).encode("utf-8")
#                 st.download_button(
#                     label="📥 Download Results as CSV",
#                     data=csv_data,
#                     file_name=f"search_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
#                     mime="text/csv"
#                 )
                
#                 # Show sample results
#                 with st.expander("👀 Preview Results", expanded=False):
#                     st.dataframe(df.head(10))
                
#             else:
#                 progress_bar.progress(100)
#                 status_placeholder.warning("📭 No articles found matching your criteria.")
                
#                 # Show suggestions
#                 st.info("""
#                 **No results found. Try:**
#                 - Expanding your date range
#                 - Using broader keywords
#                 - Checking different journals
#                 - Including preprints
#                 """)

#         except Exception as e:
#             st.error(f"❌ Unexpected error: {e}")
#             # Clear progress indicators on error
#             if 'status_placeholder' in locals():
#                 status_placeholder.empty()
#             if 'progress_bar' in locals():
#                 progress_bar.empty()

#     # --- Subscribe toggle ---
#     subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

#     # Subscription section only renders when checked
#     if subscribe:
#         col1, col2 = st.columns(2)
#         with col1:
#             freq_choice = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"])
#         with col2:
#             subscriber_email = st.text_input("📧 Email to receive updates")

#         if freq_choice == "custom":
#             custom_days = st.number_input("🔧 Custom interval (days):", min_value=1, step=1)
#             frequency = f"every {custom_days} days"
#         else:
#             frequency = freq_choice

#         st.markdown("✅ Confirm your subscription")
#         st.info(f"""
#                 **Email**: {subscriber_email or "Not provided"}  
#                 **Journals**: {', '.join(selected_journals) if selected_journals else "None selected"}  
#                 **Keywords**: {raw_keywords if raw_keywords else "None"}  
#                 **Frequency**: {frequency}
#                 """)

#         if st.button("📩 Confirm and Subscribe"):
#             if not subscriber_email:
#                 st.error("❌ Please provide an email address.")
#             elif not selected_journals and not include_preprints:
#                 st.error("❌ Please select at least one journal or include preprints.")
#             else:
#                 formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)] if selected_journals else []
#                 csv_bytes = df.to_csv(index=False).encode("utf-8") if "df" in locals() and not df.empty else generate_placeholder_csv()

#                 result = store_user_subscription(
#                     email=subscriber_email,
#                     journals=formatted_journals,
#                     keywords=raw_keywords,
#                     start_date=start_date,
#                     end_date=end_date,
#                     frequency=frequency,
#                 )
#                 st.success(f"📬 Subscribed! You'll receive {frequency} updates at {subscriber_email}.")
#                 st.write("🛠️ Supabase insert result:", result)

#                 if result["status"] == "success":
#                     unsubscribe_token = result["unsubscribe_token"]
#                     unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"

#                     email_body = f"""Hi {subscriber_email},

#                 You have successfully subscribed to automatic PubMed updates.

#                 📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'Preprints only'}
#                 🔑 Keywords: {raw_keywords or 'None'}
#                 🔁 Frequency: {frequency}
#                 📅 Date Range: {start_date} to {end_date}

#                 If you wish to unsubscribe, click the link below:
#                 🔓 {unsubscribe_link}

#                 – PubMed Tracker Team
#                     """
                    
#                     try:
#                         send_email(
#                             to_email=subscriber_email,
#                             subject="📬 Journal Tracker: Subscription Confirmed",
#                             body=email_body
#                         )
#                         st.success("✅ A confirmation email has been sent.")
#                     except Exception as e:
#                         st.warning(f"⚠️ Subscription saved, but email failed: {e}")



import streamlit as st
from tracking_main import (
    fetch_pubmed_articles_by_date,
    fetch_preprints,
    load_pubmed_journal_abbreviations,
    format_boolean_keywords_for_pubmed,
    build_pubmed_query,
    generate_placeholder_csv,
    merge_and_highlight_articles,
    compile_keyword_filter,
    standardize_date_format,
    standardize_doi_format
    )

import pandas as pd
import os
from store_subscription import store_user_subscription

from datetime import datetime, timedelta
from dotenv import load_dotenv
from itsdangerous import URLSafeSerializer, BadSignature

from email_dispatcher import send_email

load_dotenv()

EMAIL_ADDRESS = st.secrets["EMAIL_ADDRESS"]
EMAIL_PASSWORD = st.secrets["EMAIL_PASSWORD"]

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

BASE_URL = "https://journaltracker.streamlit.app"

# Streamlit app configuration
st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
st.title("📚 PubMed Journal Tracker")

# Check for unsubscribe token in the URL
if 'token' in st.query_params:
    token = st.query_params['token']  # Get the token from the query parameters
    
    # Call the unsubscribe handling function
    from unsubscribe import handle_unsubscribe  # Make sure this imports the unsubscribe logic
    handle_unsubscribe(token)  # Call the function in unsubscribe.py
    
else:
    # If no token is found, display the main application interface

    st.markdown("""
    Use this tool to search PubMed by journal, date range, and keywords.  
    You can also subscribe to automatic updates.
    """)

    journal_dict = load_pubmed_journal_abbreviations()
    full_to_abbrev = {v: k for k, v in journal_dict.items()}
    journal_options = list(full_to_abbrev.keys())

    email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance.")
    
    # Journal selection
    selected_journals = st.multiselect("📘 Select journal(s) (Optional):", options=journal_options)
    
    # Include preprints checkbox moved here
    include_preprints = st.checkbox("📑 Include preprints", help="Currently supports bioRxiv and medRxiv.")

    date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
    today = datetime.today().date()

    if date_option == "Custom":
        col1, col2 = st.columns(2)
        with col1:
            start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
        with col2:
            end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
    else:
        days = {"Past Week": 7, "Past Month": 30, "Past Year": 365}[date_option]
        start_date = str(today - timedelta(days=days))
        end_date = str(today)

    raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
        help="""🔎 **Search Tips**  
        - Use **AND**, **OR**, **NOT**  
        - Wrap phrases in **parentheses**: *(cadmium exposure)*  
        - Wildcards: `metagenom*`, `wom?n`
        """)

    # Manual Search Button
    if st.button("🔍 Search"):
        try:
            # Modified validation - allow search if either journals are selected OR preprints are included
            if not selected_journals and not include_preprints:
                st.error("❌ Please select at least one journal or include preprints.")
                st.stop()

            # Only process journal validation if journals are selected
            formatted_journals = []
            if selected_journals:
                formatted_journals = [full_to_abbrev.get(j) for j in selected_journals if j in full_to_abbrev]
                if not formatted_journals:
                    st.error("❌ Invalid journal names.")
                    st.stop()

            # Validate keywords early if provided
            if raw_keywords.strip():
                if raw_keywords.count("(") != raw_keywords.count(")"):
                    st.warning("⚠️ Unbalanced parentheses.")
                    st.stop()
                user_keywords = raw_keywords.strip()
                compiled_filter = compile_keyword_filter(raw_keywords)
                pubmed_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
                keywords = pubmed_keywords
            else:
                keywords = None

            # Create placeholders for status updates
            status_placeholder = st.empty()
            progress_bar = st.progress(0)
            
            # Initialize search
            status_placeholder.info("🚀 Starting search...")
            progress_bar.progress(10)
            
            all_articles = []
            total_journals = len(selected_journals) if selected_journals else 0
            preprint_servers = ["biorxiv", "medrxiv"] if include_preprints else []
            total_sources = total_journals + len(preprint_servers)
            
            current_step = 0
            
            # Search PubMed journals (only if journals are selected)
            if selected_journals:
                for i, journal in enumerate(selected_journals):
                    current_step += 1
                    status_placeholder.info(f"🔍 Searching {journal}... ({current_step}/{total_sources})")
                    progress_bar.progress(int((current_step / total_sources) * 80))
                    
                    try:
                        articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, pubmed_keywords)
                        for article in articles:
                            article["Journal"] = journal
                        all_articles.extend(articles)
                        
                        # Show intermediate results
                        if articles:
                            status_placeholder.success(f"✅ Found {len(articles)} articles in {journal}")
                        else:
                            status_placeholder.info(f"📭 No articles found in {journal}")
                            
                    except Exception as e:
                        st.error(f"❌ Error searching {journal}: {str(e)}")
                        continue

            # Search preprints if requested
            if include_preprints:
                for server in preprint_servers:
                    current_step += 1
                    status_placeholder.info(f"🔍 Searching {server} preprints... ({current_step}/{total_sources})")
                    progress_bar.progress(int((current_step / total_sources) * 80))
                    
                    try:
                        preprints = fetch_preprints(
                            server=server,
                            start_date=start_date,
                            end_date=end_date,
                            keywords=raw_keywords
                        )
                        for article in preprints:
                            article["Journal"] = server
                            article["Source"] = "Preprint"
                        all_articles.extend(preprints)
                        
                        # Show intermediate results
                        if preprints:
                            status_placeholder.success(f"✅ Found {len(preprints)} preprints in {server}")
                        else:
                            status_placeholder.info(f"📭 No preprints found in {server}")
                            
                    except Exception as e:
                        st.error(f"❌ Error searching {server}: {str(e)}")
                        continue

            # Finalize results
            progress_bar.progress(90)
            status_placeholder.info("📊 Processing results...")

            if all_articles:
                # Add Source field for PubMed articles
                for article in all_articles:
                    if "Source" not in article:
                        article["Source"] = "PubMed"

                # Standardize date formats before processing
                all_articles = standardize_date_format(all_articles)
                all_articles = standardize_doi_format(all_articles)
                
                merged = merge_and_highlight_articles(all_articles, [], raw_keywords)
                df = pd.DataFrame(merged)
                
                progress_bar.progress(100)
                status_placeholder.success(f"🎉 Search completed! Found {len(df)} articles total.")
                
                # Combined search details and results in a single expander
                with st.expander("📊 Search Summary & Results", expanded=True):
                    # Search parameters in a compact format
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**📘 Journals:** {', '.join(selected_journals[:2]) + ('...' if len(selected_journals) > 2 else '') if selected_journals else 'None'}")
                        st.write(f"**📅 Date Range:** {start_date} to {end_date}")
                    with col2:
                        st.write(f"**🔑 Keywords:** {raw_keywords[:30] + '...' if len(raw_keywords) > 30 else raw_keywords if raw_keywords else 'None'}")
                        st.write(f"**📑 Preprints:** {'Yes' if include_preprints else 'No'}")
                    
                    # Results breakdown
                    st.write("---")
                    st.write("**📈 Results by Source:**")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        if "Source" in df.columns:
                            source_counts = df["Source"].value_counts()
                            for source, count in source_counts.items():
                                st.write(f"• **{source}:** {count} articles")
                    
                    with col2:
                        if "Journal" in df.columns:
                            journal_counts = df["Journal"].value_counts()
                            st.write("**By Journal/Platform:**")
                            for journal, count in journal_counts.head(5).items():  # Show top 5
                                st.write(f"• {journal}: {count}")
                            if len(journal_counts) > 5:
                                st.write(f"• ... and {len(journal_counts) - 5} more")
                    
                    # Show PubMed query if applicable
                    if keywords and selected_journals:
                        with st.expander("🔍 Technical Details", expanded=False):
                            query_preview = build_pubmed_query(
                                journal=selected_journals[0],
                                start_date=start_date,
                                end_date=end_date,
                                keywords=keywords
                            )
                            st.caption("🔧 Actual PubMed API query used for search:")
                            st.code(query_preview, language="none")
                
                # Show sample results first
                with st.expander("👀 Preview Results", expanded=False):
                    st.dataframe(df.head(10))
                
                # Download button moved below preview
                csv_data = df.to_csv(index=False).encode("utf-8")
                st.download_button(
                    label="📥 Download Results as CSV",
                    data=csv_data,
                    file_name=f"search_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv"
                )
                
            else:
                progress_bar.progress(100)
                status_placeholder.warning("📭 No articles found matching your criteria.")
                
                # Show suggestions
                st.info("""
                **No results found. Try:**
                - Expanding your date range
                - Using broader keywords
                - Checking different journals
                - Including preprints
                """)

        except Exception as e:
            st.error(f"❌ Unexpected error: {e}")
            # Clear progress indicators on error
            if 'status_placeholder' in locals():
                status_placeholder.empty()
            if 'progress_bar' in locals():
                progress_bar.empty()

    # --- Subscribe toggle ---
    subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

    # Subscription section only renders when checked
    if subscribe:
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

        st.markdown("✅ Confirm your subscription")
        st.info(f"""
                **Email**: {subscriber_email or "Not provided"}  
                **Journals**: {', '.join(selected_journals) if selected_journals else "None selected"}  
                **Keywords**: {raw_keywords if raw_keywords else "None"}  
                **Frequency**: {frequency}
                """)

        if st.button("📩 Confirm and Subscribe"):
            if not subscriber_email:
                st.error("❌ Please provide an email address.")
            elif not selected_journals and not include_preprints:
                st.error("❌ Please select at least one journal or include preprints.")
            else:
                formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)] if selected_journals else []
                csv_bytes = df.to_csv(index=False).encode("utf-8") if "df" in locals() and not df.empty else generate_placeholder_csv()

                result = store_user_subscription(
                    email=subscriber_email,
                    journals=formatted_journals,
                    keywords=raw_keywords,
                    start_date=start_date,
                    end_date=end_date,
                    frequency=frequency,
                )
                st.success(f"📬 Subscribed! You'll receive {frequency} updates at {subscriber_email}.")
                st.write("🛠️ Supabase insert result:", result)

                if result["status"] == "success":
                    unsubscribe_token = result["unsubscribe_token"]
                    unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"

                    email_body = f"""Hi {subscriber_email},

                You have successfully subscribed to automatic PubMed updates.

                📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'Preprints only'}
                🔑 Keywords: {raw_keywords or 'None'}
                🔁 Frequency: {frequency}
                📅 Date Range: {start_date} to {end_date}

                If you wish to unsubscribe, click the link below:
                🔓 {unsubscribe_link}

                – PubMed Tracker Team
                    """
                    
                    try:
                        send_email(
                            to_email=subscriber_email,
                            subject="📬 Journal Tracker: Subscription Confirmed",
                            body=email_body
                        )
                        st.success("✅ A confirmation email has been sent.")
                    except Exception as e:
                        st.warning(f"⚠️ Subscription saved, but email failed: {e}")