# # # app.py

# # import streamlit as st
# # from tracking_main import (
# #     fetch_pubmed_articles_by_date,
# #     format_journal_abbreviation,
# #     fetch_preprints,
# #     load_pubmed_journal_abbreviations,
# #     format_boolean_keywords_for_pubmed,
# #     build_pubmed_query,
# #     generate_placeholder_csv,
# #     merge_and_highlight_articles,
# #     compile_keyword_filter,
# #     standardize_date_format,
# #     standardize_doi_format
# #     )

# # import pandas as pd
# # import os
# # import io
# # import logging
# # from store_subscription import store_user_subscription, generate_unsubscribe_token

# # from RateLimit import PubMedRateLimit

# # from datetime import datetime, timedelta
# # from dotenv import load_dotenv
# # from itsdangerous import URLSafeSerializer, URLSafeTimedSerializer, BadSignature

# # from email_dispatcher import send_email, generate_download_token, get_next_update_timeframe

# # load_dotenv()

# # # UPDATED: Use new domain-based email configuration
# # BREVO_API_KEY = st.secrets["BREVO_API_KEY"]
# # BREVO_SENDER_EMAIL = st.secrets["BREVO_SENDER_EMAIL"]
# # BREVO_SENDER_NAME = st.secrets["BREVO_SENDER_NAME"]

# # SUPABASE_URL = os.getenv("SUPABASE_URL")
# # SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# # BASE_URL = "https://journaltracker.streamlit.app"

# # # Initialize rate limiter
# # if 'rate_limiter' not in st.session_state:
# #     st.session_state.rate_limiter = PubMedRateLimit()

# # rate_limiter = st.session_state.rate_limiter

# # # function to validate journal selection
# # def validate_and_format_journals(selected_journals):
# #     """
# #     Validate and format selected journals using the enhanced logic
# #     """
# #     journal_dict = load_pubmed_journal_abbreviations()
# #     formatted_journals = []
    
# #     for journal in selected_journals:
# #         try:
# #             formatted_journal = format_journal_abbreviation(journal, journal_dict)
# #             formatted_journals.append(formatted_journal)
# #         except ValueError as e:
# #             st.error(f"Journal validation error for '{journal}': {str(e)}")
# #             return None
    
# #     return formatted_journals

# # # Search Execution Function
# # def execute_subscription_search(journals, keywords, include_preprints, frequency):
# #     """Execute search based on subscription parameters"""
# #     # Calculate date range based on frequency
# #     today = datetime.today().date()
    
# #     if frequency == "weekly":
# #         days_back = 7
# #     elif frequency == "monthly":
# #         days_back = 30
# #     elif frequency.startswith("every"):
# #         # Extract number from "every X days"
# #         days_back = int(frequency.split()[1])
# #     else:
# #         days_back = 7  # Default fallback
    
# #     start_date = str(today - timedelta(days=days_back))
# #     end_date = str(today)
    
# #     all_articles = []
    
# #     try:
# #         # Search PubMed journals
# #         if journals:
# #             for journal in journals:
# #                 articles = fetch_pubmed_articles_by_date(
# #                     journal, start_date, end_date, 
# #                     format_boolean_keywords_for_pubmed(keywords) if keywords else None
# #                 )
# #                 for article in articles:
# #                     article["Journal"] = journal
# #                     article["Source"] = "PubMed"
# #                 all_articles.extend(articles)
        
# #         # Search preprints
# #         if include_preprints:
# #             for server in ["biorxiv", "medrxiv"]:
# #                 preprints = fetch_preprints(
# #                     server=server,
# #                     start_date=start_date,
# #                     end_date=end_date,
# #                     keywords=keywords
# #                 )
# #                 for article in preprints:
# #                     article["Journal"] = server
# #                     article["Source"] = "Preprint"
# #                 all_articles.extend(preprints)
        
# #         if all_articles:
# #             # Process results
# #             all_articles = standardize_date_format(all_articles)
# #             all_articles = standardize_doi_format(all_articles)
# #             merged = merge_and_highlight_articles(all_articles, [], keywords)
# #             return pd.DataFrame(merged)
# #         else:
# #             return pd.DataFrame()
            
# #     except Exception as e:
# #         logging.error(f"Error in subscription search: {e}")
# #         return pd.DataFrame()
    

# # # Streamlit app configuration
# # st.set_page_config(page_title="PubMed Journal Tracker", layout="centered")
# # st.title("📚 PubMed Journal Tracker")

# # # Check for unsubscribe token in the URL
# # if 'token' in st.query_params:
# #     token = st.query_params['token']
# #     action = st.query_params.get('action', 'unsubscribe')  # Default to unsubscribe
    
# #     if action == 'download':
# #         # Handle CSV download
# #         from app_csv_downloader import handle_download
# #         handle_download(token)
# #     else:
# #         # Handle unsubscribe
# #         from app_unsubscribe import handle_unsubscribe
# #         handle_unsubscribe(token)

# # else:
# #     # If no token is found, display the main application interface

# #     st.markdown("""
# #     Use this tool to search PubMed by journal, date range, and keywords.  
# #     You can also subscribe to automatic updates.
# #     """)

# #     journal_dict = load_pubmed_journal_abbreviations()
# #     full_to_abbrev = {v: k for k, v in journal_dict.items()}

# #     # Quick fix: Add missing journals
# #     missing_journals = ["The Lancet", "BMJ"]
# #     for journal in missing_journals:
# #         if journal not in full_to_abbrev:
# #             full_to_abbrev[journal] = journal

# #     journal_options = list(full_to_abbrev.keys())
# #     journal_options.sort()

# #     email = st.text_input("📧 Enter your email (Optional):", help="Used for NCBI API compliance.")
    
# #     # Journal selection
# #     selected_journals = st.multiselect("📘 Select journal(s) (Optional):", options=journal_options)

# #     # Include preprints checkbox moved here
# #     include_preprints = st.checkbox("📑 Include preprints", help="Currently supports bioRxiv and medRxiv.")

# #     date_option = st.selectbox("📅 Select date range:", ["Past Week", "Past Month", "Past Year", "Custom"])
# #     today = datetime.today().date()

# #     if date_option == "Custom":
# #         col1, col2 = st.columns(2)
# #         with col1:
# #             start_date = st.text_input("Start date (YYYY-MM or YYYY-MM-DD):")
# #         with col2:
# #             end_date = st.text_input("End date (YYYY-MM or YYYY-MM-DD):")
# #     else:
# #         days = {"Past Week": 7, "Past Month": 30, "Past Year": 365}[date_option]
# #         start_date = str(today - timedelta(days=days))
# #         end_date = str(today)

# #     raw_keywords = st.text_area("❓ Enter your search keyword (Optional):", height=100,
# #         help="""🔎 **Search Tips**  
# #         - Use **AND**, **OR**, **NOT**  
# #         - Wrap phrases in **parentheses**: *(cadmium exposure)*  
# #         - Wildcards: `metagenom*`, `wom?n`
# #         """)

# #    # Enhanced search section with detailed progress tracking
# #     if st.button("🔍 Search"):
# #         # Rate limiting check first
# #         if not rate_limiter.can_make_request():
# #             st.stop()  # This stops execution if rate limited

# #         try:
# #             # Modified validation - allow search if either journals are selected OR preprints are included
# #             if not selected_journals and not include_preprints:
# #                 st.error("❌ Please select at least one journal or include preprints.")
# #                 st.stop()

# #             # Only process journal validation if journals are selected
# #             if selected_journals:
# #                 formatted_journals = validate_and_format_journals(selected_journals)
# #                 if formatted_journals is None:  # Validation failed
# #                     st.stop()
# #             else:
# #                 formatted_journals = []

# #             # Validation
# #             if not formatted_journals and not include_preprints:
# #                 st.error("❌ Please select at least one journal or enable preprints.")
# #                 st.stop()

# #             # Validate keywords early if provided
# #             if raw_keywords.strip():
# #                 if raw_keywords.count("(") != raw_keywords.count(")"):
# #                     st.warning("⚠️ Unbalanced parentheses.")
# #                     st.stop()
# #                 user_keywords = raw_keywords.strip()
# #                 compiled_filter = compile_keyword_filter(raw_keywords)
# #                 pubmed_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
# #                 keywords = pubmed_keywords
# #             else:
# #                 keywords = None

# #             # Create placeholders for status updates
# #             status_placeholder = st.empty()
# #             progress_bar = st.progress(0)
            
# #             # Initialize search
# #             status_placeholder.info("🚀 Starting search...")
# #             progress_bar.progress(10)
            
# #             all_articles = []
# #             total_journals = len(selected_journals) if selected_journals else 0
# #             preprint_servers = ["biorxiv", "medrxiv"] if include_preprints else []
# #             total_sources = total_journals + len(preprint_servers)
            
# #             current_step = 0
            
# #             # Search PubMed journals (only if journals are selected) - WITH RATE LIMITING
# #             if selected_journals:
# #                 for i, journal in enumerate(selected_journals):
# #                     current_step += 1
# #                     status_placeholder.info(f"🔍 Searching {journal}... ({current_step}/{total_sources})")
# #                     progress_bar.progress(int((current_step / total_sources) * 80))
                    
# #                     try:
# #                         # Use the safe search method from rate limiter
# #                         query = build_pubmed_query(journal, start_date, end_date, keywords)
# #                         search_results = rate_limiter.safe_pubmed_search(query, max_results=100)
                        
# #                         if search_results and search_results.get("IdList"):
# #                             # Use safe fetch method for getting article details
# #                             articles_xml = rate_limiter.safe_pubmed_fetch(search_results["IdList"])
                            
# #                             if articles_xml:
# #                                 # Your existing article processing logic
# #                                 # You'll need to adapt this part to work with the XML data
# #                                 # For now, let's use your existing function but add rate limiting
# #                                 articles = fetch_pubmed_articles_by_date(journal, start_date, end_date, keywords, rate_limiter)

# #                                 for article in articles:
# #                                     article["Journal"] = journal
# #                                 all_articles.extend(articles)
                                
# #                                 # Show intermediate results
# #                                 if articles:
# #                                     status_placeholder.success(f"✅ Found {len(articles)} articles in {journal}")
# #                                 else:
# #                                     status_placeholder.info(f"📭 No articles found in {journal}")
# #                             else:
# #                                 status_placeholder.warning(f"⚠️ Failed to fetch details for {journal}")
# #                         else:
# #                             status_placeholder.info(f"📭 No articles found in {journal}")
                                
# #                     except Exception as e:
# #                         st.error(f"❌ Error searching {journal}: {str(e)}")
# #                         continue

# #             # Search preprints if requested (no rate limiting needed for preprints)
# #             if include_preprints:
# #                 for server in preprint_servers:
# #                     current_step += 1
# #                     status_placeholder.info(f"🔍 Searching {server} preprints... ({current_step}/{total_sources})")
# #                     progress_bar.progress(int((current_step / total_sources) * 80))
                    
# #                     try:
# #                         preprints = fetch_preprints(
# #                             server=server,
# #                             start_date=start_date,
# #                             end_date=end_date,
# #                             keywords=raw_keywords
# #                         )
# #                         for article in preprints:
# #                             article["Journal"] = server
# #                             article["Source"] = "Preprint"
# #                         all_articles.extend(preprints)
                        
# #                         # Show intermediate results
# #                         if preprints:
# #                             status_placeholder.success(f"✅ Found {len(preprints)} preprints in {server}")
# #                         else:
# #                             status_placeholder.info(f"📭 No preprints found in {server}")
                            
# #                     except Exception as e:
# #                         st.error(f"❌ Error searching {server}: {str(e)}")
# #                         continue

# #             # Rest of your existing code remains the same...
# #             # Finalize results
# #             progress_bar.progress(90)
# #             status_placeholder.info("📊 Processing results...")

# #             if all_articles:
# #                 # Add Source field for PubMed articles
# #                 for article in all_articles:
# #                     if "Source" not in article:
# #                         article["Source"] = "PubMed"

# #                 # Standardize date formats before processing
# #                 all_articles = standardize_date_format(all_articles)
# #                 all_articles = standardize_doi_format(all_articles)
                
# #                 merged = merge_and_highlight_articles(all_articles, [], raw_keywords)
# #                 df = pd.DataFrame(merged)
                
# #                 progress_bar.progress(100)
# #                 status_placeholder.success(f"🎉 Search completed! Found {len(df)} articles total.")
                
# #                 # Display results
# #                 with st.expander("📊 Search Summary & Results", expanded=True):
# #                     # Search parameters in a compact format
# #                     col1, col2 = st.columns(2)
# #                     with col1:
# #                         st.write(f"**📘 Journals:** {', '.join(selected_journals[:2]) + ('...' if len(selected_journals) > 2 else '') if selected_journals else 'None'}")
# #                         st.write(f"**📅 Date Range:** {start_date} to {end_date}")
# #                     with col2:
# #                         st.write(f"**🔑 Keywords:** {raw_keywords[:30] + '...' if len(raw_keywords) > 30 else raw_keywords if raw_keywords else 'None'}")
# #                         st.write(f"**📑 Preprints:** {'Yes' if include_preprints else 'No'}")
                    
# #                     # Results breakdown
# #                     st.write("---")
# #                     st.write("**📈 Results by Source:**")
                    
# #                     col1, col2 = st.columns(2)
# #                     with col1:
# #                         if "Source" in df.columns:
# #                             source_counts = df["Source"].value_counts()
# #                             for source, count in source_counts.items():
# #                                 st.write(f"• **{source}:** {count} articles")
                    
# #                     with col2:
# #                         if "Journal" in df.columns:
# #                             journal_counts = df["Journal"].value_counts()
# #                             st.write("**By Journal/Platform:**")
# #                             for journal, count in journal_counts.head(5).items():  # Show top 5
# #                                 st.write(f"• {journal}: {count}")
# #                             if len(journal_counts) > 5:
# #                                 st.write(f"• ... and {len(journal_counts) - 5} more")
                    
# #                     # Show PubMed query if applicable
# #                     if keywords and selected_journals:
# #                         with st.expander("🔍 Technical Details", expanded=False):
# #                             query_preview = build_pubmed_query(
# #                                 journal=selected_journals[0],
# #                                 start_date=start_date,
# #                                 end_date=end_date,
# #                                 keywords=keywords
# #                             )
# #                             st.caption("🔧 Actual PubMed API query used for search:")
# #                             st.code(query_preview, language="none")
                
# #                 # Show sample results first
# #                 with st.expander("👀 Preview Results", expanded=False):
# #                     st.dataframe(df.head(10))
                
# #                 # Download button moved below preview
# #                 csv_data = df.to_csv(index=False).encode("utf-8")
# #                 st.download_button(
# #                     label="📥 Download Results as CSV",
# #                     data=csv_data,
# #                     file_name=f"JournalTracker_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
# #                     mime="text/csv"
# #                 )

# #             else:
# #                 progress_bar.progress(100)
# #                 status_placeholder.warning("📭 No articles found matching your criteria.")
                
# #                 # Show suggestions
# #                 st.info("""
# #                 **No results found. Try:**
# #                 - Expanding your date range
# #                 - Using broader keywords
# #                 - Checking different journals
# #                 - Including preprints
# #                 """)

# #         except Exception as e:
# #             st.error(f"❌ Unexpected error: {e}")
# #             # Clear progress indicators on error
# #             if 'status_placeholder' in locals():
# #                 status_placeholder.empty()
# #             if 'progress_bar' in locals():
# #                 progress_bar.empty()


# #     # --- Subscribe toggle ---
# #     subscribe = st.checkbox("📬 Subscribe to automatic updates", key="subscribe_toggle")

# #     # Subscription section only renders when checked
# #     if subscribe:
# #         col1, col2 = st.columns(2)
# #         with col1:
# #             freq_choice = st.selectbox("🔁 Update Frequency", ["weekly", "monthly", "custom"])
# #         with col2:
# #             subscriber_email = st.text_input("📧 Email to receive updates")

# #         if freq_choice == "custom":
# #             custom_days = st.number_input("🔧 Custom interval (days):", min_value=1, step=1)
# #             frequency = f"every {custom_days} days"
# #         else:
# #             frequency = freq_choice

# #         st.markdown("✅ Confirm your subscription")

# #         # source summary
# #         sources_summary = []
# #         if selected_journals:
# #             sources_summary.append(f"Journals: {', '.join(selected_journals)}")
# #         if include_preprints:
# #             sources_summary.append("Preprints: bioRxiv, medRxiv")

# #         sources_text = " | ".join(sources_summary) if sources_summary else "None selected"

# #         st.info(f"""
# #                 **Email**: {subscriber_email or "Not provided"}  
# #                 **Sources**: {sources_text}
# #                 **Keywords**: {raw_keywords if raw_keywords else "None"}  
# #                 **Frequency**: {frequency}
# #                 """)

# #         if st.button("📩 Confirm and Subscribe"):
# #             if not subscriber_email:
# #                 st.error("❌ Please provide an email address.")
# #             elif not selected_journals and not include_preprints:
# #                 st.error("❌ Please select at least one journal or include preprints.")
# #             else:
# #                 formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)] if selected_journals else []
                
# #                 # Execute search immediately for subscription
# #                 if 'df' in locals() and not df.empty:
# #                     # User just performed a search - use those results
# #                     csv_bytes = df.to_csv(index=False).encode("utf-8")
# #                     has_results = True
# #                 else:
# #                     # User is subscribing without searching - execute search now
# #                     search_results = execute_subscription_search(
# #                         journals=formatted_journals,
# #                         keywords=raw_keywords,
# #                         include_preprints=include_preprints,
# #                         frequency=frequency
# #                     )
# #                     if search_results is not None and not search_results.empty:
# #                         csv_bytes = search_results.to_csv(index=False).encode("utf-8")
# #                         has_results = True
# #                     else:
# #                         csv_bytes = None
# #                         has_results = False

# #                 # Store subscription first

# #                 result = store_user_subscription(
# #                     email=subscriber_email,
# #                     journals=formatted_journals,
# #                     keywords=raw_keywords,
# #                     # start_date=start_date, # might now need these two parameters for now
# #                     # end_date=end_date,
# #                     frequency=frequency,
# #                     include_preprints=include_preprints
# #                 )

# #                 if result["status"] == "success":

# #                     subscription_data = result["data"]  # Get the subscription data
# #                     subscription_id = subscription_data["id"]  # Extract the ID

# #                     # Generate proper unsubscribe token
# #                     unsubscribe_token = generate_unsubscribe_token(subscription_id)
                    
# #                     # Add debugging
# #                     st.write(f"🔧 Debug - Generated token (first 20 chars): {unsubscribe_token[:20] if unsubscribe_token else 'None'}...")    

# #                     if unsubscribe_token:
# #                         unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"
# #                     else:
# #                         st.error("❌ Failed to generate unsubscribe token")
# #                         st.stop()

# #                     # Build source description
# #                     source_list = []
# #                     if formatted_journals:
# #                         source_list.extend(formatted_journals)
# #                     if include_preprints:
# #                         source_list.extend(['bioRxiv', 'medRxiv'])
# #                     source_description = ', '.join(source_list) if source_list else 'No sources selected'

# #                     if has_results and csv_bytes is not None:
# #                         # Generate download token for results
# #                         download_token = generate_download_token(
# #                             csv_bytes, 
# #                             subscriber_email, 
# #                             subscription_id=subscription_id  # Pass the subscription ID
# #                         )
                        
# #                         if download_token:
# #                             download_link = f"{BASE_URL}?token={download_token}&action=download"
# #                             # ... rest of your email code
# #                         else:
# #                             st.error("❌ Failed to generate download link")
# #                             st.stop()

# #                         # Calculate result count
# #                         result_count = len(pd.read_csv(io.StringIO(csv_bytes.decode('utf-8'))))
                        
# #                         email_body = f"""Hi {subscriber_email},

# #                 You have successfully subscribed to automatic PubMed updates.

# #                 📊 SUBSCRIPTION DETAILS:
# #                 📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
# #                 📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
# #                 🔍 All Sources: {source_description}
# #                 🔑 Keywords: {raw_keywords or 'None'}
# #                 🔁 Frequency: {frequency}

# #                 📥 YOUR CURRENT RESULTS ({result_count} articles found):
# #                 Your search results are available for download (expires in 3 hours):
# #                 🔗 {download_link}

# #                 You will receive your next update in {get_next_update_timeframe(frequency)}.

# #                 🔓 UNSUBSCRIBE:
# #                 {unsubscribe_link}

# #                 – PubMed Tracker Team
# #                         """
                        
# #                         email_subject = f"📬 Journal Tracker: Subscription Confirmed ({result_count} results)"
# #                         success_message = f"✅ Subscription confirmed! {result_count} results sent to {subscriber_email}"
# #                     else:
# #                         # No results found
# #                         email_body = f"""Hi {subscriber_email},

# #                 You have successfully subscribed to automatic PubMed updates.

# #                 📊 SUBSCRIPTION DETAILS:
# #                 📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
# #                 📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
# #                 🔍 All Sources: {source_description}
# #                 🔑 Keywords: {raw_keywords or 'None'}
# #                 🔁 Frequency: {frequency}

# #                 📭 CURRENT SEARCH STATUS:
# #                 No articles found matching your criteria for the selected time period.

# #                 You will receive your next update in {get_next_update_timeframe(frequency)} (only if results are found).

# #                 🔓 UNSUBSCRIBE:
# #                 {unsubscribe_link}

# #                 – PubMed Tracker Team
# #                         """
                        
# #                         email_subject = "📬 Journal Tracker: Subscription Confirmed (No results)"
# #                         success_message = f"✅ Subscription confirmed! Confirmation email sent to {subscriber_email}"
                    
# #                     try:
# #                         # ✅ UPDATED: Pass email configuration explicitly
# #                         send_email(
# #                             to_email=subscriber_email,
# #                             subject=email_subject,
# #                             body=email_body,
# #                             sender_email=BREVO_SENDER_EMAIL,
# #                             sender_name=BREVO_SENDER_NAME,
# #                             api_key=BREVO_API_KEY
# #                         )
# #                         st.success(success_message)
# #                     except Exception as e:
# #                         st.error(f"❌ Subscription failed: {e}")
                
# #                 elif result["status"] == "error":
# #                     st.error(f"❌ {result['message']}")
# #                     if "Maximum of 3" in result["message"]:
# #                         st.info("💡 **Tip:** You can manage your existing subscriptions by using the unsubscribe link in any of your previous emails.")
                                
# #                 else:
# #                     st.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")

# # # Add Footer
# # st.markdown("---")
# # st.markdown(
# #     """
# #     <div style='text-align: center; color: #666; font-size: 0.8em; padding: 20px 0;'>
# #         Please report issues to <a href='https://github.com/rzhan186/Journal_tracker/issues' target='_blank'>GitHub</a>
# #     </div>
# #     """, 
# #     unsafe_allow_html=True
# # )

# # app.py

# import streamlit as st
# import pandas as pd
# from datetime import datetime, timedelta
# import time
# import re
# from tracking_main import (
#     fetch_pubmed_articles_by_date,
#     format_journal_abbreviation,
#     fetch_preprints,
#     load_pubmed_journal_abbreviations,
#     build_pubmed_query,
#     merge_and_highlight_articles,
#     standardize_date_format,
#     standardize_doi_format,
#     compile_keyword_filter,
#     format_boolean_keywords_for_pubmed,
#     execute_subscription_search,
#     OptimizedPubMedFetcher  # Add this import
# )
# from RateLimit import PubMedRateLimit
# from store_subscription import store_user_subscription, generate_unsubscribe_token
# from email_dispatcher import send_email, generate_download_token, get_next_update_timeframe
# import os

# # Configuration
# BASE_URL = os.getenv("BASE_URL", "https://journaltracker.streamlit.app/")
# BREVO_API_KEY = os.getenv("BREVO_API_KEY")
# BREVO_SENDER_EMAIL = os.getenv("BREVO_SENDER_EMAIL")
# BREVO_SENDER_NAME = os.getenv("BREVO_SENDER_NAME", "PubMed Journal Tracker")

# # Initialize rate limiter
# rate_limiter = PubMedRateLimit()

# if not os.getenv("NCBI_API_KEY"):
#     st.warning("⚠️ NCBI API key not configured. Performance may be limited.")

# # Page configuration
# st.set_page_config(
#     page_title="PubMed Journal Tracker",
#     page_icon="📚",
#     layout="wide",
#     initial_sidebar_state="expanded"
# )

# # Custom CSS for better styling  
# st.markdown("""  
# <style>  
#     .main-header {  
#         text-align: center;  
#         color: #2E86AB;  
#         margin-bottom: 2rem;  
#     }  
#     .search-container {  
#         background-color: #f8f9fa;  
#         padding: 1.5rem;  
#         border-radius: 10px;  
#         margin-bottom: 2rem;  
#         margin-top: 0rem;  
#     }  
#     .result-container {  
#         background-color: #ffffff;  
#         padding: 1rem;  
#         border-radius: 8px;  
#         border: 1px solid #e0e0e0;  
#     }  
#     .progress-container {  
#         background-color: #f0f2f6;  
#         padding: 1rem;  
#         border-radius: 8px;  
#         margin: 0rem 0 1rem 0;  
#     }  
#     .stats-card {  
#         background-color: #ffffff;  
#         padding: 1.5rem;  
#         border-radius: 10px;  
#         border: 1px solid #e0e0e0;  
#         margin-bottom: 1rem;  
#         text-align: center;  
#     }  
#     .feature-card {  
#         background-color: #f8f9fa;  
#         padding: 1.5rem;  
#         border-radius: 10px;  
#         border-left: 4px solid #2E86AB;  
#         margin-bottom: 1rem;  
#     }  
#     .info-card {  
#         background-color: #e8f4f8;  
#         padding: 1rem;  
#         border-radius: 8px;  
#         margin-bottom: 1rem;  
#     }  
#     .stProgress .st-bo {  
#         background-color: #2E86AB;  
#     }  
# </style>  
# """, unsafe_allow_html=True)  

# # Check for unsubscribe token in the URL  
# if 'token' in st.query_params:  
#     token = st.query_params['token']  
#     action = st.query_params.get('action', 'unsubscribe')  
    
#     if action == 'download':  
#         from app_csv_downloader import handle_download  
#         handle_download(token)  
#     else:  
#         from app_unsubscribe import handle_unsubscribe  
#         handle_unsubscribe(token)  
# else:  
#     # Main application interface  
#     st.markdown("<h1 class='main-header'>📚 PubMed Journal Tracker</h1>", unsafe_allow_html=True)  
    
#     # # Two column layout for main content  
#     # col1, col2 = st.columns([1, 1])  
    
#     # with col1:  
#     #     st.markdown("""  
#     #     <div class='info-card'>  
#     #         <h4>🔍 Search & Track Publications</h4>  
#     #         <p>Search PubMed by journal, date range, and keywords. Subscribe to automatic updates and never miss important research.</p>  
#     #     </div>  
#     #     """, unsafe_allow_html=True)  
    
#     # with col2:  
#     #     st.markdown("""  
#     #     <div class='feature-card'>  
#     #         <h4>✨ Features</h4>  
#     #         <ul>  
#     #             <li>🎯 Multi-journal search</li>  
#     #             <li>📑 Preprint integration</li>  
#     #             <li>🔔 Email notifications</li>  
#     #             <li>📊 Progress tracking</li>  
#     #             <li>📥 CSV export</li>  
#     #         </ul>  
#     #     </div>  
#     #     """, unsafe_allow_html=True)  

#     # Sidebar for search parameters  
#     with st.sidebar:  
#         st.header("🔍 Search Parameters")  
        
#         # Date range selection  
#         st.subheader("📅 Date Range")  
#         col1, col2 = st.columns(2)  
#         with col1:  
#             start_date = st.date_input(  
#                 "Start Date",  
#                 value=datetime.now() - timedelta(days=30),  
#                 help="Select the start date for your search"  
#             )  
#         with col2:  
#             end_date = st.date_input(  
#                 "End Date",  
#                 value=datetime.now(),  
#                 help="Select the end date for your search"  
#             )  
        
#         # Validate date range  
#         if start_date >= end_date:  
#             st.error("❌ Start date must be before end date")  
        
#         # Journal selection  
#         st.subheader("📘 Journal Selection")  
        
#         # Load journal abbreviations  
#         journal_dict = load_pubmed_journal_abbreviations()  
#         full_to_abbrev = {v: k for k, v in journal_dict.items()}  
        
#         # Add missing journals  
#         missing_journals = ["The Lancet", "BMJ"]  
#         for journal in missing_journals:  
#             if journal not in full_to_abbrev:  
#                 full_to_abbrev[journal] = journal  
        
#         journal_options = sorted(full_to_abbrev.keys())  
        
#         # Multi-select for journals  
#         selected_journals = st.multiselect(  
#             "Select Journals",  
#             options=journal_options,  
#             help="Choose one or more journals to search"  
#         )  
        
#         # Preprint option  
#         include_preprints = st.checkbox(  
#             "📑 Search Preprints (bioRxiv, medRxiv)",  
#             value=False,  
#             help="Include preprint servers in your search"  
#         )  
        
#         # Keywords input  
#         st.subheader("🔑 Keywords")  
#         raw_keywords = st.text_area(  
#             "Enter Keywords (optional)",  
#             help="Enter keywords or boolean expressions (e.g., 'cancer AND treatment', 'diabetes OR obesity')",  
#             placeholder="e.g., machine learning, AI AND healthcare"  
#         )  
        
#         # Search button moved here  
#         search_button = st.button("🔍 Search", use_container_width=True, type="primary")  
        
#         st.markdown("---")  
        
#         # Subscription settings  
#         st.subheader("📬 Subscription Settings")  
#         frequency = st.selectbox(  
#             "Update Frequency",  
#             ["Daily", "Weekly", "Monthly"],  
#             index=1,  
#             help="How often would you like to receive updates?"  
#         )  
        
#         subscriber_email = st.text_input(  
#             "Email Address",  
#             help="Enter your email to receive automatic updates"  
#         )  
        
#         # Subscribe button moved here  
#         subscribe_button = st.button("📩 Confirm and Subscribe", use_container_width=True)  

#     # Main content area  
#     # Validation helper function  
#     def validate_and_format_journals(selected_journals):  
#         """Validate and format selected journals"""  
#         if not selected_journals:  
#             return []  
        
#         formatted_journals = []  
#         invalid_journals = []  
        
#         for journal in selected_journals:  
#             abbreviated = full_to_abbrev.get(journal)  
#             if abbreviated:  
#                 formatted_journals.append(abbreviated)  
#             else:  
#                 invalid_journals.append(journal)  
        
#         if invalid_journals:  
#             st.error(f"❌ Invalid journals: {', '.join(invalid_journals)}")  
#             return None  
        
#         return formatted_journals  

#     # Enhanced search section with detailed progress tracking  
#     if search_button:  
#         # Rate limiting check first  
#         if not rate_limiter.can_make_request():  
#             st.stop()  

#         try:  
#             # Modified validation - allow search if either journals are selected OR preprints are included  
#             if not selected_journals and not include_preprints:  
#                 st.error("❌ Please select at least one journal or include preprints.")  
#                 st.stop()  

#             # Only process journal validation if journals are selected  
#             if selected_journals:  
#                 formatted_journals = validate_and_format_journals(selected_journals)  
#                 if formatted_journals is None:  # Validation failed  
#                     st.stop()  
#             else:  
#                 formatted_journals = []  

#             # Validation  
#             if not formatted_journals and not include_preprints:  
#                 st.error("❌ Please select at least one journal or enable preprints.")  
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

#             # ========================================  
#             # ENHANCED PROGRESS TRACKING SETUP  
#             # ========================================  
            
#             # Create progress container  
#             progress_container = st.container()  
#             with progress_container:  
#                 st.markdown("<div class='progress-container'>", unsafe_allow_html=True)  
#                 st.markdown("### 🔍 Search Progress")  
                
#                 # Two columns for better layout  
#                 progress_col1, progress_col2 = st.columns([3, 1])  
                
#                 with progress_col1:  
#                     progress_bar = st.progress(0)  
#                     status_text = st.empty()  
                    
#                 with progress_col2:  
#                     progress_stats = st.empty()  
#                     article_counter = st.empty()  
                
#                 st.markdown("</div>", unsafe_allow_html=True)  
            
#             # Enhanced Progress Tracker Class  
#             class DetailedProgressTracker:  
#                 def __init__(self):  
#                     self.total_sources = 0  
#                     self.completed_sources = 0  
#                     self.total_articles_found = 0  
#                     self.processed_articles = 0  
#                     self.current_source = ""  
#                     self.current_source_articles = 0  
#                     self.errors = []  
#                     self.start_time = time.time()  
                    
#                 def set_total_sources(self, total):  
#                     self.total_sources = total  
#                     self.update_display()  
                    
#                 def start_source(self, source_name):  
#                     self.current_source = source_name  
#                     self.current_source_articles = 0  
#                     self.update_display()  
                    
#                 def set_articles_found_for_source(self, count):  
#                     self.current_source_articles = count  
#                     self.total_articles_found += count  
#                     self.update_display()  
                    
#                 def increment_processed(self, count=1):  
#                     self.processed_articles += count  
#                     self.update_display()  
                    
#                 def add_error(self, error_msg):  
#                     self.errors.append(error_msg)  
                    
#                 def complete_source(self):  
#                     self.completed_sources += 1  
#                     self.update_display()  
                    
#                 def update_display(self):  
#                     # Calculate progress percentage  
#                     if self.total_sources > 0:  
#                         source_progress = (self.completed_sources / self.total_sources) * 100  
                        
#                         # Add partial progress for current source  
#                         if self.current_source_articles > 0 and self.processed_articles > 0:  
#                             current_source_progress = min(  
#                                 (self.processed_articles / self.total_articles_found) * 100,   
#                                 100  
#                             )  
#                             # Weight the current source progress  
#                             partial_progress = (current_source_progress / self.total_sources)  
#                             total_progress = min(source_progress + partial_progress, 100)  
#                         else:  
#                             total_progress = source_progress  
                        
#                         # Update progress bar  
#                         progress_bar.progress(int(total_progress))  
                        
#                         # Update status text with detailed information  
#                         elapsed_time = time.time() - self.start_time  
                        
#                         if self.completed_sources == self.total_sources:  
#                             status_text.success(f"✅ Search completed! Found {self.total_articles_found} articles in {elapsed_time:.1f}s")  
#                         else:  
#                             if self.current_source_articles > 0:  
#                                 status_text.info(f"🔍 Processing {self.current_source} | Articles: {self.current_source_articles} found")  
#                             else:  
#                                 status_text.info(f"🔍 Searching {self.current_source}...")  
                        
#                         # Update progress statistics  
#                         progress_stats.metric(  
#                             label="Sources",  
#                             value=f"{self.completed_sources}/{self.total_sources}",  
#                             delta=f"{self.current_source}" if self.current_source else None  
#                         )  
                        
#                         # Update article counter  
#                         article_counter.metric(  
#                             label="Articles Found",  
#                             value=self.total_articles_found,  
#                             delta=f"Processing..." if self.processed_articles < self.total_articles_found else "Complete"  
#                         )  
            
#             # Initialize the enhanced progress tracker  
#             tracker = DetailedProgressTracker()  
            
#             # Set up sources for tracking  
#             sources_to_search = []  
#             if selected_journals:  
#                 sources_to_search.extend(selected_journals)  
#             if include_preprints:  
#                 sources_to_search.extend(["biorxiv", "medrxiv"])  
            
#             tracker.set_total_sources(len(sources_to_search))  
            
#             # Initialize optimized fetcher  
#             optimized_fetcher = OptimizedPubMedFetcher(rate_limiter)  
            
#             # ========================================  
#             # ENHANCED SEARCH EXECUTION  
#             # ========================================  
            
#             all_articles = []  
            
#             # Search PubMed journals with enhanced progress tracking  
#             if selected_journals:  
#                 for journal in selected_journals:  
#                     tracker.start_source(journal)  
                    
#                     try:  
#                         # Build query to get count first - use higher limit for accurate count  
#                         query = build_pubmed_query(journal, start_date, end_date, keywords)  
                        
#                         # Get actual article count using fetch_article_ids_from_pubmed  
#                         from tracking_main import fetch_article_ids_from_pubmed  
#                         pmid_list, actual_count = fetch_article_ids_from_pubmed(query, rate_limiter)  
                        
#                         if pmid_list:  
#                             tracker.set_articles_found_for_source(actual_count)  
                            
#                             # Create progress callback for this journal  
#                             def create_progress_callback():  
#                                 def progress_callback(count):  
#                                     tracker.increment_processed(count)  
#                                 return progress_callback  
                            
#                             # Use the OPTIMIZED fetcher instead of the old one  
#                             articles = optimized_fetcher.fetch_pubmed_articles_optimized(  
#                                 journal, start_date, end_date, keywords,   
#                                 progress_callback=create_progress_callback()  
#                             )  
                            
#                             # Add journal information  
#                             for article in articles:  
#                                 article["Journal"] = journal  
#                                 article["Source"] = "PubMed"  
                            
#                             all_articles.extend(articles)  
                            
#                         else:  
#                             tracker.set_articles_found_for_source(0)  
                            
#                         tracker.complete_source()  
                        
#                     except Exception as e:  
#                         error_msg = f"Error searching {journal}: {str(e)}"  
#                         tracker.add_error(error_msg)  
#                         st.error(f"❌ {error_msg}")  
#                         tracker.complete_source()  
#                         continue  

#             # Search preprints with enhanced progress tracking  
#             # if include_preprints:  
#             #     preprint_servers = ["bioRxiv", "medRxiv"]  
                
#             #     for server in preprint_servers:  
#             #         tracker.start_source(server)  
                    
#             #         try:  
#             #             preprints = fetch_preprints(  
#             #                 server=server,  
#             #                 start_date=start_date,  
#             #                 end_date=end_date,  
#             #                 keywords=raw_keywords  
#             #             )  
                        
#             #             tracker.set_articles_found_for_source(len(preprints))  
                        
#             #             for article in preprints:  
#             #                 article["Journal"] = server  
#             #                 article["Source"] = "Preprint"  
#             #                 tracker.increment_processed()  
                            
#             #             all_articles.extend(preprints)  
#             #             tracker.complete_source()  
                        
#             #         except Exception as e:  
#             #             error_msg = f"Error searching {server}: {str(e)}"  
#             #             tracker.add_error(error_msg)  
#             #             st.error(f"❌ {error_msg}")  
#             #             tracker.complete_source()  
#             #             continue  

#             if include_preprints:  
#                 preprint_servers = ["biorxiv", "medrxiv"]  
                
#                 for server in preprint_servers:  
#                     tracker.start_source(server)  
                    
#                     try:  
#                         # Convert date objects to strings if needed
#                         search_start_date = str(start_date) if hasattr(start_date, 'strftime') else start_date
#                         search_end_date = str(end_date) if hasattr(end_date, 'strftime') else end_date
                        
#                         # ENHANCED: Remove max_results limit to get ALL preprints
#                         preprints = fetch_preprints(  
#                             server=server,  
#                             start_date=search_start_date,  
#                             end_date=search_end_date,  
#                             keywords=raw_keywords,  # Use raw keywords, not formatted
#                             max_results=None  # Remove limit to get all results
#                         )  
                        
#                         tracker.set_articles_found_for_source(len(preprints))  
                        
#                         for article in preprints:  
#                             article["Journal"] = server  
#                             article["Source"] = "Preprint"  
#                             tracker.increment_processed()  
                            
#                         all_articles.extend(preprints)  
#                         tracker.complete_source()  
                        
#                     except Exception as e:  
#                         error_msg = f"Error searching {server}: {str(e)}"  
#                         tracker.add_error(error_msg)  
#                         st.error(f"❌ {error_msg}")  
#                         tracker.complete_source()  
#                         continue

#             # ========================================  
#             # FINAL PROCESSING WITH PROGRESS  
#             # ========================================  
            
#             if all_articles:  
#                 # Show processing status  
#                 status_text.info("📊 Processing and formatting results...")  
                
#                 # Add Source field for PubMed articles  
#                 for article in all_articles:  
#                     if "Source" not in article:  
#                         article["Source"] = "PubMed"  

#                 # Standardize formats  
#                 all_articles = standardize_date_format(all_articles)  
#                 all_articles = standardize_doi_format(all_articles)  
                
#                 # Merge and highlight  
#                 merged = merge_and_highlight_articles(all_articles, [], raw_keywords)  
#                 df = pd.DataFrame(merged)  
                
#                 # Fix DOI format to ensure proper links  
#                 if 'DOI' in df.columns:  
#                     df['DOI'] = df['DOI'].apply(lambda x:   
#                         f"https://doi.org/{x.replace('https://doi.org/', '')}"   
#                         if x and x != "No DOI available" and not x.startswith("https://doi.org/")   
#                         else x  
#                     )  
                
#                 # Final completion  
#                 tracker.update_display()  
                
#                 # Show any errors that occurred  
#                 if tracker.errors:  
#                     with st.expander("⚠️ Errors During Search", expanded=False):  
#                         for error in tracker.errors:  
#                             st.warning(error)  
                
#                 # Display results  
#                 st.markdown("### 📊 Search Results")  
                
#                 # Summary metrics  
#                 col1, col2, col3, col4 = st.columns(4)  
                
#                 with col1:  
#                     st.markdown(f"""  
#                     <div class='stats-card'>  
#                         <h3>📄 {len(df)}</h3>  
#                         <p>Total Articles</p>  
#                     </div>  
#                     """, unsafe_allow_html=True)  
                
#                 with col2:  
#                     unique_journals = df['Journal'].nunique() if 'Journal' in df.columns else 0  
#                     st.markdown(f"""  
#                     <div class='stats-card'>  
#                         <h3>📚 {unique_journals}</h3>  
#                         <p>Journals</p>  
#                     </div>  
#                     """, unsafe_allow_html=True)  
                
#                 with col3:  
#                     pubmed_count = len(df[df['Source'] == 'PubMed']) if 'Source' in df.columns else 0  
#                     st.markdown(f"""  
#                     <div class='stats-card'>  
#                         <h3>🔬 {pubmed_count}</h3>  
#                         <p>PubMed</p>  
#                     </div>  
#                     """, unsafe_allow_html=True)  
                
#                 with col4:  
#                     preprint_count = len(df[df['Source'] == 'Preprint']) if 'Source' in df.columns else 0  
#                     st.markdown(f"""  
#                     <div class='stats-card'>  
#                         <h3>📑 {preprint_count}</h3>  
#                         <p>Preprints</p>  
#                     </div>  
#                     """, unsafe_allow_html=True)  
                
#                 st.markdown("---")  
                
#                 # Display the results table  
#                 st.dataframe(  
#                     df,  
#                     use_container_width=True,  
#                     hide_index=True,  
#                     column_config={  
#                         "Title": st.column_config.TextColumn("Title", width="large"),  
#                         "Abstract": st.column_config.TextColumn("Abstract", width="large"),  
#                         "DOI": st.column_config.LinkColumn("DOI", width="medium"),  
#                         "Publication Date": st.column_config.DateColumn("Publication Date", width="small"),  
#                         "Journal": st.column_config.TextColumn("Journal", width="medium"),  
#                         "Source": st.column_config.TextColumn("Source", width="small")  
#                     }  
#                 )  
                
#                 # Download button  
#                 csv = df.to_csv(index=False)  
#                 st.download_button(  
#                     label="📥 Download Results as CSV",  
#                     data=csv,  
#                     file_name=f"pubmed_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",  
#                     mime="text/csv",  
#                     use_container_width=True  
#                 )  
                
#             else:  
#                 # No results found  
#                 status_text.warning("📭 No articles found matching your criteria.")  
                
#                 # Show suggestions  
#                 st.info("""  
#                 **No results found. Try:**  
#                 - Expanding your date range  
#                 - Using broader keywords  
#                 - Checking different journals  
#                 - Including preprints  
#                 """)  
                
#                 # Show any errors that occurred  
#                 if tracker.errors:  
#                     with st.expander("⚠️ Errors During Search", expanded=False):  
#                         for error in tracker.errors:  
#                             st.warning(error)  

#         except Exception as e:  
#             st.error(f"❌ Unexpected error: {e}")  
#             # Clear progress indicators on error  
#             if 'status_text' in locals():  
#                 status_text.empty()  
#             if 'progress_bar' in locals():  
#                 progress_bar.empty()  

#     # Enhanced subscription confirmation with progress tracking  
#     if subscribe_button:  
#         if not subscriber_email:  
#             st.error("❌ Please provide an email address.")  
#         elif not selected_journals and not include_preprints:  
#             st.error("❌ Please select at least one journal or include preprints.")  
#         else:  
#             # ========================================  
#             # SUBSCRIPTION PROGRESS TRACKING  
#             # ========================================  
            
#             # Create subscription progress container  
#             sub_progress_container = st.container()  
#             with sub_progress_container:  
#                 st.markdown("<div class='progress-container'>", unsafe_allow_html=True)  
#                 st.markdown("### 📬 Subscription Progress")  
#                 sub_col1, sub_col2 = st.columns([3, 1])  
                
#                 with sub_col1:  
#                     sub_progress_bar = st.progress(0)  
#                     sub_status_text = st.empty()  
                    
#                 with sub_col2:  
#                     sub_stats = st.empty()  
                
#                 st.markdown("</div>", unsafe_allow_html=True)  
            
#             try:  
#                 # Step 1: Validate and format journals  
#                 sub_status_text.info("🔍 Validating subscription parameters...")  
#                 sub_progress_bar.progress(10)  
                
#                 formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)] if selected_journals else []  
                
#                 # Step 2: Execute search for subscription  
#                 sub_status_text.info("🔍 Executing search for subscription...")  
#                 sub_progress_bar.progress(20)  
                
#                 # Check if we have existing search results  
#                 if 'df' in locals() and not df.empty:  
#                     # User just performed a search - use those results  
#                     csv_bytes = df.to_csv(index=False).encode("utf-8")  
#                     has_results = True  
#                     result_count = len(df)  
#                     sub_status_text.success(f"✅ Using current search results ({result_count} articles)")  
#                     sub_progress_bar.progress(50)  
#                 else:  
#                     # User is subscribing without searching - execute search now  
#                     sub_status_text.info("🔍 Running fresh search for subscription...")  
#                     search_results = execute_subscription_search(  
#                         journals=formatted_journals,  
#                         keywords=raw_keywords,  
#                         include_preprints=include_preprints,  
#                         frequency=frequency  
#                     )  
                    
#                     if search_results is not None and not search_results.empty:  
#                         csv_bytes = search_results.to_csv(index=False).encode("utf-8")  
#                         has_results = True  
#                         result_count = len(search_results)  
#                         sub_status_text.success(f"✅ Search completed ({result_count} articles found)")  
#                     else:  
#                         csv_bytes = None  
#                         has_results = False  
#                         result_count = 0  
#                         sub_status_text.info("📭 No articles found for current criteria")  
                    
#                     sub_progress_bar.progress(50)  

#                 # Step 3: Store subscription  
#                 sub_status_text.info("💾 Storing subscription in database...")  
#                 sub_progress_bar.progress(60)  
                
#                 result = store_user_subscription(  
#                     email=subscriber_email,  
#                     journals=formatted_journals,  
#                     keywords=raw_keywords,  
#                     frequency=frequency,  
#                     include_preprints=include_preprints  
#                 )  
                
#                 sub_progress_bar.progress(70)  

#                 if result["status"] == "success":  
#                     subscription_data = result["data"]  
#                     subscription_id = subscription_data["id"]  
                    
#                     # Step 4: Generate tokens  
#                     sub_status_text.info("🔐 Generating secure tokens...")  
#                     sub_progress_bar.progress(80)  
                    
#                     unsubscribe_token = generate_unsubscribe_token(subscription_id)  
                    
#                     if unsubscribe_token:  
#                         unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"  
#                     else:  
#                         st.error("❌ Failed to generate unsubscribe token")  
#                         st.stop()  

#                     # Step 5: Prepare email content  
#                     sub_status_text.info("📧 Preparing confirmation email...")  
#                     sub_progress_bar.progress(85)  
                    
#                     # Build source description  
#                     source_list = []
#                     if formatted_journals:
#                         source_list.extend(formatted_journals)
#                     if include_preprints:
#                         source_list.extend(['biorxiv', 'medrxiv'])
#                     source_description = ', '.join(source_list) if source_list else 'No sources selected'

#                     # Step 6: Handle results and email sending
#                     if has_results and csv_bytes is not None:
#                         sub_status_text.info("📥 Generating download link...")
#                         sub_progress_bar.progress(90)
                        
#                         # Generate download token for results
#                         download_token = generate_download_token(
#                             csv_bytes, 
#                             subscriber_email, 
#                             subscription_id=subscription_id
#                         )
                        
#                         if download_token:
#                             download_link = f"{BASE_URL}?token={download_token}&action=download"
                            
#                             # Email content with results
#                             email_body = f"""Hi {subscriber_email},

# You have successfully subscribed to automatic PubMed updates.

# 📊 SUBSCRIPTION DETAILS:
# 📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
# 📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
# 🔍 All Sources: {source_description}
# 🔑 Keywords: {raw_keywords or 'None'}
# 🔁 Frequency: {frequency}

# 📥 YOUR CURRENT RESULTS ({result_count} articles found):
# Your search results are available for download (expires in 3 hours):
# 🔗 {download_link}

# You will receive your next update in {get_next_update_timeframe(frequency)}.

# 🔓 UNSUBSCRIBE:
# {unsubscribe_link}

# – PubMed Tracker Team
#                             """
                            
#                             email_subject = f"📬 Journal Tracker: Subscription Confirmed ({result_count} results)"
#                             success_message = f"✅ Subscription confirmed! {result_count} results sent to {subscriber_email}"
                            
#                         else:
#                             st.error("❌ Failed to generate download link")
#                             st.stop()
#                     else:
#                         # Email content without results
#                         email_body = f"""Hi {subscriber_email},

# You have successfully subscribed to automatic PubMed updates.

# 📊 SUBSCRIPTION DETAILS:
# 📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
# 📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
# 🔍 All Sources: {source_description}
# 🔑 Keywords: {raw_keywords or 'None'}
# 🔁 Frequency: {frequency}

# 📭 CURRENT SEARCH STATUS:
# No articles found matching your criteria for the selected time period.

# You will receive your next update in {get_next_update_timeframe(frequency)} (only if results are found).

# 🔓 UNSUBSCRIBE:
# {unsubscribe_link}

# – PubMed Tracker Team
#                         """
                        
#                         email_subject = "📬 Journal Tracker: Subscription Confirmed (No results)"
#                         success_message = f"✅ Subscription confirmed! Confirmation email sent to {subscriber_email}"
                    
#                     # Step 7: Send email
#                     sub_status_text.info("📨 Sending confirmation email...")
#                     sub_progress_bar.progress(95)
                    
#                     try:
#                         send_email(
#                             to_email=subscriber_email,
#                             subject=email_subject,
#                             body=email_body,
#                             sender_email=BREVO_SENDER_EMAIL,
#                             sender_name=BREVO_SENDER_NAME,
#                             api_key=BREVO_API_KEY
#                         )
                        
#                         # Final completion
#                         sub_progress_bar.progress(100)
#                         sub_status_text.success("✅ Subscription confirmed and email sent!")
                        
#                         # Show final statistics
#                         sub_stats.metric(
#                             label="Subscription Status",
#                             value="Active",
#                             delta=f"{result_count} articles" if has_results else "Ready for updates"
#                         )
                        
#                         st.success(success_message)
                        
#                     except Exception as e:
#                         sub_status_text.error(f"❌ Email sending failed: {e}")
#                         sub_progress_bar.progress(100)
#                         st.error(f"❌ Subscription created but email failed: {e}")
                
#                 elif result["status"] == "error":
#                     sub_status_text.error(f"❌ Subscription failed: {result['message']}")
#                     sub_progress_bar.progress(100)
#                     st.error(f"❌ {result['message']}")
#                     if "Maximum of 3" in result["message"]:
#                         st.info("💡 **Tip:** You can manage your existing subscriptions by using the unsubscribe link in any of your previous emails.")
#                 else:
#                     sub_status_text.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")
#                     sub_progress_bar.progress(100)
#                     st.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")
            
#             except Exception as e:
#                 sub_status_text.error(f"❌ Unexpected error: {e}")
#                 sub_progress_bar.progress(100)
#                 st.error(f"❌ Subscription failed: {e}")

#     # Footer
#     st.markdown("---")
#     st.markdown(
#         """
#         <div style='text-align: center; color: #666; font-size: 0.8em; padding: 20px 0;'>
#             Please report issues to <a href='https://github.com/rzhan186/Journal_tracker/issues' target='_blank'>GitHub</a>
#         </div>
#         """, 
#         unsafe_allow_html=True
# )

import streamlit as st
import pandas as pd
from datetime import datetime, timedelta
import time
import os
from tracking_main import (
    fetch_pubmed_articles_by_date,
    format_journal_abbreviation,
    fetch_preprints,
    load_pubmed_journal_abbreviations,
    build_pubmed_query,
    merge_and_highlight_articles,
    standardize_date_format,
    standardize_doi_format,
    compile_keyword_filter,
    format_boolean_keywords_for_pubmed,
    execute_subscription_search,
    OptimizedPubMedFetcher
)
from RateLimit import PubMedRateLimit
from store_subscription import store_user_subscription, generate_unsubscribe_token
from email_dispatcher import send_email, generate_download_token, get_next_update_timeframe

# --- CONFIGURATION ---
BASE_URL = os.getenv("BASE_URL", "https://journaltracker.streamlit.app/")
BREVO_API_KEY = os.getenv("BREVO_API_KEY")
BREVO_SENDER_EMAIL = os.getenv("BREVO_SENDER_EMAIL")
BREVO_SENDER_NAME = os.getenv("BREVO_SENDER_NAME", "PubMed Journal Tracker")

rate_limiter = PubMedRateLimit()
if not os.getenv("NCBI_API_KEY"):
    st.warning("⚠️ NCBI API key not configured. Performance may be limited.")

st.set_page_config(
    page_title="PubMed Journal Tracker",
    page_icon="📚",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .main-header {text-align: center; color: #2E86AB; margin-bottom: 2rem;}
    .progress-container {background-color: #f0f2f6; padding: 1rem; border-radius: 8px; margin: 0rem 0 1rem 0;}
    .stats-card {background-color: #ffffff; padding: 1.5rem; border-radius: 10px; border: 1px solid #e0e0e0; margin-bottom: 1rem; text-align: center;}
    .info-card {background-color: #e8f4f8; padding: 1rem; border-radius: 8px; margin-bottom: 1rem;}
    .stProgress .st-bo {background-color: #2E86AB;}
</style>""", unsafe_allow_html=True)

# -- Helper to always format dates as strings
def ensure_ymd(date_obj):
    try:
        return date_obj.strftime("%Y-%m-%d")
    except AttributeError:
        return str(date_obj)

# Unsubscribe / download via token in URL
if 'token' in st.query_params:
    token = st.query_params['token']
    action = st.query_params.get('action', 'unsubscribe')

    if action == 'download':
        from app_csv_downloader import handle_download
        handle_download(token)
    else:
        from app_unsubscribe import handle_unsubscribe
        handle_unsubscribe(token)

else:
    st.markdown("<h1 class='main-header'>📚 PubMed Journal Tracker</h1>", unsafe_allow_html=True)
    
    with st.sidebar:
        st.header("🔍 Search Parameters")
        st.subheader("📅 Date Range")
        col1, col2 = st.columns(2)
        with col1:
            start_date = st.date_input("Start Date", value=datetime.now() - timedelta(days=30))
        with col2:
            end_date = st.date_input("End Date", value=datetime.now())
        if start_date >= end_date:
            st.error("❌ Start date must be before end date")
        
        # Journal selection
        st.subheader("📘 Journal Selection")
        journal_dict = load_pubmed_journal_abbreviations()
        full_to_abbrev = {v: k for k, v in journal_dict.items()}
        # Optionally add some common missing journals
        for journal in ["The Lancet", "BMJ"]:
            if journal not in full_to_abbrev:
                full_to_abbrev[journal] = journal
        journal_options = sorted(full_to_abbrev.keys())
        selected_journals = st.multiselect(
            "Select Journals", options=journal_options
        )

        # Preprint checkbox
        include_preprints = st.checkbox(
            "📑 Search Preprints (bioRxiv, medRxiv)",
            value=False
        )

        # Keywords
        st.subheader("🔑 Keywords")
        raw_keywords = st.text_area(
            "Enter Keywords (optional)",
            placeholder="e.g., machine learning, AI AND healthcare"
        )

        search_button = st.button("🔍 Search", use_container_width=True, type="primary")
        st.markdown("---")
        st.subheader("📬 Subscription Settings")
        frequency = st.selectbox("Update Frequency", ["Daily", "Weekly", "Monthly"], index=1)
        subscriber_email = st.text_input("Email Address")
        subscribe_button = st.button("📩 Confirm and Subscribe", use_container_width=True)

    def validate_and_format_journals(selected_journals):
        if not selected_journals:
            return []
        formatted_journals = []
        for journal in selected_journals:
            abbreviated = full_to_abbrev.get(journal)
            if abbreviated:
                formatted_journals.append(abbreviated)
            else:
                st.error(f"❌ Invalid journal: {journal}")
                return None
        return formatted_journals

    # ---- SEARCH LOGIC ----
    if search_button:
        if not rate_limiter.can_make_request():
            st.stop()
        try:
            if not selected_journals and not include_preprints:
                st.error("❌ Please select at least one journal or include preprints.")
                st.stop()
            if selected_journals:
                formatted_journals = validate_and_format_journals(selected_journals)
                if formatted_journals is None:
                    st.stop()
            else:
                formatted_journals = []
            if not formatted_journals and not include_preprints:
                st.error("❌ Please select at least one journal or enable preprints.")
                st.stop()
            # Keywords
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

            # -- PROGRESS BAR SETUP --
            progress_container = st.container()
            with progress_container:
                st.markdown("<div class='progress-container'>", unsafe_allow_html=True)
                st.markdown("### 🔍 Search Progress")
                progress_col1, progress_col2 = st.columns([3, 1])
                with progress_col1:
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                with progress_col2:
                    progress_stats = st.empty()
                    article_counter = st.empty()
                st.markdown("</div>", unsafe_allow_html=True)

            # Progress Tracker
            class DetailedProgressTracker:
                def __init__(self):
                    self.total_sources = 0
                    self.completed_sources = 0
                    self.total_articles_found = 0
                    self.processed_articles = 0
                    self.current_source = ""
                    self.current_source_articles = 0
                    self.errors = []
                    self.start_time = time.time()
                def set_total_sources(self, total):
                    self.total_sources = total
                    self.update_display()
                def start_source(self, source_name):
                    self.current_source = source_name
                    self.current_source_articles = 0
                    self.update_display()
                def set_articles_found_for_source(self, count):
                    self.current_source_articles = count
                    self.total_articles_found += count
                    self.update_display()
                def increment_processed(self, count=1):
                    self.processed_articles += count
                    self.update_display()
                def add_error(self, error_msg):
                    self.errors.append(error_msg)
                def complete_source(self):
                    self.completed_sources += 1
                    self.update_display()
                def update_display(self):
                    if self.total_sources > 0:
                        source_progress = (self.completed_sources / self.total_sources) * 100
                        if self.current_source_articles > 0 and self.processed_articles > 0:
                            current_source_progress = min(
                                (self.processed_articles / self.total_articles_found) * 100, 100
                            )
                            partial_progress = (current_source_progress / self.total_sources)
                            total_progress = min(source_progress + partial_progress, 100)
                        else:
                            total_progress = source_progress
                        progress_bar.progress(int(total_progress))
                        elapsed_time = time.time() - self.start_time
                        if self.completed_sources == self.total_sources:
                            status_text.success(
                                f"✅ Search completed! Found {self.total_articles_found} articles in {elapsed_time:.1f}s"
                            )
                        else:
                            if self.current_source_articles > 0:
                                status_text.info(f"🔍 Processing {self.current_source} | Articles: {self.current_source_articles} found")
                            else:
                                status_text.info(f"🔍 Searching {self.current_source}...")
                        progress_stats.metric(
                            label="Sources",
                            value=f"{self.completed_sources}/{self.total_sources}",
                            delta=f"{self.current_source}" if self.current_source else None
                        )
                        article_counter.metric(
                            label="Articles Found",
                            value=self.total_articles_found,
                            delta=f"Processing..." if self.processed_articles < self.total_articles_found else "Complete"
                        )

            tracker = DetailedProgressTracker()
            
            # -- List sources --
            sources_to_search = []
            if selected_journals:
                sources_to_search.extend(selected_journals)
            if include_preprints:
                sources_to_search.extend(["biorxiv", "medrxiv"])
            tracker.set_total_sources(len(sources_to_search))

            optimized_fetcher = OptimizedPubMedFetcher(rate_limiter)
            all_articles = []

            # === PUBMED SEARCH ===
            if selected_journals:
                for journal in selected_journals:
                    tracker.start_source(journal)
                    try:
                        query = build_pubmed_query(journal, start_date, end_date, keywords)
                        from tracking_main import fetch_article_ids_from_pubmed
                        pmid_list, actual_count = fetch_article_ids_from_pubmed(query, rate_limiter)
                        if pmid_list:
                            tracker.set_articles_found_for_source(actual_count)
                            def create_progress_callback():
                                def progress_callback(count):
                                    tracker.increment_processed(count)
                                return progress_callback
                            articles = optimized_fetcher.fetch_pubmed_articles_optimized(
                                journal, start_date, end_date, keywords, progress_callback=create_progress_callback()
                            )
                            for article in articles:
                                article["Journal"] = journal
                                article["Source"] = "PubMed"
                            all_articles.extend(articles)
                        else:
                            tracker.set_articles_found_for_source(0)
                        tracker.complete_source()
                    except Exception as e:
                        error_msg = f"Error searching {journal}: {str(e)}"
                        tracker.add_error(error_msg)
                        st.error(f"❌ {error_msg}")
                        tracker.complete_source()
                        continue

            # === PREPRINT SEARCH ===
            if include_preprints:
                preprint_servers = ["biorxiv", "medrxiv"]
                search_start_date = ensure_ymd(start_date)
                search_end_date = ensure_ymd(end_date)
                for server in preprint_servers:
                    tracker.start_source(server)
                    try:
                        preprints = fetch_preprints(
                            server=server,
                            start_date=search_start_date,
                            end_date=search_end_date,
                            keywords=raw_keywords if raw_keywords.strip() else None,
                            max_results=None
                        )
                        tracker.set_articles_found_for_source(len(preprints))
                        for article in preprints:
                            article["Journal"] = server
                            article["Source"] = "Preprint"
                            tracker.increment_processed()
                        all_articles.extend(preprints)
                        tracker.complete_source()
                    except Exception as e:
                        error_msg = f"Error searching {server}: {str(e)}"
                        tracker.add_error(error_msg)
                        st.error(f"❌ {error_msg}")
                        tracker.complete_source()
                        continue

            # === FINAL PROCESSING & RESULT DISPLAY ===
            if all_articles:
                status_text.info("📊 Processing and formatting results...")
                for article in all_articles:
                    if "Source" not in article:
                        article["Source"] = "PubMed"
                all_articles = standardize_date_format(all_articles)
                all_articles = standardize_doi_format(all_articles)
                merged = merge_and_highlight_articles(all_articles, [], raw_keywords)
                df = pd.DataFrame(merged)
                if 'DOI' in df.columns:
                    df['DOI'] = df['DOI'].apply(lambda x:
                        f"https://doi.org/{x.replace('https://doi.org/', '')}"
                        if x and x != "No DOI available" and not x.startswith("https://doi.org/")
                        else x
                    )
                tracker.update_display()
                if tracker.errors:
                    with st.expander("⚠️ Errors During Search", expanded=False):
                        for error in tracker.errors:
                            st.warning(error)
                st.markdown("### 📊 Search Results")
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.markdown(
                        f"<div class='stats-card'><h3>📄 {len(df)}</h3><p>Total Articles</p></div>", unsafe_allow_html=True)
                with col2:
                    unique_journals = df['Journal'].nunique() if 'Journal' in df.columns else 0
                    st.markdown(
                        f"<div class='stats-card'><h3>📚 {unique_journals}</h3><p>Journals</p></div>", unsafe_allow_html=True)
                with col3:
                    pubmed_count = len(df[df['Source'] == 'PubMed']) if 'Source' in df.columns else 0
                    st.markdown(
                        f"<div class='stats-card'><h3>🔬 {pubmed_count}</h3><p>PubMed</p></div>", unsafe_allow_html=True)
                with col4:
                    preprint_count = len(df[df['Source'] == 'Preprint']) if 'Source' in df.columns else 0
                    st.markdown(
                        f"<div class='stats-card'><h3>📑 {preprint_count}</h3><p>Preprints</p></div>", unsafe_allow_html=True)
                st.markdown("---")
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True
                )
                csv = df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Results as CSV",
                    data=csv,
                    file_name=f"pubmed_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            else:
                status_text.warning("📭 No articles found matching your criteria.")
                st.info("""**No results found. Try:**
- Expanding your date range
- Using broader keywords
- Checking different journals
- Including preprints
""")
                if tracker.errors:
                    with st.expander("⚠️ Errors During Search", expanded=False):
                        for error in tracker.errors:
                            st.warning(error)
        except Exception as e:
            st.error(f"❌ Unexpected error: {e}")
            if 'status_text' in locals():
                status_text.empty()
            if 'progress_bar' in locals():
                progress_bar.empty()

    # --- SUBSCRIPTION CONFIRMATION ---
    if subscribe_button:
        if not subscriber_email:
            st.error("❌ Please provide an email address.")
        elif not selected_journals and not include_preprints:
            st.error("❌ Please select at least one journal or include preprints.")
        else:
            sub_progress_container = st.container()
            with sub_progress_container:
                st.markdown("<div class='progress-container'>", unsafe_allow_html=True)
                st.markdown("### 📬 Subscription Progress")
                sub_col1, sub_col2 = st.columns([3, 1])
                with sub_col1:
                    sub_progress_bar = st.progress(0)
                    sub_status_text = st.empty()
                with sub_col2:
                    sub_stats = st.empty()
                st.markdown("</div>", unsafe_allow_html=True)
            try:
                sub_status_text.info("🔍 Validating subscription parameters...")
                sub_progress_bar.progress(10)
                formatted_journals = [full_to_abbrev.get(name) for name in selected_journals if full_to_abbrev.get(name)] if selected_journals else []
                sub_status_text.info("🔍 Executing search for subscription...")
                sub_progress_bar.progress(20)
                if 'df' in locals() and not df.empty:
                    csv_bytes = df.to_csv(index=False).encode("utf-8")
                    has_results = True
                    result_count = len(df)
                    sub_status_text.success(f"✅ Using current search results ({result_count} articles)")
                    sub_progress_bar.progress(50)
                else:
                    sub_status_text.info("🔍 Running fresh search for subscription...")
                    search_results = execute_subscription_search(
                        journals=formatted_journals,
                        keywords=raw_keywords,
                        include_preprints=include_preprints,
                        frequency=frequency
                    )
                    if search_results is not None and not search_results.empty:
                        csv_bytes = search_results.to_csv(index=False).encode("utf-8")
                        has_results = True
                        result_count = len(search_results)
                        sub_status_text.success(f"✅ Search completed ({result_count} articles found)")
                    else:
                        csv_bytes = None
                        has_results = False
                        result_count = 0
                        sub_status_text.info("📭 No articles found for current criteria")
                    sub_progress_bar.progress(50)
                sub_status_text.info("💾 Storing subscription in database...")
                sub_progress_bar.progress(60)
                result = store_user_subscription(
                    email=subscriber_email,
                    journals=formatted_journals,
                    keywords=raw_keywords,
                    frequency=frequency,
                    include_preprints=include_preprints
                )
                sub_progress_bar.progress(70)
                if result["status"] == "success":
                    subscription_data = result["data"]
                    subscription_id = subscription_data["id"]
                    sub_status_text.info("🔐 Generating secure tokens...")
                    sub_progress_bar.progress(80)
                    unsubscribe_token = generate_unsubscribe_token(subscription_id)
                    if unsubscribe_token:
                        unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"
                    else:
                        st.error("❌ Failed to generate unsubscribe token")
                        st.stop()
                    source_list = []
                    if formatted_journals:
                        source_list.extend(formatted_journals)
                    if include_preprints:
                        source_list.extend(['biorxiv', 'medrxiv'])
                    source_description = ', '.join(source_list) if source_list else 'No sources selected'
                    if has_results and csv_bytes is not None:
                        sub_status_text.info("📥 Generating download link...")
                        sub_progress_bar.progress(90)
                        download_token = generate_download_token(
                            csv_bytes,
                            subscriber_email,
                            subscription_id=subscription_id
                        )
                        if download_token:
                            download_link = f"{BASE_URL}?token={download_token}&action=download"
                            email_body = f"""Hi {subscriber_email},

You have successfully subscribed to automatic PubMed updates.

📊 SUBSCRIPTION DETAILS:
📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
🔍 All Sources: {source_description}
🔑 Keywords: {raw_keywords or 'None'}
🔁 Frequency: {frequency}

📥 YOUR CURRENT RESULTS ({result_count} articles found):
🔗 {download_link}

You will receive your next update in {get_next_update_timeframe(frequency)}.

🔓 UNSUBSCRIBE:
{unsubscribe_link}

– PubMed Tracker Team
                            """
                            email_subject = f"📬 Journal Tracker: Subscription Confirmed ({result_count} results)"
                            success_message = f"✅ Subscription confirmed! {result_count} results sent to {subscriber_email}"
                        else:
                            st.error("❌ Failed to generate download link")
                            st.stop()
                    else:
                        email_body = f"""Hi {subscriber_email},

You have successfully subscribed to automatic PubMed updates.

📊 SUBSCRIPTION DETAILS:
📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
🔍 All Sources: {source_description}
🔑 Keywords: {raw_keywords or 'None'}
🔁 Frequency: {frequency}

📭 CURRENT SEARCH STATUS:
No articles found matching your criteria for the selected time period.

You will receive your next update in {get_next_update_timeframe(frequency)} (only if results are found).

🔓 UNSUBSCRIBE:
{unsubscribe_link}

– PubMed Tracker Team
                        """
                        email_subject = "📬 Journal Tracker: Subscription Confirmed (No results)"
                        success_message = f"✅ Subscription confirmed! Confirmation email sent to {subscriber_email}"
                    sub_status_text.info("📨 Sending confirmation email...")
                    sub_progress_bar.progress(95)
                    try:
                        send_email(
                            to_email=subscriber_email,
                            subject=email_subject,
                            body=email_body,
                            sender_email=BREVO_SENDER_EMAIL,
                            sender_name=BREVO_SENDER_NAME,
                            api_key=BREVO_API_KEY
                        )
                        sub_progress_bar.progress(100)
                        sub_status_text.success("✅ Subscription confirmed and email sent!")
                        sub_stats.metric(
                            label="Subscription Status",
                            value="Active",
                            delta=f"{result_count} articles" if has_results else "Ready for updates")
                        st.success(success_message)
                    except Exception as e:
                        sub_status_text.error(f"❌ Email sending failed: {e}")
                        sub_progress_bar.progress(100)
                        st.error(f"❌ Subscription created but email failed: {e}")
                elif result["status"] == "error":
                    sub_status_text.error(f"❌ Subscription failed: {result['message']}")
                    sub_progress_bar.progress(100)
                    st.error(f"❌ {result['message']}")
                    if "Maximum of 3" in result["message"]:
                        st.info("💡 **Tip:** You can manage your existing subscriptions with the unsubscribe link in your previous emails.")
                else:
                    sub_status_text.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")
                    sub_progress_bar.progress(100)
                    st.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")
            except Exception as e:
                sub_status_text.error(f"❌ Unexpected error: {e}")
                sub_progress_bar.progress(100)
                st.error(f"❌ Subscription failed: {e}")

    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center; color: #666; font-size: 0.8em; padding: 20px 0;'>
            Please report issues to <a href='https://github.com/rzhan186/Journal_tracker/issues' target='_blank'>GitHub</a>
        </div>
        """,
        unsafe_allow_html=True
    )