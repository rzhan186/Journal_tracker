import streamlit as st
import pandas as pd
from datetime import datetime, timedelta
import time
import re
from tracking_main import (
    fetch_pubmed_articles_by_date,
    fetch_article_ids_from_pubmed,
    format_journal_abbreviation,
    fetch_preprints,
    fetch_preprints_with_progress,
    load_pubmed_journal_abbreviations,
    build_pubmed_query,
    merge_and_highlight_articles,
    standardize_date_format,
    standardize_doi_format,
    compile_keyword_filter,
    format_boolean_keywords_for_pubmed,
    execute_subscription_search,
    generate_bibtex_from_dataframe,
    OptimizedPubMedFetcher
)
from RateLimit import PubMedRateLimit
from store_subscription import store_user_subscription, generate_unsubscribe_token
from email_dispatcher import send_email, generate_download_token, get_next_update_timeframe
import os

# Configuration
BASE_URL = os.getenv("BASE_URL", "https://journaltracker.streamlit.app/")
BREVO_API_KEY = os.getenv("BREVO_API_KEY")
BREVO_SENDER_EMAIL = os.getenv("BREVO_SENDER_EMAIL")
BREVO_SENDER_NAME = os.getenv("BREVO_SENDER_NAME", "PubMed Journal Tracker")

# Initialize rate limiter
rate_limiter = PubMedRateLimit()

if not os.getenv("NCBI_API_KEY"):
    st.warning("⚠️ NCBI API key not configured. Performance may be limited.")


def validate_and_format_journals(selected_journals, full_to_abbrev):
    """Validate and format selected journals, returning abbreviated names."""
    if not selected_journals:
        return []

    formatted_journals = []
    invalid_journals = []

    for journal in selected_journals:
        abbreviated = full_to_abbrev.get(journal)
        if abbreviated:
            formatted_journals.append(abbreviated)
        else:
            invalid_journals.append(journal)

    if invalid_journals:
        st.error(f"❌ Invalid journals: {', '.join(invalid_journals)}")
        return None

    return formatted_journals


class DetailedProgressTracker:
    def __init__(self, progress_bar, status_text, progress_stats, article_counter):
        self.progress_bar = progress_bar
        self.status_text = status_text
        self.progress_stats = progress_stats
        self.article_counter = article_counter
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
        if self.current_source_articles == 0:
            self.current_source_articles = count
            self.total_articles_found += count
        else:
            self.total_articles_found = self.total_articles_found - self.current_source_articles + count
            self.current_source_articles = count
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
        if self.total_sources == 0:
            return

        source_progress = (self.completed_sources / self.total_sources) * 100

        if self.current_source_articles > 0 and self.processed_articles > 0:
            current_source_progress = min(
                (self.processed_articles / self.total_articles_found) * 100, 100
            )
            total_progress = min(source_progress + (current_source_progress / self.total_sources), 100)
        else:
            total_progress = source_progress

        self.progress_bar.progress(int(total_progress))

        elapsed_time = time.time() - self.start_time

        if self.completed_sources == self.total_sources:
            self.status_text.success(f"✅ Search completed! Found {self.total_articles_found} articles in {elapsed_time:.1f}s")
        elif self.current_source_articles > 0:
            self.status_text.info(f"🔍 Processing {self.current_source} | Articles: {self.current_source_articles} found")
        else:
            self.status_text.info(f"🔍 Searching {self.current_source}...")

        self.progress_stats.metric(
            label="Sources",
            value=f"{self.completed_sources}/{self.total_sources}",
            delta=f"{self.current_source}" if self.current_source else None
        )

        if self.completed_sources == self.total_sources:
            self.article_counter.metric(label="Articles Found", value=self.total_articles_found, delta="Complete")
        else:
            self.article_counter.metric(
                label="Articles Found",
                value=self.total_articles_found,
                delta="Processing..." if self.current_source else "Complete"
            )

# Page configuration
st.set_page_config(
    page_title="PubMed Journal Tracker",
    page_icon="📚",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state variables
if 'search_results' not in st.session_state:
    st.session_state.search_results = None
if 'search_completed' not in st.session_state:
    st.session_state.search_completed = False
if 'search_timestamp' not in st.session_state:
    st.session_state.search_timestamp = None


# Custom CSS for better styling  
st.markdown("""  
<style>  
    .main-header {  
        text-align: center;  
        color: #2E86AB;  
        margin-bottom: 2rem;  
    }  
    .search-container {  
        background-color: #f8f9fa;  
        padding: 1.5rem;  
        border-radius: 10px;  
        margin-bottom: 2rem;  
        margin-top: 0rem;  
    }  
    .result-container {  
        background-color: #ffffff;  
        padding: 1rem;  
        border-radius: 8px;  
        border: 1px solid #e0e0e0;  
    }  
    .progress-container {  
        background-color: #f0f2f6;  
        padding: 1rem;  
        border-radius: 8px;  
        margin: 0rem 0 1rem 0;  
    }  
    .stats-card {  
        background-color: #ffffff;  
        padding: 1.5rem;  
        border-radius: 10px;  
        border: 1px solid #e0e0e0;  
        margin-bottom: 1rem;  
        text-align: center;  
    }  
    .feature-card {  
        background-color: #f8f9fa;  
        padding: 1.5rem;  
        border-radius: 10px;  
        border-left: 4px solid #2E86AB;  
        margin-bottom: 1rem;  
    }  
    .info-card {  
        background-color: #e8f4f8;  
        padding: 1rem;  
        border-radius: 8px;  
        margin-bottom: 1rem;  
    }  
    .stProgress .st-bo {  
        background-color: #2E86AB;  
    }  
</style>  
""", unsafe_allow_html=True)  

# Check for unsubscribe token in the URL  
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
    # Main application interface  
    st.markdown("<h1 class='main-header'>📚 PubMed Journal Tracker</h1>", unsafe_allow_html=True)  
    
    # Sidebar for search parameters  
    with st.sidebar:  
        st.header("🔍 Search Parameters")  
        
        # Date range selection  
        st.subheader("📅 Date Range")  
        col1, col2 = st.columns(2)  
        with col1:  
            start_date = st.date_input(  
                "Start Date",  
                value=datetime.now() - timedelta(days=30),  
                help="Select the start date for your search"  
            )  
        with col2:  
            end_date = st.date_input(  
                "End Date",  
                value=datetime.now(),  
                help="Select the end date for your search"  
            )  
        
        # Validate date range  
        if start_date >= end_date:  
            st.error("❌ Start date must be before end date")  
        
        # Journal selection  
        st.subheader("📘 Journal Selection")  
        
        # Load journal abbreviations  
        journal_dict = load_pubmed_journal_abbreviations()  
        full_to_abbrev = {v: k for k, v in journal_dict.items()}  
        
        # Add missing journals  
        missing_journals = ["The Lancet", "BMJ"]  
        for journal in missing_journals:  
            if journal not in full_to_abbrev:  
                full_to_abbrev[journal] = journal  
        
        journal_options = sorted(full_to_abbrev.keys())  
        
        # Multi-select for journals  
        selected_journals = st.multiselect(  
            "Select Journals",  
            options=journal_options,  
            help="Choose one or more journals to search"  
        )  
        
        # Preprint option  
        include_preprints = st.checkbox(  
            "📑 Search Preprints (bioRxiv, medRxiv)",  
            value=False,  
            help="Include preprint servers in your search"  
        )  
        
        # Keywords input  
        st.subheader("🔑 Keywords")  
        raw_keywords = st.text_area(  
            "Enter Keywords (optional)",  
            help="Enter keywords or boolean expressions (e.g., 'cancer AND treatment', 'diabetes OR obesity')",  
            placeholder="e.g., machine learning, AI AND healthcare"  
        )  
        
        # Search button moved here  
        search_button = st.button("🔍 Search", use_container_width=True, type="primary")  
        
        st.markdown("---")  
        
        # Subscription settings  
        st.subheader("📬 Subscription Settings")  
        frequency = st.selectbox(  
            "Update Frequency",  
            ["Daily", "Weekly", "Monthly"],  
            index=1,  
            help="How often would you like to receive updates?"  
        )  
        
        subscriber_email = st.text_input(  
            "Email Address",  
            help="Enter your email to receive automatic updates"  
        )  
        
        # Subscribe button moved here  
        subscribe_button = st.button("📩 Confirm and Subscribe", use_container_width=True)  

    # Search section

    if search_button:  
        # Rate limiting check first  
        if not rate_limiter.can_make_request():  
            st.stop()  

        try:  
            # Modified validation - allow search if either journals are selected OR preprints are included  
            if not selected_journals and not include_preprints:  
                st.error("❌ Please select at least one journal or include preprints.")  
                st.stop()  

            # Only process journal validation if journals are selected
            if selected_journals:
                formatted_journals = validate_and_format_journals(selected_journals, full_to_abbrev)
                if formatted_journals is None:  # Validation failed
                    st.stop()  
            else:  
                formatted_journals = []  

            # Validation  
            if not formatted_journals and not include_preprints:  
                st.error("❌ Please select at least one journal or enable preprints.")  
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

            # ========================================  
            # ENHANCED PROGRESS TRACKING SETUP  
            # ========================================  
            
            # Create progress container  
            progress_container = st.container()  
            with progress_container:  
                st.markdown("---")  
                st.markdown("### 🔍 Search Progress")  
                
                # Two columns for better layout  
                progress_col1, progress_col2 = st.columns([3, 1])  
                
                with progress_col1:  
                    progress_bar = st.progress(0)  
                    status_text = st.empty()  
                    
                with progress_col2:  
                    progress_stats = st.empty()  
                    article_counter = st.empty()  
                
                st.markdown("</div>", unsafe_allow_html=True)  
            
            tracker = DetailedProgressTracker(progress_bar, status_text, progress_stats, article_counter)

            # Set up sources for tracking  
            sources_to_search = []  
            if selected_journals:  
                sources_to_search.extend(selected_journals)  
            if include_preprints:  
                sources_to_search.extend(["biorxiv", "medrxiv"])  
            
            tracker.set_total_sources(len(sources_to_search))  
            
            # Initialize optimized fetcher  
            optimized_fetcher = OptimizedPubMedFetcher(rate_limiter)  
            
            
            # ========================================  
            # ENHANCED SEARCH EXECUTION  
            # ========================================  
            
            all_articles = []  
            
            # Search PubMed journals with enhanced progress tracking  
            if selected_journals:
                for journal in selected_journals:
                    tracker.start_source(journal)

                    try:
                        query = build_pubmed_query(journal, start_date, end_date, keywords)
                        pmid_list, actual_count = fetch_article_ids_from_pubmed(query, rate_limiter)  
                        
                        if pmid_list:  
                            tracker.set_articles_found_for_source(actual_count)  
                            
                            # Create progress callback for this journal  
                            def create_progress_callback():  
                                def progress_callback(count):  
                                    tracker.increment_processed(count)  
                                return progress_callback  
                            
                            # Use the OPTIMIZED fetcher instead of the old one  
                            articles = optimized_fetcher.fetch_pubmed_articles_optimized(  
                                journal, start_date, end_date, keywords,   
                                progress_callback=create_progress_callback()  
                            )  
                            
                            # Add journal information  
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

            # In the preprint search loop - add this check
            if include_preprints:
                for server in ["biorxiv", "medrxiv"]:
                    tracker.start_source(server)
                    
                    try:
                        # Convert date objects to strings if needed
                        search_start_date = str(start_date) if hasattr(start_date, 'strftime') else start_date
                        search_end_date = str(end_date) if hasattr(end_date, 'strftime') else end_date
                        
                        # Create a simple progress callback for preprints
                        def preprint_progress_callback(found_count):
                            # Update the tracker with found articles
                            tracker.set_articles_found_for_source(found_count)
                            
                        # Use raw_keywords directly, not processed keywords
                        preprints = fetch_preprints_with_progress(  
                            server=server,  
                            start_date=search_start_date,  
                            end_date=search_end_date,  
                            keywords=raw_keywords.strip() if raw_keywords else None,
                            max_results=None,
                            progress_callback=preprint_progress_callback
                        )
                        
                        # Update tracker with final count
                        tracker.set_articles_found_for_source(len(preprints))  
                        
                        for article in preprints:  
                            article["Journal"] = server  
                            article["Source"] = "Preprint"
                            
                        all_articles.extend(preprints)  
                        tracker.complete_source()  
                        
                    except Exception as e:  
                        error_msg = f"Error searching {server}: {str(e)}"  
                        tracker.add_error(error_msg)  
                        st.error(f"❌ {error_msg}")  
                        tracker.complete_source()  
                        continue

            # ========================================  
            # FINAL PROCESSING WITH PROGRESS  
            # ========================================  
            
            if all_articles:  
                # Show processing status  
                status_text.info("📊 Processing and formatting results...")  
                
                # Add Source field for PubMed articles  
                for article in all_articles:  
                    if "Source" not in article:  
                        article["Source"] = "PubMed"  

                # Standardize formats  
                all_articles = standardize_date_format(all_articles)  
                all_articles = standardize_doi_format(all_articles)  
                
                # Merge and highlight  
                merged = merge_and_highlight_articles(all_articles, [], raw_keywords)  
                df = pd.DataFrame(merged)  
                
                # Store results in session state to prevent loss on download
                st.session_state.search_results = df
                st.session_state.search_completed = True
                st.session_state.search_timestamp = datetime.now()

                # Fix DOI format to ensure proper links  
                if 'DOI' in df.columns:  
                    df['DOI'] = df['DOI'].apply(lambda x:   
                        f"https://doi.org/{x.replace('https://doi.org/', '')}"   
                        if x and x != "No DOI available" and not x.startswith("https://doi.org/")   
                        else x  
                    )  
                
                # Final completion  
                tracker.update_display()  
                
                # Show any errors that occurred  
                if tracker.errors:  
                    with st.expander("⚠️ Errors During Search", expanded=False):  
                        for error in tracker.errors:  
                            st.warning(error)  
                
                # Display results  
                status_text.success(f"✅ Search completed! Found {len(df)} articles")

            else:  
                # No results found  
                status_text.warning("📭 No articles found matching your criteria.")  
                
                # Show suggestions  
                st.info("""  
                **No results found. Try:**  
                - Expanding your date range  
                - Using broader keywords  
                - Checking different journals  
                - Including preprints  
                """)  
                
                # Show any errors that occurred  
                if tracker.errors:  
                    with st.expander("⚠️ Errors During Search", expanded=False):  
                        for error in tracker.errors:  
                            st.warning(error)  

        except Exception as e:  
            st.error(f"❌ Unexpected error: {e}")  
            # Clear progress indicators on error  
            if 'status_text' in locals():  
                status_text.empty()  
            if 'progress_bar' in locals():  
                progress_bar.empty()  


    # ========================================
    # PERSISTENT RESULTS DISPLAY SECTION  
    # ========================================

    # Display results if they exist in session state (regardless of whether search was just run)
    if st.session_state.get('search_results') is not None and not st.session_state.search_results.empty:
        df = st.session_state.search_results
        
        st.markdown("---")  
        st.markdown("### 📊 Search Results")  
        
        # Show when the search was performed
        if hasattr(st.session_state, 'search_timestamp'):
            st.caption(f"Search performed: {st.session_state.search_timestamp.strftime('%Y-%m-%d %H:%M:%S')}")

        # Download buttons in two columns - ALWAYS AVAILABLE when results exist
        col1, col2 = st.columns(2)

        with col1:
            csv = df.to_csv(index=False)  
            st.download_button(  
                label="📥 Download CSV",  
                data=csv,  
                file_name=f"pubmed_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",  
                mime="text/csv",  
                use_container_width=True,
                key="persistent_csv_download"
            )

        with col2:
            bibtex_content = generate_bibtex_from_dataframe(df)
            st.download_button(
                label="📚 Download BibTeX",
                data=bibtex_content,
                file_name=f"pubmed_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.bib",
                mime="application/x-bibtex",
                use_container_width=True,
                key="persistent_bibtex_download"
            )

        # Display the dataframe
        st.dataframe(  
            df,  
            use_container_width=True,  
            hide_index=True,  
            column_config={  
                "Title": st.column_config.TextColumn("Title", width="large"),  
                "Abstract": st.column_config.TextColumn("Abstract", width="large"),  
                "DOI": st.column_config.LinkColumn("DOI", width="medium"),  
                "Publication Date": st.column_config.DateColumn("Publication Date", width="small"),  
                "Journal": st.column_config.TextColumn("Journal", width="small"),  
                "Source": st.column_config.TextColumn("Source", width="small")  
            }  
        )

        # Clear results button
        col1, col2, col3 = st.columns(3)
        with col2:
            if st.button("🗑️ Clear Results", use_container_width=True):
                st.session_state.search_results = None
                st.session_state.search_completed = False
                if hasattr(st.session_state, 'search_timestamp'):
                    delattr(st.session_state, 'search_timestamp')
                st.rerun()

    elif st.session_state.get('search_completed'):
        st.info("📭 No results found from your last search. Try refining your search criteria.")

    # Enhanced subscription confirmation with progress tracking  
    if subscribe_button:  
        if not subscriber_email:  
            st.error("❌ Please provide an email address.")  
        elif not selected_journals and not include_preprints:  
            st.error("❌ Please select at least one journal or include preprints.")  
        else:  
            # ========================================  
            # SUBSCRIPTION PROGRESS TRACKING  
            # ========================================  
            
            # Create subscription progress container  
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
                # Step 1: Validate and format journals  
                sub_status_text.info("🔍 Validating subscription parameters...")  
                sub_progress_bar.progress(10)  
                
                formatted_journals = validate_and_format_journals(selected_journals, full_to_abbrev) or [] if selected_journals else []
                
                # Step 2: Execute search for subscription  
                sub_status_text.info("🔍 Executing search for subscription...")  
                sub_progress_bar.progress(20)  
                
                # Check if we have existing search results  
                if 'df' in locals() and not df.empty:  
                    # User just performed a search - use those results  
                    csv_bytes = df.to_csv(index=False).encode("utf-8")  
                    has_results = True  
                    result_count = len(df)  
                    sub_status_text.success(f"✅ Using current search results ({result_count} articles)")  
                    sub_progress_bar.progress(50)  
                else:  
                    # User is subscribing without searching - execute search now  
                    sub_status_text.info("🔍 Running fresh search for subscription... Don't close the app")  
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

                # Step 3: Store subscription  
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
                    
                    # Step 4: Generate tokens  
                    sub_status_text.info("🔐 Generating secure tokens...")  
                    sub_progress_bar.progress(80)  
                    
                    unsubscribe_token = generate_unsubscribe_token(subscription_id)  
                    
                    if unsubscribe_token:  
                        unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"  
                    else:  
                        st.error("❌ Failed to generate unsubscribe token")  
                        st.stop()  

                    # Step 5: Prepare email content  
                    sub_status_text.info("📧 Preparing confirmation email...")  
                    sub_progress_bar.progress(85)  
                    
                    # Build source description  
                    source_list = []
                    if formatted_journals:
                        source_list.extend(formatted_journals)
                    if include_preprints:
                        source_list.extend(['biorxiv', 'medrxiv'])
                    source_description = ', '.join(source_list) if source_list else 'No sources selected'

                    # Step 6: Handle results and email sending
                    if has_results and csv_bytes is not None:
                        sub_status_text.info("📥 Generating download link...")
                        sub_progress_bar.progress(90)
                        
                        # Generate download token for results
                        download_token = generate_download_token(
                            csv_bytes, 
                            subscriber_email, 
                            subscription_id=subscription_id
                        )
                        
                        if download_token:
                            download_link = f"{BASE_URL}?token={download_token}&action=download"
                            
                            # Email content with results
                            email_body = f"""Hi {subscriber_email},

You have successfully subscribed to automatic PubMed updates.

📊 SUBSCRIPTION DETAILS:
📘 Journals: {', '.join(formatted_journals) if formatted_journals else 'None'}
📑 Preprints: {('bioRxiv, medRxiv' if include_preprints else 'None')}
🔍 All Sources: {source_description}
🔑 Keywords: {raw_keywords or 'None'}
🔁 Frequency: {frequency}

📥 YOUR CURRENT RESULTS ({result_count} articles found):
Your search results are available for download (expires in 3 hours):
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
                        # Email content without results
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
                    
                    # Step 7: Send email
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
                        
                        # Final completion
                        sub_progress_bar.progress(100)
                        sub_status_text.success("✅ Subscription confirmed and email sent!")
                        
                        # Show final statistics
                        sub_stats.metric(
                            label="Subscription Status",
                            value="Active",
                            delta=f"{result_count} articles" if has_results else "Ready for updates"
                        )
                        
                        st.success(success_message)
                        
                    except Exception as e:
                        sub_status_text.error(f"❌ Email sending failed: {e}")
                        sub_progress_bar.progress(100)
                        st.error(f"❌ Subscription created but email failed: {e}")
                
                elif result["status"] == "error":
                    sub_status_text.error(f"❌ Subscription failed: {result['message']}")
                    sub_progress_bar.progress(100)
                    # st.error(f"❌ {result['message']}")
                    # if "Maximum of 3" in result["message"]:
                    #     st.info("💡 **Tip:** You can manage your existing subscriptions by using the unsubscribe link in any of your previous emails.")
                else:
                    sub_status_text.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")
                    sub_progress_bar.progress(100)
                    st.error(f"❌ Subscription failed: {result.get('message', 'Unknown error')}")
            
            except Exception as e:
                sub_status_text.error(f"❌ Unexpected error: {e}")
                sub_progress_bar.progress(100)
                st.error(f"❌ Subscription failed: {e}")

    # Footer
    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center; color: #666; font-size: 0.8em; padding: 20px 0;'>
            Please report issues to <a href='https://github.com/rzhan186/Journal_tracker/issues' target='_blank'>GitHub</a>
        </div>
        """, 
        unsafe_allow_html=True
)

