# scheduler.py
import os
import logging
from datetime import datetime, timedelta
from dotenv import load_dotenv
import pandas as pd
import io
import time

# Your existing imports
from store_subscription import supabase
from email_dispatcher import send_email, generate_download_token
from tracking_main import (
    fetch_pubmed_articles_by_date,
    fetch_preprints,
    format_boolean_keywords_for_pubmed,
    merge_and_highlight_articles,
    standardize_date_format,
    standardize_doi_format
)

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

BASE_URL = "https://journaltracker.streamlit.app"

def get_days_from_frequency(frequency):
    """Convert frequency string to number of days"""
    if frequency == "weekly":
        return 7
    elif frequency == "monthly":
        return 30
    elif frequency.startswith("every"):
        # Extract number from "every X days"
        try:
            return int(frequency.split()[1])
        except:
            return 7  # fallback
    else:
        return 7  # default

def is_subscription_due(subscription):
    """Check if a subscription is due for an update"""
    try:
        last_sent = subscription.get('last_sent')
        frequency = subscription.get('frequency', 'weekly')
        
        if not last_sent:
            # Never sent before - send now
            return True
        
        # Parse last sent date
        last_sent_date = datetime.fromisoformat(last_sent.replace('Z', '+00:00'))
        days_interval = get_days_from_frequency(frequency)
        
        # Check if enough time has passed
        next_send = last_sent_date + timedelta(days=days_interval)
        return datetime.now() >= next_send
        
    except Exception as e:
        logging.error(f"Error checking subscription due date: {e}")
        return False

def execute_search_for_subscription(subscription):
    """Execute search based on subscription parameters"""
    try:
        # Extract subscription parameters
        journals = subscription.get('journals', [])
        keywords = subscription.get('keywords')
        include_preprints = subscription.get('include_preprints', False)
        frequency = subscription.get('frequency', 'weekly')
        
        # Calculate date range
        days_back = get_days_from_frequency(frequency)
        today = datetime.now().date()
        start_date = str(today - timedelta(days=days_back))
        end_date = str(today)
        
        logging.info(f"Searching for {subscription['email']}: {start_date} to {end_date}")
        
        all_articles = []
        
        # Search PubMed journals
        if journals:
            for journal in journals:
                articles = fetch_pubmed_articles_by_date(
                    journal, start_date, end_date,
                    format_boolean_keywords_for_pubmed(keywords) if keywords else None
                )
                for article in articles:
                    article["Journal"] = journal
                    article["Source"] = "PubMed"
                all_articles.extend(articles)
        
        # Search preprints
        if include_preprints:
            for server in ["biorxiv", "medrxiv"]:
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
        
        if all_articles:
            # Process results
            all_articles = standardize_date_format(all_articles)
            all_articles = standardize_doi_format(all_articles)
            merged = merge_and_highlight_articles(all_articles, [], keywords)
            return pd.DataFrame(merged)
        else:
            return pd.DataFrame()
            
    except Exception as e:
        logging.error(f"Error executing search for subscription {subscription['id']}: {e}")
        return pd.DataFrame()

# def send_subscription_email(subscription, results_df):
#     """Send subscription update email"""
#     try:
#         email = subscription['email']
        
#         # Generate unsubscribe link (you'll need to import your serializer)
#         from itsdangerous import URLSafeTimedSerializer
        
#         UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
#         serializer = URLSafeTimedSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")
        
#         unsubscribe_data = {
#             'email': email,
#             'timestamp': datetime.now().isoformat()
#         }
#         unsubscribe_token = serializer.dumps(unsubscribe_data)
#         unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"
        
#         # Prepare email content
#         journals = subscription.get('journals', [])
#         keywords = subscription.get('keywords', 'None')
#         include_preprints = subscription.get('include_preprints', False)
#         frequency = subscription.get('frequency', 'weekly')
        
#         if not results_df.empty:
#             # Generate CSV and download link
#             csv_bytes = results_df.to_csv(index=False).encode("utf-8")
#             download_token = generate_download_token(csv_bytes, email)
#             download_link = f"{BASE_URL}?token={download_token}&action=download"
#             result_count = len(results_df)
            
#             subject = f"📬 Journal Tracker: {result_count} New Articles Found"
            
#             email_body = f"""Hi {email},

# Your {frequency} PubMed update is ready!

# 📊 SEARCH PARAMETERS:
# 📘 Journals: {', '.join(journals) if journals else 'None'}
# 📑 Preprints: {'bioRxiv, medRxiv' if include_preprints else 'None'}
# 🔑 Keywords: {keywords}

# 📥 YOUR RESULTS ({result_count} articles found):
# Download your results (expires in 24 hours):
# 🔗 {download_link}

# 🔓 UNSUBSCRIBE: {unsubscribe_link}

# – PubMed Tracker Team
#             """
#         else:
#             # No results found - don't send email
#             logging.info(f"No results for {email}, skipping email")
#             return False
        
#         # Send email
#         send_email(email, subject, email_body)
#         logging.info(f"Subscription email sent to {email}")
#         return True
        
#     except Exception as e:
#         logging.error(f"Error sending subscription email: {e}")
#         return False

def send_subscription_email(subscription, results_df):
    """Send subscription update email with improved formatting and error handling"""
    try:
        email = subscription['email']
        
        # Generate unsubscribe link
        from itsdangerous import URLSafeTimedSerializer
        
        UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")
        serializer = URLSafeTimedSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")
        
        unsubscribe_data = {
            'email': email,
            'timestamp': datetime.now().isoformat()
        }
        unsubscribe_token = serializer.dumps(unsubscribe_data)
        unsubscribe_link = f"{BASE_URL}?token={unsubscribe_token}"
        
        # Prepare email content
        journals = subscription.get('journals', [])
        keywords = subscription.get('keywords', 'None')
        include_preprints = subscription.get('include_preprints', False)
        frequency = subscription.get('frequency', 'weekly')
        
        # ✅ IMPROVED: Always update last_sent, even if no results
        if not results_df.empty:
            # Generate CSV and download link
            csv_bytes = results_df.to_csv(index=False).encode("utf-8")
            download_token = generate_download_token(csv_bytes, email)
            download_link = f"{BASE_URL}?token={download_token}&action=download"
            result_count = len(results_df)
            
            # ✅ IMPROVED: Better subject line
            subject = f"📚 Journal Tracker: {result_count} New Article{'s' if result_count != 1 else ''} Found"
            
            # ✅ IMPROVED: Professional HTML-formatted email
            email_body = f"""Hello!

Your {frequency} Journal Tracker update is ready with {result_count} new article{'s' if result_count != 1 else ''}.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
📊 SUBSCRIPTION DETAILS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

📘 Journals: {', '.join(journals) if journals else 'None selected'}
📑 Preprints: {'bioRxiv, medRxiv' if include_preprints else 'None'}
🔍 Keywords: {keywords if keywords and keywords != 'None' else 'None'}
🔄 Frequency: {frequency.title()}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
📥 YOUR RESULTS ({result_count} articles)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Your research articles are ready for download:
🔗 {download_link}

⏰ Download expires in 24 hours

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

📧 Questions? Reply to this email
🔓 Unsubscribe: {unsubscribe_link}

Best regards,
The Journal Tracker Team
journaltracker.site"""

            # ✅ IMPROVED: Better error handling for email sending
            try:
                send_email(email, subject, email_body)
                logging.info(f"✅ Subscription email sent to {email} ({result_count} articles)")
                return True
            except Exception as email_error:
                logging.error(f"❌ Failed to send email to {email}: {email_error}")
                return False
                
        else:
            # ✅ IMPROVED: Still update last_sent even if no results to prevent constant checking
            logging.info(f"📭 No results for {email}, updating last_sent anyway to prevent constant checking")
            return True  # Changed from False to True
        
    except Exception as e:
        logging.error(f"❌ Error in send_subscription_email for {subscription.get('email', 'unknown')}: {e}")
        return False





def update_subscription_last_sent(subscription_id):
    """Update the last_sent timestamp for a subscription"""
    try:
        supabase.table("subscriptions").update({
            "last_sent": datetime.now().isoformat()
        }).eq("id", subscription_id).execute()
        
        logging.info(f"Updated last_sent for subscription {subscription_id}")
        
    except Exception as e:
        logging.error(f"Error updating last_sent: {e}")

# def process_due_subscriptions():
#     """Main function to process all due subscriptions"""
#     try:
#         logging.info("Starting subscription processing...")
        
#         # Get all active subscriptions
#         response = supabase.table("subscriptions").select("*").eq("active", True).execute()
        
#         total_subscriptions = len(response.data)
#         processed = 0
#         sent = 0
        
#         logging.info(f"Found {total_subscriptions} active subscriptions")
        
#         for subscription in response.data:
#             try:
#                 if is_subscription_due(subscription):
#                     logging.info(f"Processing subscription for {subscription['email']}")
                    
#                     # Execute search
#                     results = execute_search_for_subscription(subscription)
                    
#                     # Send email if results found
#                     if send_subscription_email(subscription, results):
#                         sent += 1
                        
#                     # Update last_sent timestamp
#                     update_subscription_last_sent(subscription['id'])
#                     processed += 1
#                 else:
#                     logging.info(f"Subscription for {subscription['email']} not due yet")
                    
#             except Exception as e:
#                 logging.error(f"Error processing subscription {subscription['id']}: {e}")
#                 continue
        
#         logging.info(f"Processing complete: {processed} processed, {sent} emails sent")
        
#     except Exception as e:
#         logging.error(f"Error in process_due_subscriptions: {e}")


def process_due_subscriptions(dry_run=False):
    """Main function to process all due subscriptions"""
    try:
        mode = "🧪 DRY RUN" if dry_run else "📬 LIVE MODE"
        logging.info(f"Starting subscription processing ({mode})...")
        
        response = supabase.table("subscriptions").select("*").eq("active", True).execute()
        
        total_subscriptions = len(response.data)
        processed = 0
        sent = 0
        
        logging.info(f"Found {total_subscriptions} active subscriptions")
        
        for subscription in response.data:
            try:
                if is_subscription_due(subscription):
                    email = subscription['email']
                    frequency = subscription.get('frequency', 'weekly')
                    
                    if dry_run:
                        logging.info(f"🧪 DRY RUN: Would send {frequency} update to {email}")
                        processed += 1
                        continue  # Skip actual processing
                    
                    # ✅ Only execute real processing in live mode
                    logging.info(f"📬 Processing subscription for {email}")
                    
                    # Execute search
                    results = execute_search_for_subscription(subscription)
                    
                    # Send email if results found
                    if send_subscription_email(subscription, results):
                        sent += 1
                        
                    # Update last_sent timestamp
                    update_subscription_last_sent(subscription['id'])
                    processed += 1
                else:
                    if not dry_run:  # Only log in live mode to reduce noise
                        logging.info(f"Subscription for {subscription['email']} not due yet")
                        
            except Exception as e:
                logging.error(f"Error processing subscription {subscription['id']}: {e}")
                continue
        
        if dry_run:
            logging.info(f"🧪 DRY RUN COMPLETE: Would process {processed} subscriptions")
        else:
            logging.info(f"📬 LIVE PROCESSING COMPLETE: {processed} processed, {sent} emails sent")
        
    except Exception as e:
        logging.error(f"Error in process_due_subscriptions: {e}")

# dry run mode
if __name__ == "__main__":
    import sys
    dry_run = "--dry-run" in sys.argv or "-d" in sys.argv
    process_due_subscriptions(dry_run=dry_run)