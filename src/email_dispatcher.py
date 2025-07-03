# email_dispatcher.py

from supabase import create_client
from dotenv import load_dotenv
import os
import logging
import sib_api_v3_sdk
from sib_api_v3_sdk.rest import ApiException
import streamlit as st
import uuid
from datetime import datetime, timedelta

# Load secrets
load_dotenv()

# Supabase configuration (keep existing)
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# ✅ NEW: Brevo email configuration
def get_email_config():
    """Get email configuration from environment or Streamlit secrets"""
    try:
        # Try Streamlit secrets first (for deployed app)
        return {
            'api_key': st.secrets["BREVO_API_KEY"],
            'sender_email': st.secrets["BREVO_SENDER_EMAIL"],
            'sender_name': st.secrets["BREVO_SENDER_NAME"]
        }
    except:
        # Fall back to environment variables (for local development)
        return {
            'api_key': os.getenv("BREVO_API_KEY"),
            'sender_email': os.getenv("BREVO_SENDER_EMAIL"),
            'sender_name': os.getenv("BREVO_SENDER_NAME")
        }

# Validate email configuration
email_config = get_email_config()
if not all(email_config.values()):
    raise ValueError("❌ BREVO_API_KEY, BREVO_SENDER_EMAIL, or BREVO_SENDER_NAME not found")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Initialize Supabase client
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

def send_email(to_email, subject, body, sender_email=None, sender_name=None, api_key=None):
    """
    UPDATED: Send email using Brevo API with new domain
    """
    try:
        if not to_email:
            raise ValueError("Recipient email cannot be None.")
        
        # Use provided parameters or fall back to configuration
        config = get_email_config()
        api_key = api_key or config['api_key']
        sender_email = sender_email or config['sender_email']
        sender_name = sender_name or config['sender_name']
        
        # Configure Brevo API
        configuration = sib_api_v3_sdk.Configuration()
        configuration.api_key['api-key'] = api_key
        api_instance = sib_api_v3_sdk.TransactionalEmailsApi(sib_api_v3_sdk.ApiClient(configuration))
        
        # Convert plain text body to HTML (preserve line breaks)
        html_content = body.replace('\n', '<br>\n')
        
        # Create email
        send_smtp_email = sib_api_v3_sdk.SendSmtpEmail(
            to=[{"email": to_email}],
            sender={
                "name": sender_name,
                "email": sender_email
            },
            subject=subject,
            text_content=body,        # Plain text version
            html_content=html_content # HTML version
        )
        
        # Send email
        api_response = api_instance.send_transac_email(send_smtp_email)
        logging.info(f"✅ Email sent to {to_email} successfully. Message ID: {api_response.message_id}")
        return True, api_response.message_id
        
    except ApiException as e:
        error_msg = f"Brevo API error: {e}"
        logging.error(f"❌ Failed to send email to {to_email}: {error_msg}")
        raise Exception(error_msg)
    except Exception as e:
        error_msg = f"Email sending failed: {str(e)}"
        logging.error(f"❌ Failed to send email to {to_email}: {error_msg}")
        raise Exception(error_msg)

######################################################################
# CSV download functions

import base64
from datetime import datetime
from itsdangerous import URLSafeTimedSerializer

DOWNLOAD_SECRET = os.getenv("DOWNLOAD_SECRET")
download_serializer = URLSafeTimedSerializer(DOWNLOAD_SECRET, salt="csv-download")

def generate_download_token(csv_data, email, subscription_id=None):
    """Generate token and store CSV in file storage instead of database"""
    try:
        # Handle both bytes and string data
        if isinstance(csv_data, bytes):
            csv_bytes = csv_data
        elif isinstance(csv_data, str):
            csv_bytes = csv_data.encode('utf-8')
        else:
            csv_bytes = str(csv_data).encode('utf-8')
        
        # Generate a unique token
        token = str(uuid.uuid4()).replace('-', '')[:16]
        
        # Upload to Supabase Storage instead of storing in database
        file_path = f"csv_downloads/{token}.csv"
        
        try:
            # Upload to Supabase Storage
            supabase.storage.from_('csv-files').upload(
                file_path, 
                csv_bytes,
                file_options={"content-type": "text/csv"}
            )
            
            # Store only the token and metadata in database
            expires_at = (datetime.now() + timedelta(hours=3)).isoformat()
            
            # Store minimal data in database
            if subscription_id:
                supabase.table('subscriptions').update({
                    'csv_token': token,
                    'csv_file_path': file_path,  # Store file path instead of data
                    'csv_expires_at': expires_at
                }).eq('id', subscription_id).execute()
            
            logging.info(f"✅ CSV file uploaded to storage: {file_path}")
            return token
            
        except Exception as storage_error:
            logging.error(f"❌ Failed to upload to storage: {storage_error}")
            return None
            
    except Exception as e:
        logging.error(f"❌ Error generating download token: {str(e)}")
        return None

def get_csv_from_token(token):
    """Retrieve CSV data from file storage with automatic cleanup"""
    try:
        # Get file path from database
        result = supabase.table('subscriptions').select('csv_file_path, csv_expires_at, email').eq('csv_token', token).execute()
        
        if not result.data:
            return None, None
        
        record = result.data[0]
        
        # Check expiration
        if record.get('csv_expires_at'):
            expires_at = datetime.fromisoformat(record['csv_expires_at'])
            if datetime.now() > expires_at:
                # **CLEANUP: Remove expired file from storage**
                try:
                    supabase.storage.from_('csv-files').remove([record['csv_file_path']])
                    logging.info(f"🗑️ Cleaned up expired file: {record['csv_file_path']}")
                except Exception as cleanup_error:
                    logging.warning(f"⚠️ Failed to cleanup expired file: {cleanup_error}")
                
                # Clear database record
                supabase.table('subscriptions').update({
                    'csv_token': None,
                    'csv_file_path': None,
                    'csv_expires_at': None
                }).eq('csv_token', token).execute()
                
                return None, None
        
        # Download from storage
        file_path = record['csv_file_path']
        response = supabase.storage.from_('csv-files').download(file_path)
        
        return response, record['email']
        
    except Exception as e:
        logging.error(f"❌ Error retrieving CSV: {str(e)}")
        return None, None

def cleanup_expired_csv_files():
    """Scheduled cleanup function for expired CSV files"""
    try:
        current_time = datetime.now().isoformat()
        
        # Get all expired records
        expired_records = supabase.table('subscriptions').select('csv_file_path, csv_token').lt('csv_expires_at', current_time).execute()
        
        if expired_records.data:
            cleanup_count = 0
            for record in expired_records.data:
                if record.get('csv_file_path'):
                    try:
                        # Remove file from storage
                        supabase.storage.from_('csv-files').remove([record['csv_file_path']])
                        cleanup_count += 1
                    except Exception as e:
                        logging.warning(f"⚠️ Failed to cleanup {record['csv_file_path']}: {e}")
            
            # Clear database records
            supabase.table('subscriptions').update({
                'csv_token': None,
                'csv_file_path': None,
                'csv_expires_at': None
            }).lt('csv_expires_at', current_time).execute()
            
            logging.info(f"🗑️ Cleaned up {cleanup_count} expired CSV files")
        
    except Exception as e:
        logging.error(f"❌ Error in cleanup_expired_csv_files: {str(e)}")


def get_next_update_timeframe(frequency):
    """Convert frequency to human-readable timeframe"""
    if frequency == "weekly":
        return "1 week"
    elif frequency == "monthly":
        return "1 month"
    elif frequency.startswith("every"):
        days = frequency.split()[1]
        return f"{days} days"
    else:
        return "as scheduled"