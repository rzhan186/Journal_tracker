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
# CSV download functions (unchanged)

import base64
from datetime import datetime
from itsdangerous import URLSafeTimedSerializer

DOWNLOAD_SECRET = os.getenv("DOWNLOAD_SECRET")
download_serializer = URLSafeTimedSerializer(DOWNLOAD_SECRET, salt="csv-download")

# def generate_download_token(csv_data, email):
#     """Generate a secure token for CSV download that expires in 24 hours"""
#     # Handle both bytes and string data
#     if isinstance(csv_data, str):
#         csv_data = csv_data.encode('utf-8')
    
#     payload = {
#         'csv_data': base64.b64encode(csv_data).decode('utf-8'),
#         'email': email,
#         'timestamp': datetime.now().isoformat()
#     }
#     return download_serializer.dumps(payload)

# def get_csv_from_token(token):
#     """Retrieve CSV data from token (with 24-hour expiration)"""
#     try:
#         # Token expires after 24 hours (86400 seconds)
#         payload = download_serializer.loads(token, max_age=86400)
#         csv_data = base64.b64decode(payload['csv_data'].encode('utf-8'))
#         return csv_data, payload['email']
#     except Exception as e:
#         logging.error(f"Failed to retrieve CSV from token: {str(e)}")
#         return None, None

def generate_download_token(csv_data, email, subscription_id=None):
    """Generate a secure token for CSV download that expires in 24 hours"""
    try:
        # Handle both bytes and string data
        if isinstance(csv_data, str):
            csv_data = csv_data.encode('utf-8')
        elif isinstance(csv_data, bytes):
            pass  # Already bytes
        else:
            csv_data = str(csv_data).encode('utf-8')
        
        # Generate a short unique token
        token = str(uuid.uuid4())[:16]  # Short 16-character token
        
        # Calculate expiration time (24 hours from now)
        expires_at = (datetime.now() + timedelta(hours=24)).isoformat()
        
        # Store in downloads table (or use your preferred storage method)
        result = supabase.table('csv_downloads').insert({
            'token': token,
            'subscription_id': subscription_id,  # This can be None
            'csv_data': base64.b64encode(csv_data).decode('utf-8'),
            'email': email,
            'created_at': datetime.now().isoformat(),
            'expires_at': expires_at
        }).execute()
        
        if result.data:
            logging.info(f"✅ CSV download token generated: {token}")
            return token
        else:
            logging.error("❌ Failed to store CSV download data")
            return None
            
    except Exception as e:
        logging.error(f"❌ Error generating download token: {str(e)}")
        return None
    
def get_csv_from_token(token):
    """Retrieve CSV data from subscription table using token"""
    try:
        # Query subscription table for the token
        result = supabase.table('subscriptions').select('*').eq('csv_token', token).execute()
        
        if not result.data:
            logging.error(f"❌ Token not found: {token}")
            return None, None
        
        record = result.data[0]
        
        # Check if token has expired
        if record.get('csv_expires_at'):
            expires_at = datetime.fromisoformat(record['csv_expires_at'])
            if datetime.now() > expires_at:
                logging.error(f"❌ Token expired: {token}")
                # Clear expired token and data
                supabase.table('subscriptions').update({
                    'csv_token': None,
                    'current_csv_data': None,
                    'csv_expires_at': None
                }).eq('csv_token', token).execute()
                return None, None
        
        csv_data = record.get('current_csv_data')
        if not csv_data:
            logging.error(f"❌ No CSV data found for token: {token}")
            return None, None
            
        email = record['email']
        
        logging.info(f"✅ CSV data retrieved for token: {token}")
        return csv_data.encode('utf-8'), email
        
    except Exception as e:
        logging.error(f"❌ Error retrieving CSV from token: {str(e)}")
        return None, None

def cleanup_expired_csv_tokens():
    """Clean up expired CSV tokens and data"""
    try:
        current_time = datetime.now().isoformat()
        result = supabase.table('subscriptions').update({
            'csv_token': None,
            'current_csv_data': None,
            'csv_expires_at': None
        }).lt('csv_expires_at', current_time).execute()
        
        if result.data:
            logging.info(f"✅ Cleaned up {len(result.data)} expired CSV tokens")
    except Exception as e:
        logging.error(f"❌ Error cleaning up expired CSV tokens: {str(e)}")


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