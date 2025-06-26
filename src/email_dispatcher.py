from supabase import create_client
from dotenv import load_dotenv
import os
import smtplib
from email.message import EmailMessage
import logging

# Load secrets
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")  
EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")  # Use your email
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")  # Use your app password if 2FA is enabled

# Check if email credentials are loaded
if not EMAIL_ADDRESS or not EMAIL_PASSWORD:
    raise ValueError("❌ EMAIL_ADDRESS or EMAIL_PASSWORD not found in .env file")

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Initialize Supabase client
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

def send_email(to_email, subject, body):
    """Sends an email using SMTP."""
    try:
        if not to_email:
            raise ValueError("Recipient email cannot be None.")
        
        msg = EmailMessage()
        msg['Subject'] = subject
        msg['From'] = EMAIL_ADDRESS
        msg['To'] = to_email
        msg.set_content(body)

        # Connect to the server and send email
        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
            smtp.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
            smtp.send_message(msg)
        logging.info(f"Email sent to {to_email} successfully.")  # Logging email sent success
         
    except Exception as e:
        logging.error(f"Failed to send email to {to_email}: {str(e)}")  # Improved error logging
        raise  # Reraise exception for further handling if needed


######################################################################
# function to provide csv in the email

import base64
from datetime import datetime
from itsdangerous import URLSafeTimedSerializer

DOWNLOAD_SECRET = os.getenv("DOWNLOAD_SECRET")
download_serializer = URLSafeTimedSerializer(DOWNLOAD_SECRET, salt="csv-download")

def generate_download_token(csv_data, email):
    """Generate a secure token for CSV download that expires in 24 hours"""
    # Encode CSV data and email info
    payload = {
        'csv_data': base64.b64encode(csv_data).decode('utf-8'),
        'email': email,
        'timestamp': datetime.now().isoformat()
    }
    return download_serializer.dumps(payload)

def get_csv_from_token(token):
    """Retrieve CSV data from token (with 24-hour expiration)"""
    try:
        # Token expires after 24 hours (86400 seconds)
        payload = download_serializer.loads(token, max_age=86400)
        csv_data = base64.b64decode(payload['csv_data'].encode('utf-8'))
        return csv_data, payload['email']
    except Exception as e:
        return None, None
    

