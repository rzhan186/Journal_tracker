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



# def process_subscriptions(): # this function is not intendted for use in app.py 
#     """Processes subscriptions by fetching from the Supabase and sending emails."""
#     logging.info("📨 Checking subscriptions...")
#     try:
#         # Execute the query to retrieve subscriptions
#         res = supabase.table("subscriptions").select("*").execute()

#         # Check if the response is successful
#         if not res.data:
#             logging.warning("No subscriptions found.")  # Log if no data is returned
#             return
        
#         # Process each subscription
#         for user in res.data:
#             email = user['email']
#             journals = user['journals']
#             frequency = user['frequency']
#             body = f"""Hi {email},

# This is your {frequency} update for the following journals: {journals}.

# 🔎 This is just a test. The real search and article results would be inserted here later.

# – PubMed Tracker Team
#             """
#             send_email(email, f"Your {frequency} PubMed update", body)  # Send the email
#             logging.info(f"✅ Sent test email to {email}")  # Log sent emails
            
#     except Exception as e:
#         logging.error(f"Failed to process subscriptions: {str(e)}")  # Improved error logging

# if __name__ == "__main__":
#     process_subscriptions()  # Call the function to process subscriptions