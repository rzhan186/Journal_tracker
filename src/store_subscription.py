import os
import logging
from supabase import create_client
from itsdangerous import URLSafeSerializer, BadSignature
from dotenv import load_dotenv
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

# Load environment variables
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")

EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")

if not UNSUBSCRIBE_SECRET or not EMAIL_PASSWORD:
    raise ValueError("❌ Ensure UNSUBSCRIBE_SECRET and EMAIL_PASSWORD are set in .env")

# Initialize Supabase client and URL serializer
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

def generate_unsubscribe_token(email, journals, keywords, frequency):
    """Generates a secure unsubscribe token for the user."""
    data = {
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "frequency": frequency,
    }
    return serializer.dumps(data)

def send_confirmation_email(email, unsubscribe_token):
    """Sends a confirmation email to the subscribed user."""
    subject = "Subscription Confirmation"
    unsubscribe_link = f"https://journaltracker.streamlit.app/?token={unsubscribe_token}"  # Use your deployed app URL
    body = f"""
    Thank you for your subscription!

    You will receive updates based on your preferences.

    If you'd like to unsubscribe, click the link below:
    {unsubscribe_link}
    """

    msg = MIMEMultipart()
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = email
    msg['Subject'] = subject
    msg.attach(MIMEText(body, 'plain'))

    try:
        with smtplib.SMTP('smtp.gmail.com', 587) as server:
            server.starttls()
            server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
            server.send_message(msg)
            logging.info("✅ Confirmation email sent successfully to %s", email)
    except Exception as e:
        logging.error("❌ Failed to send email: %s", e)

def store_user_subscription(email, journals, keywords, start_date, end_date, frequency):
    """Store user subscription in the Supabase database."""
    # Generate unsubscribe token
    unsubscribe_token = generate_unsubscribe_token(email, journals, keywords, frequency)
    logging.info("Generated unsubscribe token: %s", unsubscribe_token)  # Log the generated token

    # Insert into Supabase
    response = supabase.table("subscriptions").insert({
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "frequency": frequency,
        "unsubscribe_token": unsubscribe_token,  # Save the token
        "active": True  # Active subscription status
    }).execute()

    logging.info("Attempting to store subscription for email: %s", email)

    # Checking for insertion success and handling response
    if response.data:  # Checking if we got data back, which means success
        send_confirmation_email(email, unsubscribe_token)
        logging.info("Subscription stored successfully with ID: %s", response.data[0]['id'])
    else:
        # Handle error -- if response.data is None, an error occurred
        logging.error("❌ Error storing subscription: %s", response.error)  # Use response.error directly as a fallback if there's a response but no data.

    return response

def verify_unsubscribe_token(token):
    """Verifies the provided unsubscribe token."""
    try:
        return serializer.loads(token)
    except BadSignature:
        return None