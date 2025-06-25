import os
import logging
from supabase import create_client
from itsdangerous import URLSafeSerializer, BadSignature, SignatureExpired
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

# Set up logging
logging.basicConfig(level=logging.INFO)

def generate_unsubscribe_token(email, journals, keywords, frequency):
    """Generates a secure unsubscribe token for the user."""
    data = {
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "frequency": frequency,
    }
    token = serializer.dumps(data)
    logging.info("Generated unsubscribe token for email: %s", email)
    return token

def send_confirmation_email(email, unsubscribe_token):
    """Sends a confirmation email to the subscribed user."""
    subject = "Subscription Confirmation"
    unsubscribe_link = f"https://journaltracker.streamlit.app/unsubscribe/?token={unsubscribe_token}"  # Use your deployed app URL
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
        logging.error("❌ Failed to send email to %s: %s", email, e)

def store_user_subscription(email, journals, keywords, start_date, end_date, frequency):
    """Store user subscription in the Supabase database."""
    # Generate unsubscribe token
    unsubscribe_token = generate_unsubscribe_token(email, journals, keywords, frequency)

    # Insert into Supabase
    logging.info("Storing subscription for email: %s", email)
    response = supabase.table("subscriptions").insert({
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "frequency": frequency,
        "unsubscribe_token": unsubscribe_token,  # Save the token
        "active": True  # Active subscription status
    }).execute()

    # Checking for insertion success and handling response
    if response.data:  # Checking if we got data back, which means success
        send_confirmation_email(email, unsubscribe_token)
        logging.info("Subscription stored successfully for email: %s with ID: %s", email, response.data[0]['id'])
    else:
        # Handle error -- if response.data is None, an error occurred
        logging.error("❌ Error storing subscription for email %s: %s", email, response.error)  # Log the error message

    return response

def verify_unsubscribe_token(token):
    """Verifies the provided unsubscribe token."""
    try:
        return serializer.loads(token)
    except BadSignature:
        logging.error("❌ Invalid unsubscribe token: %s", token)
        return None
    
def get_user_subscription(email):
    """Fetch the user's subscription from the database based on the provided email."""
    logging.info(f"Fetching subscription for: {email}")
    response = supabase.table("subscriptions").select("*").eq("email", email).execute()
    
    # Use the response status code to check for errors
    if response.status_code != 200:
        logging.error(f"Error fetching subscription: {response.error}")  # Log the error message if any
        return None  # Return None on failure

    # Return the first subscription found if available, else return None
    return response.data[0] if response.data else None

def verify_unsubscribe_token(token):
    secret = os.getenv("UNSUBSCRIBE_SECRET")
    serializer = URLSafeSerializer(secret, salt="unsubscribe")
    try:
        user_info = serializer.loads(token)  # Decodes the token
        # Fetch the user's subscription from your database (using the user's email address)
        # Return subscription details if found
        return get_user_subscription(user_info['email'])  # Implement get_user_subscription appropriately
    except (BadSignature, SignatureExpired):
        return None  # Token is invalid or expired