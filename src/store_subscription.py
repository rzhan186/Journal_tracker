import os
from supabase import create_client
from itsdangerous import URLSafeSerializer, BadSignature
from dotenv import load_dotenv
import logging

# Load environment variables
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")

if not UNSUBSCRIBE_SECRET:
    raise ValueError("❌ UNSUBSCRIBE_SECRET not set in .env")

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

def store_user_subscription(email, journals, keywords, start_date, end_date, frequency):
    """Store user subscription in the Supabase database."""
    # Generate unsubscribe token
    unsubscribe_token = generate_unsubscribe_token(email, journals, keywords, frequency)
    
    # Insert into Supabase
    response = supabase.table("subscriptions").insert({
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "frequency": frequency,
        "unsubscribe_token": unsubscribe_token,  # Save the token
        "active": True  # Active subscription status
    }).execute()

    logging.info("Subscription stored successfully for email: %s", email)
    return response

def verify_unsubscribe_token(token):
    """Verifies the provided unsubscribe token."""
    try:
        return serializer.loads(token)
    except BadSignature:
        return None