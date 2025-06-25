# store_subscription.py
import os
from dotenv import load_dotenv
from itsdangerous import URLSafeSerializer, BadSignature, SignatureExpired
import logging
from supabase import create_client

load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")

supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

def get_user_subscription(email):
    """Fetches a user's active subscription from the database."""
    try:
        response = supabase.table("subscriptions").select("*").eq("email", email).eq("active", True).execute()
        
        if response.error:
            logging.error(f"Error fetching subscription: {response.error}")
            return None

        if response.data:
            return response.data[0]  # Return the first subscription found
        else:
            return None  # No active subscription found

    except Exception as e:
        logging.exception("Error in get_user_subscription")
        return None


def verify_unsubscribe_token(token):
    """Verifies the unsubscribe token and returns the associated subscription data."""
    try:
        email = serializer.loads(token)  # Decodes the token to get the email
        subscription = get_user_subscription(email)  # Fetch subscription by email

        if subscription:
            return subscription
        else:
            return None  # No subscription found for the email

    except (BadSignature, SignatureExpired):
        logging.warning("Invalid or expired token")
        return None
    except Exception as e:
        logging.exception("Error in verify_unsubscribe_token")
        return None