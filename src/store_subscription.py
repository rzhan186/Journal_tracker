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

def store_user_subscription(email, journals, keywords, start_date, end_date, frequency):
    """Stores a user's subscription details in the Supabase database."""
    try:
        # Generate the unsubscribe token
        unsubscribe_token = serializer.dumps(email)  # Create a token for this email

        data = {
            'email': email,
            'journals': journals,
            'keywords': keywords,
            'start_date': start_date,
            'end_date': end_date,
            'frequency': frequency,
            'active': True,
            'unsubscribe_token': unsubscribe_token,  # Store the token
        }
        
        response = supabase.table("subscriptions").insert(data).execute()
        logging.info(f"Supabase response: {response}")

        if hasattr(response, 'data') and response.data is not None:
            # Successfully stored
            logging.info(f"Subscription stored successfully for email: {email}")
            return {
                "status": "success", 
                "data": response.data,
                "unsubscribe_token": unsubscribe_token  # <== this must be here
            }
        elif hasattr(response, 'error'):
            logging.error(f"Error storing subscription: {response.error}")
            return {"status": "error", "message": str(response.error)}
        else:
            logging.error("Unexpected response structure from Supabase")
            return {"status": "error", "message": "Unexpected response structure"}

    except Exception as e:
        logging.exception("Error storing subscription")
        return {"status": "error", "message": str(e)}
    


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