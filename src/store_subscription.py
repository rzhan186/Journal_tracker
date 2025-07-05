# store_subscription.py

from supabase import create_client
import os
from dotenv import load_dotenv
import logging
import pandas as pd

# Load environment variables
load_dotenv()

# Initialize Supabase client
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

def store_user_subscription(email, journals, keywords, frequency, include_preprints=False):
    """Store subscription in Supabase database with limit check"""
    try:
        # Check existing active subscriptions for this email
        existing_response = supabase.table("subscriptions").select("id").eq("email", email).eq("active", True).execute()
        
        if existing_response.data and len(existing_response.data) >= 3:
            return {
                "status": "error",
                "data": None,
                "message": "Maximum of 3 active subscriptions per email address. Please unsubscribe from some subscriptions before creating new ones."
            }
        
        # Prepare data for insertion
        subscription_data = {
            "email": email,
            "journals": journals,
            "keywords": keywords,
            "frequency": frequency,
            "include_preprints": include_preprints,
            "active": True,
            "last_sent": None
        }
        
        # Insert into database
        response = supabase.table("subscriptions").insert(subscription_data).execute()
        
        if response.data:
            return {
                "status": "success",
                "data": response.data[0],
                "message": "Subscription created successfully"
            }
        else:
            logging.error("Failed to store subscription")
            return {
                "status": "error",
                "data": None,
                "message": "Failed to store subscription"
            }
            
    except Exception as e:
        logging.error(f"Error storing subscription: {e}")
        return {
            "status": "error",
            "data": None,
            "message": f"Error storing subscription: {str(e)}"
        }

def get_user_subscriptions_by_email(email):
    """Get all active subscriptions for a user by email"""
    try:
        response = supabase.table("subscriptions").select("*").eq("email", email).eq("active", True).execute()
        return response.data if response.data else []
    except Exception as e:
        logging.error(f"Error getting subscriptions for {email}: {e}")
        return []

def verify_unsubscribe_token(token):
    """Verify unsubscribe token and return subscription data - NO EXPIRATION"""
    try:
        from itsdangerous import URLSafeSerializer
        import os
        import logging
        from supabase import create_client
        
        # Initialize supabase client
        SUPABASE_URL = os.getenv("SUPABASE_URL")
        SUPABASE_KEY = os.getenv("SUPABASE_KEY") 
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        
        # Get the secret key
        secret_key = os.getenv("UNSUBSCRIBE_SECRET")
        if not secret_key:
            logging.error("UNSUBSCRIBE_SECRET not found in environment variables")
            return None
        
        # Create serializer WITHOUT TimestampSigner (no expiration)
        serializer = URLSafeSerializer(secret_key)
        
        # Decode the token to get subscription_id
        subscription_id = serializer.loads(token)
        logging.info(f"Decoded subscription ID: {subscription_id}")
        
        # Get subscription from database
        response = supabase.table("subscriptions").select("*").eq("id", subscription_id).execute()
        
        if response.data and len(response.data) > 0:
            subscription_data = response.data[0]  # Get the first (and should be only) result
            logging.info(f"Found subscription for email: {subscription_data.get('email', 'Unknown')}")
            return subscription_data
        else:
            logging.warning(f"No subscription found for ID: {subscription_id}")
            return None
            
    except Exception as e:
        logging.error(f"Error verifying unsubscribe token: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None

def generate_unsubscribe_token(subscription_id):
    """Generate a permanent unsubscribe token for a subscription"""
    try:
        from itsdangerous import URLSafeSerializer
        import os
        import logging
        
        secret_key = os.getenv("UNSUBSCRIBE_SECRET")
        if not secret_key:
            raise ValueError("UNSUBSCRIBE_SECRET not found")
        
        # Create serializer WITHOUT TimestampSigner (no expiration)
        serializer = URLSafeSerializer(secret_key) 
        
        # Generate token
        token = serializer.dumps(subscription_id)
        return token
        
    except Exception as e:
        logging.error(f"Error generating unsubscribe token: {e}")
        return None

def get_active_subscriptions():
    """Get all active subscriptions from database"""
    try:
        response = supabase.table("subscriptions").select("*").eq("active", True).execute()
        return response.data if response.data else []
    except Exception as e:
        logging.error(f"Error getting active subscriptions: {e}")
        return []

def update_subscription_last_sent(subscription_id):
    """Update the last_sent timestamp for a subscription"""
    try:
        from datetime import datetime
        
        response = supabase.table("subscriptions").update({
            "last_sent": datetime.now().isoformat()
        }).eq("id", subscription_id).execute()
        
        return response.data is not None
    except Exception as e:
        logging.error(f"Error updating last_sent for subscription {subscription_id}: {e}")
        return False

def deactivate_subscription(subscription_id):
    """Deactivate a subscription"""
    try:
        response = supabase.table("subscriptions").update({
            "active": False
        }).eq("id", subscription_id).execute()
        
        return response.data is not None
    except Exception as e:
        logging.error(f"Error deactivating subscription {subscription_id}: {e}")
        return False
    
