# store_subscription.py

# import os
# from dotenv import load_dotenv
# from itsdangerous import URLSafeSerializer, BadSignature, SignatureExpired
# import logging
# from supabase import create_client
# from datetime import datetime

# load_dotenv()

# SUPABASE_URL = os.getenv("SUPABASE_URL")
# SUPABASE_KEY = os.getenv("SUPABASE_KEY")
# UNSUBSCRIBE_SECRET = os.getenv("UNSUBSCRIBE_SECRET")

# supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
# serializer = URLSafeSerializer(UNSUBSCRIBE_SECRET, salt="unsubscribe")

# def store_user_subscription(email, journals, keywords, frequency, include_preprints=False):
#     """
#     Store user subscription with proper JSON handling
#     """
#     try:
#         # Ensure journals is stored as a proper list
#         if isinstance(journals, str):
#             journals_list = [j.strip() for j in journals.split(',') if j.strip()]
#         elif isinstance(journals, list):
#             journals_list = journals
#         else:
#             journals_list = []

#         subscription_data = {
#             'email': email,
#             'journals': journals_list,  # This will be stored as JSON
#             'keywords': keywords,
#             'frequency': frequency,
#             #'start_date': start_date,
#             #'end_date': end_date,
#             'include_preprints': include_preprints,  # Boolean
#             'active': True,
#             'created_at': datetime.now().isoformat()    
#         }

#         # Insert into Supabase
#         response = supabase.table("subscriptions").insert(subscription_data).execute()
        
#         if response.data:
#             return {
#                 "status": "success",
#                 "message": "Subscription stored successfully",
#                 "subscription_id": response.data[0]['id']
#             }
#         else:
#             return {
#                 "status": "error", 
#                 "message": "Failed to store subscription"
#             }
            
#     except Exception as e:
#         logging.error(f"Error storing subscription: {e}")
#         return {
#             "status": "error",
#             "message": f"Database error: {str(e)}"
#         }




# def get_user_subscription(email):
#     """Fetches a user's active subscription from the database."""
#     try:
#         response = supabase.table("subscriptions").select("*").eq("email", email).eq("active", True).execute()
        
#         if response.error:
#             logging.error(f"Error fetching subscription: {response.error}")
#             return None

#         if response.data:
#             return response.data[0]  # Return the first subscription found
#         else:
#             return None  # No active subscription found

#     except Exception as e:
#         logging.exception("Error in get_user_subscription")
#         return None


# def get_user_subscriptions_by_email(email):
#     """
#     Retrieve all active subscriptions for a given email address
#     """
#     try:
#         response = supabase.table("subscriptions").select("*") \
#             .eq("email", email) \
#             .eq("active", True) \
#             .order("created_at", desc=True) \
#             .execute()
        
#         if response.data:
#             return response.data
#         else:
#             return []
            
#     except Exception as e:
#         logging.error(f"Error retrieving subscriptions for {email}: {e}")
#         return []

# # def verify_unsubscribe_token(token):
# #     """
# #     Verify unsubscribe token and return subscription details
# #     Updated to return user info instead of specific subscription
# #     """
# #     try:
# #         # Decode the token to get the email
# #         data = serializer.loads(token)
# #         email = data.get('email')
        
# #         if not email:
# #             return None
        
# #         # Return a dict with the email for compatibility
# #         return {'email': email}
        
# #     except Exception as e:
# #         logging.error(f"Error verifying unsubscribe token: {e}")
# #         return None


# def verify_unsubscribe_token(token):
#     """Verify unsubscribe token and return subscription data - NO EXPIRATION"""
#     try:
#         from itsdangerous import URLSafeSerializer
#         import os
        
#         # Get the secret key
#         secret_key = os.getenv("UNSUBSCRIBE_SECRET")
#         if not secret_key:
#             logging.error("UNSUBSCRIBE_SECRET not found in environment variables")
#             return None
        
#         # Create serializer WITHOUT TimestampSigner (no expiration)
#         serializer = URLSafeSerializer(secret_key)
        
#         # Decode the token
#         subscription_id = serializer.loads(token)
        
#         # Get subscription from database
#         response = supabase.table("subscriptions").select("*").eq("id", subscription_id).single().execute()
        
#         if response.data:
#             return response.data
#         else:
#             logging.warning(f"No subscription found for token: {subscription_id}")
#             return None
            
#     except Exception as e:
#         logging.error(f"Error verifying unsubscribe token: {e}")
#         return None


# def generate_unsubscribe_token(subscription_id):
#     """Generate a permanent unsubscribe token for a subscription"""
#     try:
#         from itsdangerous import URLSafeSerializer
#         import os
        
#         secret_key = os.getenv("UNSUBSCRIBE_SECRET")
#         if not secret_key:
#             raise ValueError("UNSUBSCRIBE_SECRET not found")
        
#         # Create serializer WITHOUT TimestampSigner (no expiration)
#         serializer = URLSafeSerializer(secret_key) 
        
#         # Generate token
#         token = serializer.dumps(subscription_id)
#         return token
        
#     except Exception as e:
#         logging.error(f"Error generating unsubscribe token: {e}")
#         return None

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
    """Store subscription in Supabase database"""
    try:
        # Prepare data for insertion
        subscription_data = {
            "email": email,
            "journals": journals,  # This should be a list
            "keywords": keywords,
            "frequency": frequency,
            "include_preprints": include_preprints,
            "active": True,
            "last_sent": None
        }
        
        # Insert into database
        response = supabase.table("subscriptions").insert(subscription_data).execute()
        
        if response.data:
            # ✅ RETURN: Consistent format with status
            return {
                "status": "success",
                "data": response.data[0],  # The created subscription
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
        
        # Decode the token
        subscription_id = serializer.loads(token)
        
        # Get subscription from database
        response = supabase.table("subscriptions").select("*").eq("id", subscription_id).single().execute()
        
        if response.data:
            return response.data
        else:
            logging.warning(f"No subscription found for token: {subscription_id}")
            return None
            
    except Exception as e:
        logging.error(f"Error verifying unsubscribe token: {e}")
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