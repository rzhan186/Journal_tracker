# scheduler.py
import schedule 
import time
from datetime import datetime, timedelta
from store_subscription import supabase
from email_dispatcher import send_email, generate_download_token
# ... other imports

def process_due_subscriptions():
    """Check for subscriptions that are due for updates"""
    try:
        # Get all active subscriptions
        response = supabase.table("subscriptions").select("*").eq("active", True).execute()
        
        for subscription in response.data:
            if is_subscription_due(subscription):
                send_subscription_update(subscription)
                update_next_send_date(subscription)
                
    except Exception as e:
        logging.error(f"Error processing subscriptions: {e}")

# Run every hour
schedule.every().hour.do(process_due_subscriptions)

while True:
    schedule.run_pending()
    time.sleep(60)