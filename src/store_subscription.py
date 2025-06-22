from supabase import create_client, Client
from dotenv import load_dotenv
import os

# Load variables from .env file
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# Supabase project credentials
supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

def store_user_subscription(email, journals, keywords, start_date, end_date, frequency):
    try:
        # Convert lists to comma-separated strings if needed
        if isinstance(journals, list):
            journals = ", ".join(journals)
        if keywords is None:
            keywords = ""

        data = {
            "email": email,
            "journals": journals,
            "keywords": keywords,
            "start_date": start_date,
            "end_date": end_date,
            "frequency": frequency
        }

        response = supabase.table("subscriptions").insert(data).execute()
        if response.status_code == 201:
            print("✅ Subscription successfully stored.")
        else:
            print(f"⚠️ Failed to store subscription: {response.status_code}")
    except Exception as e:
        print(f"❌ Error storing subscription: {e}")
