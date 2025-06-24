from supabase import create_client, Client
from dotenv import load_dotenv
from email_dispatcher import send_email
import os
import datetime

# Load variables from .env file
load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# Supabase project credentials
supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)

def store_user_subscription(email, journals, keywords, start_date, end_date, frequency, csv_bytes):
    # 1. Store metadata in Supabase DB
    data = {
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "start_date": start_date,
        "end_date": end_date,
        "frequency": frequency
    }
    supabase.table("subscriptions").insert(data).execute()

    # 2. Upload CSV file to Supabase Storage
    safe_time = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    filename = f"{email}_{safe_time}.csv"
    supabase.storage.from_("subscription-files").upload(filename, csv_bytes, {"content-type": "text/csv"})

    # 3. Get public URL or signed URL
    signed_url_data = supabase.storage.from_("subscription-files").create_signed_url(filename, 86400)
    signed_url = signed_url_data.get("signedURL")   

    # 4. Send confirmation email
    send_email(
        to=email,
        subject="✅ Subscription Confirmed",
        body=f"""Thank you for subscribing to PubMed updates!

Journals: {', '.join(journals)}
Keywords: {keywords}
Frequency: {frequency}

📥 Download your search results: {signed_url}
""")

    return {"status": "success", "url": signed_url}
