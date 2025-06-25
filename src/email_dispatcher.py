# email_dispatcher.py
from supabase import create_client
from dotenv import load_dotenv
import os
import smtplib
from email.message import EmailMessage
import time

# Load secrets

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")  
EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS") # switch two brevo later. 
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")

if not EMAIL_ADDRESS or not EMAIL_PASSWORD:
    raise ValueError("❌ EMAIL_ADDRESS or EMAIL_PASSWORD not found in .env file")

supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

def send_email(to_email, subject, body):
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = to_email
    msg.set_content(body)

    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
        smtp.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
        smtp.send_message(msg)

def process_subscriptions():
    print("📨 Checking subscriptions...")
    res = supabase.table("subscriptions").select("*").execute()
    for user in res.data:
        email = user['email']
        journals = user['journals']
        frequency = user['frequency']
        body = f"""Hi {email},

This is your {frequency} update for the following journals: {journals}.

🔎 This is just a test. The real search and article results would be inserted here later.

– PubMed Tracker Team
        """
        send_email(email, f"Your {frequency} PubMed update", body)
        print(f"✅ Sent test email to {email}")

if __name__ == "__main__":
    process_subscriptions()
