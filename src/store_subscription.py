# unsubscribe.py
from itsdangerous import URLSafeSerializer, BadSignature
from supabase import create_client, Client
from dotenv import load_dotenv
import os
import streamlit as st

# Load environment variables
load_dotenv()

# Initialize Supabase client
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SECRET_KEY = os.getenv("UNSUBSCRIBE_SECRET")  # Must be defined in your .env

SECRET_KEY = os.getenv("UNSUBSCRIBE_SECRET")
if not SECRET_KEY:
    raise ValueError("❌ UNSUBSCRIBE_SECRET not set in .env")


supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)
serializer = URLSafeSerializer(SECRET_KEY, salt="unsubscribe")

def generate_unsubscribe_token(email, journals, keywords, frequency):
    data = {
        "email": email,
        "journals": journals,
        "keywords": keywords,
        "frequency": frequency,
    }
    return serializer.dumps(data)

def verify_unsubscribe_token(token):
    try:
        return serializer.loads(token)
    except BadSignature:
        return None

def delete_subscription(subscription):
    response = supabase.table("subscriptions") \
        .delete() \
        .match({
            "email": subscription["email"],
            "journals": subscription["journals"],
            "keywords": subscription["keywords"],
            "frequency": subscription["frequency"]
        }) \
        .execute()
    return response

# Streamlit App Interface
st.set_page_config(page_title="Unsubscribe", layout="centered")
st.title("🛑 Unsubscribe from PubMed Alerts")

token = st.query_params.get("token")

if not token:
    st.error("❌ No unsubscribe token provided.")
    st.stop()

subscription = verify_unsubscribe_token(token)

if not subscription:
    st.error("❌ Invalid or expired unsubscribe link.")
    st.stop()

st.write("You're about to unsubscribe from the following:")
st.markdown(f"""
- **Email**: {subscription['email']}
- **Journals**: {', '.join(subscription['journals'])}
- **Keywords**: {subscription['keywords'] or "None"}
- **Frequency**: {subscription['frequency']}
""")

if st.button("✅ Confirm Unsubscribe"):
    result = delete_subscription(subscription)
    st.success("✅ You have been unsubscribed.")
    st.write("📭 You will no longer receive updates for this configuration.")
else:
    st.info("Click the button above to complete the unsubscribe process.")
