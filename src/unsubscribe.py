# unsubscribe.py – Streamlit page for unsubscribing
import streamlit as st
from supabase import create_client
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

st.set_page_config(page_title="Unsubscribe | PubMed Tracker")

st.title("📭 Unsubscribe from Updates")

# Get token from URL query parameter
token = st.query_params.get("token", None)

if not token:
    st.error("❌ No unsubscribe token provided.")
    st.stop()

# Look up the subscription
response = supabase.table("subscriptions").select("*").eq("unsubscribe_token", token).limit(1).execute()

if not response.data:
    st.error("❌ Invalid or expired unsubscribe token.")
    st.stop()

subscription = response.data[0]

if not subscription.get("active", True):
    st.info("📪 This subscription is already inactive.")
    st.stop()

# Show subscription details
st.markdown("### You are about to unsubscribe from:")
st.markdown(f"""
- **Email**: `{subscription['email']}`
- **Journals**: `{', '.join(subscription['journals'])}`
- **Keywords**: `{subscription['keywords'] or 'None'}`
- **Frequency**: `{subscription['frequency']}`
""")

if st.button("🔕 Confirm Unsubscribe"):
    supabase.table("subscriptions").update({"active": False}).eq("unsubscribe_token", token).execute()
    st.success("✅ You have been unsubscribed from this update.")
