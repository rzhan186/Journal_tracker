import streamlit as st
from supabase import create_client
import os
from dotenv import load_dotenv
from store_subscription import verify_unsubscribe_token

# Load environment variables
load_dotenv()
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

st.set_page_config(page_title="Unsubscribe | PubMed Tracker")
st.title("🛑 Unsubscribe from PubMed Alerts")

# Get token from URL query parameter
token = st.query_params.get("token", None)

if not token:
    st.error("❌ No unsubscribe token provided.")
    st.stop()

# Verify the unsubscribe token
subscription = verify_unsubscribe_token(token)

if not subscription:
    st.error("❌ Invalid or expired unsubscribe link.")
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
    # Update the subscription status in the database to inactive
    supabase.table("subscriptions").update({"active": False}).eq("unsubscribe_token", token).execute()
    st.success("✅ You have been unsubscribed from this update.")
else:
    st.info("Click the button above to complete the unsubscribe process.")