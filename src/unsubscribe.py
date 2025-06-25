# unsubscribe.py

import streamlit as st
from supabase import create_client
import os
from dotenv import load_dotenv
from store_subscription import verify_unsubscribe_token
import logging

# Load environment variables
load_dotenv()
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

# Define the function to handle the unsubscribe process
def handle_unsubscribe(token):
    """Handles the unsubscription process using the provided token."""
    st.set_page_config(page_title="Unsubscribe | PubMed Tracker")
    st.title("🛑 Unsubscribe from PubMed Alerts")

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
        response = supabase.table("subscriptions").update({"active": False}).eq("unsubscribe_token", token).execute()
        
        if response.error:
            st.error("❌ Failed to unsubscribe. Please try again later.")
        else:
            st.success("✅ You have been unsubscribed from this update.")
    else:
        st.info("Click the button above to complete the unsubscribe process.")