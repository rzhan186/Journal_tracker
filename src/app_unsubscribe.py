# app_unsubscribe.py

import streamlit as st
from supabase import create_client
import os
from dotenv import load_dotenv
import logging
from store_subscription import get_user_subscriptions_by_email, verify_unsubscribe_token
import pandas as pd

# Load environment variables
load_dotenv()
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

# Initialize Supabase client
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

def handle_unsubscribe(token):
    """Handles the selective unsubscribe process using the provided token."""
    st.set_page_config(page_title="Manage Subscriptions | PubMed Tracker")
    st.title("📧 Manage Your PubMed Subscriptions")

    # Verify the unsubscribe token to get the user's email
    subscription = verify_unsubscribe_token(token)

    # Check if the token is valid
    if not subscription:
        st.error("❌ Invalid or expired unsubscribe link.")
        st.info("Please use the most recent unsubscribe link from your email, or contact support.")
        return

    user_email = subscription['email']
    st.success(f"✅ Welcome back, {user_email}")

    # Get all active subscriptions for this user
    try:
        all_subscriptions = get_user_subscriptions_by_email(user_email)
        
        if not all_subscriptions:
            st.info("📭 You don't have any active subscriptions.")
            st.markdown("---")
            st.info("🔍 **Want to create a new subscription?** [Go to main app](https://journaltracker.streamlit.app)")
            return

        st.markdown(f"### 📋 Your Active Subscriptions ({len(all_subscriptions)} total)")
        st.markdown("Select the subscriptions you want to **unsubscribe from**:")

        # Create a list to store selected subscription IDs
        selected_for_removal = []

        # Display each subscription with a checkbox
        for i, sub in enumerate(all_subscriptions):
            checkbox_key = f"unsub_{sub['id']}"
            
            # Format subscription details - PROPER JSON HANDLING
            journals_raw = sub.get('journals', [])
            
            # Handle JSON data properly
            if isinstance(journals_raw, list):
                journals_display = ', '.join(journals_raw) if journals_raw else 'None'
            else:
                journals_display = 'None'
            
            # FIX THE KEYWORDS DISPLAY
            keywords_display = sub.get('keywords') if sub.get('keywords') else 'None'
            
            # FIX THE PREPRINTS DISPLAY
            include_preprints = sub.get('include_preprints', False)
            preprints_display = 'Yes' if include_preprints else 'No'
            
            # Create columns for better layout
            col1, col2 = st.columns([1, 4])
            
            with col1:
                # Checkbox to select for unsubscription
                selected = st.checkbox(
                    "Remove", 
                    key=checkbox_key,
                    help="Check to unsubscribe from this alert"
                )
            
            with col2:
                # Subscription details in an expandable section
                with st.expander(f"🔔 Subscription #{i+1} - {sub['frequency']}", expanded=False):
                    st.markdown(f"""
                    - **📘 Journals**: {journals_display}
                    - **📑 Preprints**: {preprints_display}
                    - **🔑 Keywords**: {keywords_display}
                    - **🔁 Frequency**: {sub['frequency']}
                    - **📅 Created**: {sub.get('created_at', 'N/A')[:10]}
                    - **🆔 ID**: {sub['id']}
                    """)
                    
                    # DEBUG: Show raw data (remove this after fixing)
                    if st.checkbox(f"🔧 Show raw data", key=f"debug_{sub['id']}"):
                        st.json(sub)
            
            # Add to removal list if selected
            if selected:
                selected_for_removal.append(sub['id'])
            
            st.markdown("---")

        # Show summary of selections
        if selected_for_removal:
            st.warning(f"⚠️ You are about to unsubscribe from **{len(selected_for_removal)}** subscription(s).")
            
            # Confirmation section
            st.markdown("### ⚠️ Confirm Unsubscription")
            
            col1, col2 = st.columns(2)
            
            with col1:
                if st.button("🔕 Confirm Unsubscribe Selected", type="primary"):
                    unsubscribe_selected_subscriptions(selected_for_removal, user_email)
            
            with col2:
                if st.button("❌ Cancel", type="secondary"):
                    st.rerun()  # Refresh the page to clear selections
        
        else:
            st.info("💡 Check the boxes next to subscriptions you want to remove, then click 'Confirm Unsubscribe'.")
        
        # Footer options
        st.markdown("---")
        st.markdown("### 🔧 Other Options")
        
        col1, col2 = st.columns(2)
        with col1:
            if st.button("🔕 Unsubscribe from ALL"):
                unsubscribe_all_subscriptions(user_email)
        
        with col2:
            st.markdown("[🔍 Create New Subscription](https://journaltracker.streamlit.app)")

    except Exception as e:
        st.error("❌ An error occurred while loading your subscriptions.")
        logging.error(f"Error loading subscriptions for {user_email}: {e}")

def unsubscribe_selected_subscriptions(subscription_ids, user_email):
    """Unsubscribe from selected subscriptions"""
    try:
        # Update selected subscriptions to inactive
        for sub_id in subscription_ids:
            response = supabase.table("subscriptions").update({"active": False}) \
                .eq("id", sub_id).eq("email", user_email).execute()
        
        st.success(f"✅ Successfully unsubscribed from {len(subscription_ids)} subscription(s).")
        st.info("You will no longer receive updates for the selected subscriptions.")
        
        # Option to refresh the page
        if st.button("🔄 Refresh Page"):
            st.rerun()
            
    except Exception as e:
        st.error("❌ Failed to unsubscribe. Please try again later.")
        logging.error(f"Error unsubscribing selected subscriptions: {e}")

def unsubscribe_all_subscriptions(user_email):
    """Unsubscribe from all subscriptions for a user"""
    try:
        # Confirm before proceeding
        st.warning("⚠️ This will unsubscribe you from ALL active subscriptions.")
        
        if st.button("🔴 Yes, Unsubscribe from ALL", key="confirm_all"):
            response = supabase.table("subscriptions").update({"active": False}) \
                .eq("email", user_email).eq("active", True).execute()
            
            st.success("✅ You have been unsubscribed from all PubMed alerts.")
            st.info("You will no longer receive any automated updates.")
            
            # Option to create new subscription
            st.markdown("---")
            st.info("🔍 **Want to create a new subscription?** [Go to main app](https://journaltracker.streamlit.app)")
            
    except Exception as e:
        st.error("❌ An error occurred while unsubscribing from all subscriptions.")
        logging.error(f"Error unsubscribing all for {user_email}: {e}")

