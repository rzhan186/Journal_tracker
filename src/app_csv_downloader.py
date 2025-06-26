import streamlit as st
import datetime 
from email_dispatcher import get_csv_from_token  # Import the function

def handle_download():
    if 'token' in st.query_params:
        token = st.query_params['token']
        csv_data, email = get_csv_from_token(token)
        
        if csv_data:
            st.success("✅ Download ready!")
            st.download_button(
                label="📥 Download Your Search Results",
                data=csv_data,
                file_name=f"journal_search_results_{datetime.now().strftime('%Y%m%d')}.csv",
                mime="text/csv"
            )
            st.info("This download link has been used and will expire soon for security.")
        else:
            st.error("❌ Download link has expired or is invalid. Links expire after 24 hours.")
            st.info("Please run a new search or contact support if you need assistance.")
    else:
        st.error("❌ No download token provided.")