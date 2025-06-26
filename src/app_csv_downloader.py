# 

import streamlit as st
from datetime import datetime
from email_dispatcher import get_csv_from_token

def handle_download(token):
    """
    Handle CSV download requests with secure token validation.
    Token expires after 24 hours for security.
    """
    
    st.title("📥 Download Your Search Results")
    
    try:
        # Validate token and retrieve CSV data
        csv_data, email = get_csv_from_token(token)
        
        if csv_data and email:
            st.success("✅ Download link validated successfully!")
            st.info(f"📧 Results prepared for: {email}")
            
            # Generate filename with timestamp
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"journal_search_results_{timestamp}.csv"
            
            # Provide download button
            st.download_button(
                label="📥 Download Your Search Results (CSV)",
                data=csv_data,
                file_name=filename,
                mime="text/csv",
                help="Click to download your personalized search results"
            )
            
            # Security notice
            st.warning("🔒 **Security Notice:** This download link will expire soon for your protection. Please save your results now.")
            
            # Show file info
            file_size = len(csv_data) / 1024  # Size in KB
            st.caption(f"📄 File size: {file_size:.1f} KB | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            
            # Option to return to main app
            if st.button("🔙 Return to Main App"):
                st.switch_page("app.py")  # Adjust if your main file has a different name
                
        else:
            # Token invalid or expired
            st.error("❌ **Download Link Invalid or Expired**")
            st.markdown("""
            **Possible reasons:**
            - Link has expired (links expire after 24 hours)
            - Link has already been used
            - Invalid or corrupted link
            
            **What to do:**
            - Run a new search to generate fresh results
            - Check if you have a more recent email with a valid link  
            - Contact support if you continue having issues
            """)
            
            # Provide navigation back to main app
            if st.button("🔍 Go to Search Page"):
                st.switch_page("app.py")  # Adjust if your main file has a different name
                
    except Exception as e:
        # Handle any unexpected errors
        st.error("❌ **Download Error**")
        st.markdown(f"""
        An unexpected error occurred while processing your download:
                    
        Please try again or contact support if the issue persists.
        """)

# Navigation option
if st.button("🏠 Return to Journal Tracker"):
    st.switch_page("app.py")

# Footer with additional info
st.markdown("---")
st.caption("🔐 All downloads are secured with time-limited tokens for your privacy and security.")