import streamlit as st
import time
from datetime import datetime, timedelta
import random
import os
from Bio import Entrez

class PubMedRateLimit:
    def __init__(self):
        self.init_session_state()
        self.configure_entrez()
    
    def init_session_state(self):
        if 'rate_limit' not in st.session_state:
            st.session_state.rate_limit = {
                'last_request': None,
                'request_count': 0,
                'request_times': [],
                'blocked_until': None
            }
    
    def configure_entrez(self):
        """Configure Entrez with app email and API key"""
        # Try to get from Streamlit secrets first, then environment
        try:
            email = st.secrets.get("PUBMED_EMAIL", "journaltracker.app@gmail.com")
        except:
            email = os.getenv('PUBMED_EMAIL', 'journaltracker.app@gmail.com')
        
        Entrez.email = email
        Entrez.tool = "Journal Tracker App"
        
        # Add API key if available
        try:
            api_key = st.secrets.get("NCBI_API_KEY") or os.getenv('NCBI_API_KEY')
        except:
            api_key = os.getenv('NCBI_API_KEY')
            
        if api_key:
            Entrez.api_key = api_key
            self.max_requests_per_minute = 15  
        else:
            self.max_requests_per_minute = 8
    

    def can_make_request(self):
        now = datetime.now()
        state = st.session_state.rate_limit
        
        # Check if currently blocked
        if state['blocked_until'] and now < state['blocked_until']:
            remaining = (state['blocked_until'] - now).total_seconds()
            st.error(f"🚫 Too many requests. Please wait {remaining:.0f} seconds before searching again.")
            return False
        
        # Clean old requests (keep only last minute)
        one_minute_ago = now - timedelta(minutes=1)
        state['request_times'] = [t for t in state['request_times'] if t > one_minute_ago]
        
        # Check rate limits
        if len(state['request_times']) >= self.max_requests_per_minute:
            st.error("🚫 You're searching too quickly. Please wait a minute before trying again.")
            state['blocked_until'] = now + timedelta(minutes=1)
            return False
        
        # Check minimum interval
        if state['last_request']:
            time_since_last = (now - state['last_request']).total_seconds()
            min_interval = 2  # 2 seconds minimum
            if time_since_last < min_interval:
                st.error(f"⏳ Please wait {min_interval - time_since_last:.1f} more seconds.")
                return False
        
        # Update state (silently)
        state['last_request'] = now
        state['request_times'].append(now)
        state['request_count'] += 1
        
        return True
    
    def safe_pubmed_search(self, query, max_results=100):
        """
        Perform PubMed search with built-in delays and error handling
        """
        try:
            # Add random delay to prevent predictable timing
            delay = random.uniform(0.8, 1.5)
            time.sleep(delay)
            
            # Perform search
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            search_results = Entrez.read(handle)
            handle.close()
            
            # Small delay after search
            time.sleep(0.3)
            
            return search_results
            
        except Exception as e:
            error_msg = str(e).lower()
            if "rate limit" in error_msg or "429" in error_msg:
                st.error("🚫 PubMed is busy right now. Please wait a moment and try again.")
                # Block for 2 minutes on rate limit error
                st.session_state.rate_limit['blocked_until'] = datetime.now() + timedelta(minutes=2)
                time.sleep(2)
            else:
                st.error("❌ Search temporarily unavailable. Please try again in a moment.")
            return None

    def safe_pubmed_fetch(self, id_list):
        """
        Safely fetch article details with rate limiting
        """
        try:
            # Delay before fetch
            time.sleep(random.uniform(0.5, 1.0))
            
            # Fetch articles
            handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="text")
            articles = handle.read()
            handle.close()
            
            return articles
            
        except Exception as e:
            st.error("❌ Unable to fetch article details. Please try again.")
            return None
        
    # def show_usage_stats(self):
    #     """Display current usage statistics"""
    #     state = st.session_state.rate_limit
    #     recent_requests = len([t for t in state['request_times'] 
    #                          if t > datetime.now() - timedelta(minutes=1)])
        
    #     # Check if API key is configured
    #     has_api_key = hasattr(Entrez, 'api_key') and Entrez.api_key
    #     api_status = "✅ API Key Active" if has_api_key else "⚠️ No API Key"
        
    #     st.sidebar.info(f"""
    #     **🔍 PubMed API Status:**
    #     {api_status}
        
    #     **📊 Usage Stats:**
    #     - Total searches: {state['request_count']}
    #     - Last minute: {recent_requests}/{self.max_requests_per_minute}
    #     - Rate limit: {self.max_requests_per_minute} requests/minute
        
    #     **⚡ Performance:**
    #     - Email: {Entrez.email}
    #     - Tool: {Entrez.tool}
    #     """)
        
    #     # Show rate limit warning if getting close
    #     if recent_requests >= (self.max_requests_per_minute * 0.8):  # 80% of limit
    #         st.sidebar.warning("⚠️ Approaching rate limit!")