# Ratelimit.py

import streamlit as st
import time
from datetime import datetime, timedelta
import random
import os
from Bio import Entrez
import threading

# class PubMedRateLimit:
#     def __init__(self):
#         self.init_session_state()
#         self.configure_entrez()
    
#     def init_session_state(self):
#         if 'rate_limit' not in st.session_state:
#             st.session_state.rate_limit = {
#                 'last_request': None,
#                 'request_count': 0,
#                 'request_times': [],
#                 'blocked_until': None
#             }
    
#     def configure_entrez(self):
#         """Configure Entrez with app email and API key"""
#         # Try to get from Streamlit secrets first, then environment
#         try:
#             email = st.secrets.get("PUBMED_EMAIL", "journaltracker.app@gmail.com")
#         except:
#             email = os.getenv('PUBMED_EMAIL', 'journaltracker.app@gmail.com')
        
#         Entrez.email = email
#         Entrez.tool = "Journal Tracker App"
        
#         # Add API key if available
#         try:
#             api_key = st.secrets.get("NCBI_API_KEY") or os.getenv('NCBI_API_KEY')
#         except:
#             api_key = os.getenv('NCBI_API_KEY')
            
#         if api_key:
#             Entrez.api_key = api_key
#             self.max_requests_per_minute = 15  
#         else:
#             self.max_requests_per_minute = 8
    

#     def can_make_request(self):
#         now = datetime.now()
#         state = st.session_state.rate_limit
        
#         # Check if currently blocked
#         if state['blocked_until'] and now < state['blocked_until']:
#             remaining = (state['blocked_until'] - now).total_seconds()
#             st.error(f"🚫 Too many requests. Please wait {remaining:.0f} seconds before searching again.")
#             return False
        
#         # Clean old requests (keep only last minute)
#         one_minute_ago = now - timedelta(minutes=1)
#         state['request_times'] = [t for t in state['request_times'] if t > one_minute_ago]
        
#         # Check rate limits
#         if len(state['request_times']) >= self.max_requests_per_minute:
#             st.error("🚫 You're searching too quickly. Please wait a minute before trying again.")
#             state['blocked_until'] = now + timedelta(minutes=1)
#             return False
        
#         # Check minimum interval
#         if state['last_request']:
#             time_since_last = (now - state['last_request']).total_seconds()
#             min_interval = 2  # 2 seconds minimum
#             if time_since_last < min_interval:
#                 st.error(f"⏳ Please wait {min_interval - time_since_last:.1f} more seconds.")
#                 return False
        
#         # Update state (silently)
#         state['last_request'] = now
#         state['request_times'].append(now)
#         state['request_count'] += 1
        
#         return True
    
#     def safe_pubmed_search(self, query, max_results=100):
#         """
#         Perform PubMed search with built-in delays and error handling
#         """
#         try:
#             # Add random delay to prevent predictable timing
#             delay = random.uniform(0.8, 1.5)
#             time.sleep(delay)
            
#             # Perform search
#             handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
#             search_results = Entrez.read(handle)
#             handle.close()
            
#             # Small delay after search
#             time.sleep(0.3)
            
#             return search_results
            
#         except Exception as e:
#             error_msg = str(e).lower()
#             if "rate limit" in error_msg or "429" in error_msg:
#                 st.error("🚫 PubMed is busy right now. Please wait a moment and try again.")
#                 # Block for 2 minutes on rate limit error
#                 st.session_state.rate_limit['blocked_until'] = datetime.now() + timedelta(minutes=2)
#                 time.sleep(2)
#             else:
#                 st.error("❌ Search temporarily unavailable. Please try again in a moment.")
#             return None

#     def safe_pubmed_fetch(self, id_list):
#         """
#         Safely fetch article details with rate limiting
#         """
#         try:
#             # Delay before fetch
#             time.sleep(random.uniform(0.5, 1.0))
            
#             # Fetch articles
#             handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="text")
#             articles = handle.read()
#             handle.close()
            
#             return articles
            
#         except Exception as e:
#             st.error("❌ Unable to fetch article details. Please try again.")
#             return None
        


class PubMedRateLimit:
    def __init__(self):
        self.init_session_state()
        self.configure_entrez()
        self.last_request_time = 0
        self.request_lock = threading.Lock()
        self.request_count = 0
    
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
            self.max_requests_per_minute = 600  # 10 requests/second * 60 seconds
            self.max_requests_per_second = 10
            self.min_interval = 0.1  # 100ms between requests
            self.has_api_key = True
            print("✅ NCBI API key configured successfully")
        else:
            self.max_requests_per_minute = 180  # 3 requests/second * 60 seconds
            self.max_requests_per_second = 3
            self.min_interval = 0.34  # ~333ms between requests
            self.has_api_key = False
            print("⚠️ No NCBI API key found - using standard rate limits")
    
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
            wait_time = 1 if self.has_api_key else 2  # Shorter wait with API key
            st.error(f"🚫 You're searching too quickly. Please wait {wait_time} minute(s) before trying again.")
            state['blocked_until'] = now + timedelta(minutes=wait_time)
            return False
        
        # Check minimum interval - more lenient with API key
        if state['last_request']:
            time_since_last = (now - state['last_request']).total_seconds()
            min_interval = 0.5 if self.has_api_key else 2  # Shorter interval with API key
            if time_since_last < min_interval:
                st.error(f"⏳ Please wait {min_interval - time_since_last:.1f} more seconds.")
                return False
        
        # Update state (silently)
        state['last_request'] = now
        state['request_times'].append(now)
        state['request_count'] += 1
        
        return True
    
    def _wait_for_rate_limit(self):
        """Thread-safe rate limiting for batch operations"""
        with self.request_lock:
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            
            if time_since_last < self.min_interval:
                sleep_time = self.min_interval - time_since_last
                time.sleep(sleep_time)
            
            self.last_request_time = time.time()
            self.request_count += 1
            
            # Log progress every 50 requests
            if self.request_count % 50 == 0:
                print(f"📊 API requests made: {self.request_count}")
    
    def safe_pubmed_search(self, query, max_results=100):
        """
        Perform PubMed search with built-in delays and error handling
        """
        try:
            # Apply rate limiting
            self._wait_for_rate_limit()
            
            # Reduced delay with API key
            delay = random.uniform(0.1, 0.3) if self.has_api_key else random.uniform(0.8, 1.5)
            time.sleep(delay)
            
            # Perform search
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            search_results = Entrez.read(handle)
            handle.close()
            
            # Small delay after search
            post_delay = 0.1 if self.has_api_key else 0.3
            time.sleep(post_delay)
            
            return search_results
            
        except Exception as e:
            error_msg = str(e).lower()
            if "rate limit" in error_msg or "429" in error_msg:
                st.error("🚫 PubMed is busy right now. Please wait a moment and try again.")
                # Shorter block time with API key
                block_time = 1 if self.has_api_key else 2
                st.session_state.rate_limit['blocked_until'] = datetime.now() + timedelta(minutes=block_time)
                time.sleep(block_time * 30)  # Wait 30 seconds per minute of block
            else:
                st.error("❌ Search temporarily unavailable. Please try again in a moment.")
            return None

    def safe_pubmed_fetch(self, id_list):
        """
        Safely fetch article details with rate limiting
        """
        try:
            # Apply rate limiting
            self._wait_for_rate_limit()
            
            # Optimized delay with API key
            delay = random.uniform(0.05, 0.15) if self.has_api_key else random.uniform(0.5, 1.0)
            time.sleep(delay)
            
            # Fetch articles
            handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="text")
            articles = handle.read()
            handle.close()
            
            return articles
            
        except Exception as e:
            st.error("❌ Unable to fetch article details. Please try again.")
            return None
    
    def batch_fetch_articles(self, pmid_list, batch_size=None):
        """
        Fetch articles in batches with optimized settings based on API key availability
        """
        if not pmid_list:
            return []
        
        # Optimize batch size based on API key
        if batch_size is None:
            batch_size = 200 if self.has_api_key else 50
        
        # Split into batches
        batches = [pmid_list[i:i + batch_size] for i in range(0, len(pmid_list), batch_size)]
        all_articles = []
        
        for i, batch in enumerate(batches):
            try:
                # Apply rate limiting
                self._wait_for_rate_limit()
                
                # Convert batch to comma-separated string
                pmid_str = ",".join(batch)
                
                # Fetch batch
                handle = Entrez.efetch(db="pubmed", id=pmid_str, rettype="xml")
                articles = Entrez.read(handle)
                handle.close()
                
                # Process articles
                for article in articles.get("PubmedArticle", []):
                    all_articles.append(article)
                
                # Progress indicator
                if len(batches) > 1:
                    progress = (i + 1) / len(batches) * 100
                    print(f"📊 Batch progress: {progress:.1f}% ({i + 1}/{len(batches)})")
                
            except Exception as e:
                print(f"⚠️ Error processing batch {i + 1}: {e}")
                continue
        
        return all_articles
    
    def get_rate_limit_info(self):
        """Get current rate limit configuration info"""
        return {
            "has_api_key": self.has_api_key,
            "max_requests_per_second": self.max_requests_per_second,
            "max_requests_per_minute": self.max_requests_per_minute,
            "min_interval": self.min_interval,
            "total_requests": self.request_count
        }