# functions to test if key functionalities are working

def test_token_system():
    """Test the token generation and verification system"""
    try:
        # Get a real subscription ID from the database
        SUPABASE_URL = os.getenv("SUPABASE_URL")
        SUPABASE_KEY = os.getenv("SUPABASE_KEY")
        supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
        
        # Get the first active subscription
        response = supabase.table("subscriptions").select("id, email").eq("active", True).limit(1).execute()
        
        if not response.data:
            print("No active subscriptions found in database")
            return False
            
        test_subscription_id = response.data[0]['id']  # This will be a UUID
        test_email = response.data[0]['email']
        
        print(f"Testing with subscription ID: {test_subscription_id}")
        print(f"Testing with email: {test_email}")
        
        # Generate token
        token = generate_unsubscribe_token(test_subscription_id)
        print(f"Generated token: {token}")
        
        # Verify token
        subscription = verify_unsubscribe_token(token)
        print(f"Verified subscription: {subscription is not None}")
        
        if subscription:
            print(f"Email: {subscription.get('email')}")
            print(f"ID matches: {subscription.get('id') == test_subscription_id}")
            return True
        else:
            print("Token verification failed")
            return False
            
    except Exception as e:
        print(f"Test failed with error: {e}")
        return False
