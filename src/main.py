# import functions from other modules
from tracking_main import fetch_pubmed_articles_by_date
from tracking_main import export_fetched_articles_as_csv

def main():
    print("Welcome to the Research Tracker!")
    print("Please select an option:")
    print("1. Track journals")
    print("2. Exit")
    
    while True:
        choice = input("Enter your choice: ")
        
        if choice == '1':
            
            # Prompt user for email address 
            # email = input("Enter your email address: ").strip()
            
            # Prompt user for journal name
            journal = input("Enter the journal name: ").strip()

            # Prompt user for start date
            start_date = input("Enter the start date (YYYY-MM or YYYY-MM-DD): ").strip()

            # Prompt user for end date
            end_date = input("Enter the end date (YYYY-MM or YYYY-MM-DD): ").strip()

            try:
                # Fetch articles with provided inputs
                articles = fetch_pubmed_articles_by_date(journal, start_date, end_date)
                
                if articles:
                    print(f"About to save {len(articles)} articles as csv format.")
                    # Export articles to CSV
                    export_fetched_articles_as_csv(articles, journal, start_date, end_date)

                else:
                    print("❌ No articles found for the given criteria.")
                
            except ValueError as ve:
                print(f"Error: {ve}")
            except Exception as e:
                print(f"An error occurred: {e}")

        elif choice == '2':
            print("Exiting...")
            break
        else:
            print("Invalid choice, please try again.")


# Running the function with user input for specific dates
if __name__ == "__main__":
    main()