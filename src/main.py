# import functions from other modules
from datetime import datetime

from tracking_main import fetch_pubmed_articles_by_date
from tracking_main import export_fetched_articles_as_csv
from tracking_main import validate_date_input
from tracking_main import load_pubmed_journal_abbreviations
from tracking_main import format_journal_abbreviation

def main():
    print("Welcome to the Research Tracker!")
    print("Please select an option:")
    print("1. Track journals")
    print("2. Exit")
    
    while True:
        choice = input("Enter your choice: ").strip()
        
        if choice == '1':
            while True:
                # Prompt user for journal name
                journal = input("Enter the journal name: ").strip()
                
                journal_dict = load_pubmed_journal_abbreviations()

                try: 
                    formatted_journal = format_journal_abbreviation(journal, journal_dict)
                    break  # Exit inner loop if journal is valid
                except ValueError as e:
                    print(f"Error: {e}")
                    print("Please enter a valid journal name or abbreviation.")

            # Start date validation loop
            while True:
                start_date = input("Enter the start date (YYYY-MM or YYYY-MM-DD):\n(Press enter to search for the current month): ").strip()

                if not start_date:  # Use current month if skipped
                    # today = datetime.today()
                    # start_date = today.strftime("%Y-%m")  # Set start date to current month
                    print(f"No start date provided. Using current month: {start_date}.")
                    break  # Exit the loop since we have a valid date
                
                # Validate the start date format
                if validate_date_input(start_date):
                    break  # Exit the loop if valid
                else:
                    print("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15). Please try again.")

            while True:
                # Prompt user for end date
                end_date = input("Enter the end 🏁 date (YYYY-MM or YYYY-MM-DD):\n(Press enter to search for the current month: ").strip()

                if not end_date:  # Set end date to start date if not provided
                    end_date = start_date
                    print(f"No end date provided. Using start date as end date: {end_date}.")
                    break

                # Validate end date format
                if validate_date_input(end_date):
                    break
                else:
                    print("❌ Invalid format for end date. Please enter again.")

            # Fetch articles with provided inputs and catch potential errors
            try:
                articles = fetch_pubmed_articles_by_date(formatted_journal, start_date, end_date)

                # Informing the user about article fetching results
                if articles:
                    print(f"About to save {len(articles)} articles as CSV format.")
                    # Export articles to CSV
                    export_fetched_articles_as_csv(articles, journal, start_date, end_date)
                else:
                    print("❌ No articles found for the given criteria.")

            except Exception as e:
                print(f"An error occurred while fetching articles: {e}")

        elif choice == '2':
            print("Exiting...")
            break
        else:
            print("Invalid choice, please try again.")

# Running the function with user input for specific dates
if __name__ == "__main__":
    main()

