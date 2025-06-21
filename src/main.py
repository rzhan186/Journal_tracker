# import functions from other modules
from datetime import datetime

from tracking_main import (
    fetch_pubmed_articles_by_date,
    export_fetched_articles_as_csv,
    validate_date_input,
    load_pubmed_journal_abbreviations,
    format_journal_abbreviation,
    format_boolean_keywords_for_pubmed,
)

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

                if not start_date:
                    print(f"No start date provided. Using current month.")
                    break
                
                if validate_date_input(start_date):
                    break
                else:
                    print("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15). Please try again.")

            while True:
                end_date = input("Enter the end 🏁 date (YYYY-MM or YYYY-MM-DD):\n(Press enter to search for the current month): ").strip()

                if not end_date:
                    end_date = start_date
                    print(f"No end date provided. Using start date as end date: {end_date}.")
                    break

                if validate_date_input(end_date):
                    break
                else:
                    print("❌ Invalid format for end date. Please enter again.")

            # 🧠 Loop to allow retrying keywords if no articles are found
            while True:
                raw_keywords = input(
                    "\nEnter keyword logic for Title/Abstract search\n"
                    "(e.g., (climat* OR \"global warming\") AND (mercury OR pollution)):\n"
                    "(Press Enter to skip keyword filtering): "
                ).strip()

                formatted_keywords = None
                if raw_keywords:
                    try:
                        formatted_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
                        print(f"\n🔍 Using keyword filter:\n{formatted_keywords}\n")
                    except Exception as e:
                        print(f"⚠️ Failed to format keyword logic: {e}")
                        continue  # Retry keyword input

                # Try fetching articles
                try:
                    articles = fetch_pubmed_articles_by_date(
                        formatted_journal, start_date, end_date, formatted_keywords
                    )

                    if articles:
                        print(f"✅ Found {len(articles)} articles. Saving to CSV...")
                        export_fetched_articles_as_csv(articles, journal, start_date, end_date)
                        break  # ✅ Exit keyword loop on success
                    else:
                        print("❌ No articles found. You can try refining your keywords.")
                        retry = input("Do you want to enter new keywords? (y/n): ").strip().lower()
                        if retry != 'y':
                            break  # Exit loop even if no articles
                except Exception as e:
                    print(f"❌ An error occurred while fetching articles: {e}")
                    break  # Don't retry on serious errors

        elif choice == '2':
            print("Exiting...")
            break
        else:
            print("Invalid choice, please try again.")

if __name__ == "__main__":
    main()
