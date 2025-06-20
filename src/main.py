# import functions from other modules
from datetime import datetime

from tracking_main import fetch_pubmed_articles_by_date
from tracking_main import export_fetched_articles_as_csv
from tracking_main import validate_date_input
from tracking_main import load_pubmed_journal_abbreviations
from tracking_main import format_journal_abbreviation

# def main():
#     print("Welcome to the Research Tracker!")
#     print("Please select an option:")
#     print("1. Track journals")
#     print("2. Exit")
    
#     while True:
#         choice = input("Enter your choice: ")
        
#         if choice == '1':
            
#             while True: 
#             # Prompt user for email address 
#             # email = input("Enter your email address: ").strip()
            
#                 # Prompt user for journal name
#                 journal = input("Enter the journal name: ").strip()
#                 journal_dict = load_pubmed_journal_abbreviations()

#                 try: 
#                     formatted_journal = format_journal_abbreviation(journal, journal_dict)
#                     break
#                 except ValueError as e:
#                     print(f"Error: {e}")
#                     print("Please enter a valid journal name or abbreviation.")
            
#             # Now use formatted_journal in subsequent calls
#             # Prompt user for start date
#             start_date = input("Enter the start date (YYYY-MM or YYYY-MM-DD): \n(Press enter again if like to search the current month.) ").strip()
            
#             # Initialize end_date as None by default
#             end_date = None

#             try: 
#                 # Validate start date format right away if it's provided
#                 if start_date:
#                     if not validate_date_input(start_date):
#                         raise ValueError("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15).")
                        
#                     end_date = input("Enter the end date (YYYY-MM or YYYY-MM-DD): ").strip()

#                     # Fetch articles with provided inputs
#                     articles = fetch_pubmed_articles_by_date(formatted_journal, start_date, end_date)

#                     if articles:
#                         print(f"About to save {len(articles)} articles as csv format.")
#                         # Export articles to CSV
#                         export_fetched_articles_as_csv(articles, journal, start_date, end_date)

#                         # End the program once finished exporting
#                         break

#                     else:
#                         print("❌ No articles found for the given criteria.")
                    
#             except ValueError as ve:
#                 print(f"Error: {ve}")
#             except Exception as e:
#                 print(f"An error occurred: {e}")

#         elif choice == '2':
#             print("Exiting...")
#             break
#         else:
#             print("Invalid choice, please try again.")



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
                start_date = input("Enter the start date (YYYY-MM or YYYY-MM-DD):\n(Press enter to search the current month): ").strip()

                if not start_date:  # Use current month if skipped
                    today = datetime.today()
                    start_date = today.strftime("%Y-%m")  # Set start date to current month
                    print(f"No start date provided. Using current month: {start_date}.")
                    break  # Exit the loop since we have a valid date
                
                # Validate the start date format
                if validate_date_input(start_date):
                    break  # Exit the loop if valid
                else:
                    print("❌ Invalid format! Dates must be in YYYY-MM or YYYY-MM-DD format (e.g., 2024-01 or 2024-01-15). Please try again.")

            # Prompt user for end date
            end_date = input("Enter the end 🏁 date (YYYY-MM or YYYY-MM-DD):\n(Press enter to search the current month): ").strip()
            if not end_date:  # Set end date to start date if not provided
                end_date = start_date

            # Validate end date format
            if end_date and not validate_date_input(end_date):
                print("❌ Invalid format for end date. Please ensure it's in the right format. It will be set to the same as start date.")

            # Fetch articles with provided inputs
            articles = fetch_pubmed_articles_by_date(formatted_journal, start_date, end_date)

            # Informing the user about article fetching results
            if articles:
                print(f"About to save {len(articles)} articles as CSV format.")
                # Export articles to CSV
                export_fetched_articles_as_csv(articles, journal, start_date, end_date)
            else:
                print("❌ No articles found for the given criteria.")

        elif choice == '2':
            print("Exiting...")
            break
        else:
            print("Invalid choice, please try again.")

# Running the function with user input for specific dates
if __name__ == "__main__":
    main()
