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
    
    while True:
        print("\n📘 How it works:")
        print("1. Enter a journal name (e.g., Environmental Science & Technology)")
        print("2. Enter a date range (e.g., 2024-01 to 2024-03)")
        print("3. Optionally enter keyword logic for targeted searches")
        print("   💡 Example: (cadmium OR \"cadmium exposure\") AND rice")
        print("   ❗ Wildcards (*, ?) are NOT currently supported.")
        print("   ⏪ Type 'exit' or 'q' anytime to cancel and return.\n")
        
        print("Please select an option:")
        print("1. Track journals")
        print("2. Exit")

        choice = input("Enter your choice: ").strip().lower()
        
        if choice == '2' or choice in {'exit', 'q'}:
            print("👋 Exiting. Goodbye!")
            break

        if choice != '1':
            print("Invalid choice. Please try again.")
            continue
        
        # Journal name input
        while True:
            journal = input("Enter the journal name: ").strip()
            if journal.lower() in {'exit', 'q'}:
                break

            journal_dict = load_pubmed_journal_abbreviations()
            try:
                formatted_journal = format_journal_abbreviation(journal, journal_dict)
                break
            except ValueError as e:
                print(f"Error: {e}")
                print("Please enter a valid journal name or abbreviation.")

        if journal.lower() in {'exit', 'q'}:
            continue

        # Start date input
        while True:
            start_date = input("Enter the start date (YYYY-MM or YYYY-MM-DD):\n(Press Enter for current month): ").strip()
            if start_date.lower() in {'exit', 'q'}:
                break
            if not start_date:
                print("📅 No start date provided. Using current month.")
                break
            if validate_date_input(start_date):
                break
            print("❌ Invalid format. Please use YYYY-MM or YYYY-MM-DD.")

        if start_date.lower() in {'exit', 'q'}:
            continue

        # End date input
        while True:
            end_date = input("Enter the end date (YYYY-MM or YYYY-MM-DD):\n(Press Enter to use same as start date): ").strip()
            if end_date.lower() in {'exit', 'q'}:
                break
            if not end_date:
                end_date = start_date
                print(f"📅 No end date provided. Using start date as end date: {end_date}")
                break
            if validate_date_input(end_date):
                break
            print("❌ Invalid format. Please use YYYY-MM or YYYY-MM-DD.")

        if end_date.lower() in {'exit', 'q'}:
            continue

        # Keyword loop
        while True:
            raw_keywords = input(
                "\nEnter keyword logic for Title/Abstract search:\n"
                "(e.g., (cadmium OR \"cadmium exposure\") AND rice)\n"
                "Press Enter to skip. Wildcards are NOT supported.\n> "
            ).strip()

            if raw_keywords.lower() in {'exit', 'q'}:
                break

            formatted_keywords = None
            if raw_keywords:
                # Check for unbalanced parentheses
                if raw_keywords.count("(") != raw_keywords.count(")"):
                    print("⚠️ Unbalanced parentheses in keyword logic. Please correct it.")
                    continue
                try:
                    formatted_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
                    print(f"\n🔍 Formatted keyword logic:\n{formatted_keywords}")
                except Exception as e:
                    print(f"❌ Error formatting keywords: {e}")
                    continue

            # Show summary before fetching
            print("\n🔎 Search Summary:")
            print(f"📚 Journal: {journal}")
            print(f"📅 Date Range: {start_date} to {end_date}")
            print(f"🔍 Keywords: {formatted_keywords if formatted_keywords else '[None entered]'}")
            confirm = input("Proceed with this search? (y/n): ").strip().lower()
            if confirm != 'y':
                print("⏪ Returning to main menu...\n")
                break

            try:
                articles = fetch_pubmed_articles_by_date(
                    formatted_journal, start_date, end_date, formatted_keywords
                )

                if articles:
                    print(f"✅ Found {len(articles)} articles. Saving to CSV...")
                    export_fetched_articles_as_csv(articles, journal, start_date, end_date)
                    break
                else:
                    print("❌ No articles found. You can try refining your keywords.")
                    retry = input("🔁 Enter new keywords? (y/n): ").strip().lower()
                    if retry != 'y':
                        break

            except Exception as e:
                print(f"❌ An error occurred: {e}")
                break

if __name__ == "__main__":
    main()
