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

    user_email = input("📧 Enter your email address (optional, press Enter to skip): ").strip()

    while True:
        print("\n📘 How it works:")
        print("1. Enter one or more journal names (e.g., Environmental Science & Technology, J Hazard Mater)")
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

        # Journal names input (comma-separated)
        while True:
            journal_input = input("Enter journal name(s), separated by commas: ").strip()
            if journal_input.lower() in {'exit', 'q'}:
                break
            journal_names = [j.strip() for j in journal_input.split(',') if j.strip()]
            if not journal_names:
                print("❌ No valid journal names entered. Try again.")
                continue

            journal_dict = load_pubmed_journal_abbreviations()
            formatted_journals = {}
            invalid_journals = []
            for journal in journal_names:
                try:
                    formatted_journal = format_journal_abbreviation(journal, journal_dict)
                    formatted_journals[journal] = formatted_journal
                except ValueError as e:
                    print(f"❌ Error with '{journal}': {e}")
                    invalid_journals.append(journal)

            if invalid_journals:
                print("\n⚠️ The following journal(s) were not recognized:")
                for j in invalid_journals:
                    print(f"- {j}")
                retry = input("Would you like to re-enter journal names? (y/n): ").strip().lower()
                if retry == 'y':
                    continue
                elif not formatted_journals:
                    print("❌ No valid journals to continue. Returning to main menu.")
                    break

            if formatted_journals:
                break

        if journal_input.lower() in {'exit', 'q'}:
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
                if raw_keywords.count("(") != raw_keywords.count(")"):
                    print("⚠️ Unbalanced parentheses in keyword logic. Please correct it.")
                    continue
                try:
                    formatted_keywords = format_boolean_keywords_for_pubmed(raw_keywords)
                    print(f"\n🔍 Formatted keyword logic:\n{formatted_keywords}")
                except Exception as e:
                    print(f"❌ Error formatting keywords: {e}")
                    continue

            print("\n🔎 Search Summary:")
            print(f"📚 Journals: {', '.join(formatted_journals.keys())}")
            print(f"📅 Date Range: {start_date} to {end_date}")
            print(f"🔍 Keywords: {formatted_keywords if formatted_keywords else '[None entered]'}")
            confirm = input("Proceed with this search? (y/n): ").strip().lower()
            if confirm != 'y':
                print("⏪ Returning to main menu...\n")
                break

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            for user_journal, formatted_journal in formatted_journals.items():
                try:
                    print(f"\n🔄 Searching: {user_journal} ({formatted_journal})")
                    articles = fetch_pubmed_articles_by_date(
                        formatted_journal, start_date, end_date, formatted_keywords
                    )

                    if articles:
                        print(f"✅ Found {len(articles)} articles. Saving to CSV...")
                        export_fetched_articles_as_csv(articles, user_journal, start_date, end_date, timestamp)
                    else:
                        print(f"❌ No articles found for {user_journal}. You can try refining your keywords.")

                except Exception as e:
                    print(f"❌ Error while processing {user_journal}: {e}")

            break  # Exit keyword loop after processing all journals

if __name__ == "__main__":
    main()
