
# import functions from other modules
from tracking_main import fetch_pubmed_articles_by_date


def main():
    print("Welcome to the Research Tracker!")
    print("Please select an option:")
    print("1. Track journals")
    print("2. Exit")
    
    while True:
        choice = input("Enter your choice: ")
        
        if choice == '1':
            fetch_pubmed_articles_by_date()
        elif choice == '2':
            print("Exiting...")
            break
        else:
            print("Invalid choice, please try again.")


if __name__ == "__main__":
    main()