# src/main.py

def main():
    print("Welcome to the Research Tracker!")
    print("Please select an option:")
    print("1. Search for articles")
    print("2. Track journals")
    print("3. Exit")
    
    while True:
        choice = input("Enter your choice: ")
        
        if choice == '1':
            search_articles()
        elif choice == '2':
            track_journals()
        elif choice == '3':
            print("Exiting...")
            break
        else:
            print("Invalid choice, please try again.")

def search_articles():
    print("Functionality for searching articles will go here.")

def track_journals():
    print("Functionality for tracking journals will go here.")

if __name__ == "__main__":
    main()