import smtplib
from email.message import EmailMessage
import os
from dotenv import load_dotenv

load_dotenv()

EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")

def send_email(to_email, subject, body):
    try:
        msg = EmailMessage()
        msg['Subject'] = subject
        msg['From'] = EMAIL_ADDRESS
        msg['To'] = to_email
        msg.set_content(body)

        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
            smtp.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
            smtp.send_message(msg)
        print(f"Email sent to {to_email} successfully.")
    
    except Exception as e:
        print(f"Failed to send email: {str(e)}")

if __name__ == "__main__":
    test_email = "rzhan186@gmail.com"  # Replace with your email for testing
    subject = "Test Email"
    body = "This is a test email from the PubMed Tracker."
    send_email(test_email, subject, body)