# install transformer 
# pip3 install transformers

from transformers import pipeline

# Load a pre-trained model
model = pipeline("feature-extraction", model="bert-base-uncased")

def get_keywords(input_text):
    # Preliminary step: Identify the keywords in the user input
    # Here we can simply split based on spaces for extracting terms
    keywords = input_text.lower().split()  # Tokenize by spaces
    keyword_features = model(input_text)  # Get BERT features

    # Generate related terms (this is a simplistic approach; consider enhancing)
    related_terms = []
    
    for keyword in keywords:
        # You can apply token similarity logic to derive related terms based on the vectors
        # Here you might use a database of known synonyms or perform clustering on similar vectors
        related_terms.append(keyword)  # Placeholder for related term logic

    return list(set(related_terms))

# Example use
topic = input("Enter a topic: ")
keywords_query = get_keywords(topic) 
print("Generated keywords:", keywords_query)