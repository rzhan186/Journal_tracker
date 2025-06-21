import nltk
nltk.download('wordnet')
nltk.download('omw-1.4')  # Optional but helps with multilingual synonyms

from nltk.corpus import wordnet as wn

def expand_with_wordnet(term):
    """
    Returns a list of synonyms for a given English word using WordNet.
    """
    synonyms = set()

    for synset in wn.synsets(term):
        for lemma in synset.lemmas():
            word = lemma.name().replace('_', ' ')
            synonyms.add(word.lower())

    # Always include the original term
    synonyms.add(term.lower())

    return list(synonyms)

print(expand_with_wordnet("pollution"))
# Output might include: ['pollution', 'fouling', 'contamination', 'defilement']

print(expand_with_wordnet("climate"))
# Output might include: ['climate', 'clime', 'atmosphere', 'weather']

def generate_keywords_query_with_wordnet(user_input):
    """
    Generate a PubMed query by expanding each word with WordNet synonyms.
    """
    import re

    raw_words = re.findall(r'\b\w+\b', user_input.lower())
    STOPWORDS = {"the", "and", "of", "to", "in", "on", "with", "a", "an", "by", "for", "from", "at"}
    base_keywords = [word for word in raw_words if word not in STOPWORDS]

    all_keywords = set()
    for kw in base_keywords:
        expanded = expand_with_wordnet(kw)
        all_keywords.update(expanded)

    query_parts = [f'"{kw}"[Title/Abstract]' for kw in all_keywords]
    final_query = " OR ".join(query_parts)
    print(f"🔍 Final expanded query:\n{final_query}")
    return final_query

user_input = input("Enter a topic: ")
query = generate_keywords_query_with_wordnet(user_input)
print(query)