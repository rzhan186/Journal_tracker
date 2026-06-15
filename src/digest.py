import os
import time
from datetime import datetime, timedelta

import google.generativeai as genai
import streamlit as st
from Bio import Entrez

# Journal display names mapped to their PubMed-searchable names
DIGEST_JOURNALS = {
    "Nature": "Nature",
    "Science": "Science",
    "Cell": "Cell",
    "PNAS": "Proc Natl Acad Sci U S A",
}

PAPERS_PER_JOURNAL = 5


def _fetch_journal_papers(pubmed_name: str, start_date: str, end_date: str) -> list[dict]:
    """Fetch recent papers from a single journal."""
    query = (
        f'"{pubmed_name}"[Journal]'
        f' AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'
        f' AND hasabstract[text]'
        f' AND ("journal article"[Publication Type] OR "review"[Publication Type])'
        f' NOT ("editorial"[Publication Type] OR "letter"[Publication Type] OR "comment"[Publication Type])'
    )

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=PAPERS_PER_JOURNAL, sort="pub date")
        record = Entrez.read(handle)
        handle.close()
        pmid_list = record.get("IdList", [])
    except Exception as e:
        print(f"⚠️ Error searching {pubmed_name}: {e}")
        return []

    papers = []
    for pmid in pmid_list:
        try:
            time.sleep(0.35)
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            info = Entrez.read(handle)
            handle.close()

            article_data = info["PubmedArticle"][0]["MedlineCitation"]["Article"]
            title = str(article_data["ArticleTitle"])

            abstract_sections = article_data.get("Abstract", {}).get("AbstractText", [])
            if not abstract_sections:
                continue
            if isinstance(abstract_sections[0], dict):
                abstract = " ".join(s.get("text", "") for s in abstract_sections)
            else:
                abstract = str(abstract_sections[0])

            doi = None
            for id_item in info["PubmedArticle"][0]["PubmedData"]["ArticleIdList"]:
                if id_item.attributes["IdType"] == "doi":
                    doi = str(id_item)
                    break

            papers.append({
                "title": title,
                "abstract": abstract,
                "doi": f"https://doi.org/{doi}" if doi else None,
            })

        except Exception as e:
            print(f"⚠️ Error fetching PMID {pmid}: {e}")
            continue

    return papers


def _summarize_abstract(title: str, abstract: str) -> str:
    """Generate a plain-English 1-2 sentence summary using Gemini Flash."""
    try:
        try:
            api_key = st.secrets.get("GOOGLE_API_KEY")
        except Exception:
            api_key = os.getenv("GOOGLE_API_KEY")

        if not api_key:
            return abstract[:250] + "..." if len(abstract) > 250 else abstract

        genai.configure(api_key=api_key)
        model = genai.GenerativeModel("gemini-1.5-flash")
        prompt = (
            "Summarize this scientific paper in 1-2 plain English sentences "
            "for a weekly research digest. Focus on the key finding. "
            "Do not begin with 'This paper', 'The study', or 'Researchers'.\n\n"
            f"Title: {title}\n\nAbstract: {abstract}"
        )
        response = model.generate_content(prompt)
        return response.text.strip()

    except Exception as e:
        print(f"⚠️ Summarization failed: {e}")
        return abstract[:250] + "..." if len(abstract) > 250 else abstract


@st.cache_data(ttl=604800, show_spinner=False)
def get_weekly_digest() -> dict:
    """
    Fetch and summarize this week's top papers from Nature, Science, Cell, and PNAS.
    Results are cached for 7 days.

    Returns:
        dict: {journal_display_name: [{"title", "summary", "doi"}, ...]}
    """
    end_date = datetime.today().strftime("%Y/%m/%d")
    start_date = (datetime.today() - timedelta(days=7)).strftime("%Y/%m/%d")

    digest = {}

    for display_name, pubmed_name in DIGEST_JOURNALS.items():
        papers = _fetch_journal_papers(pubmed_name, start_date, end_date)
        summarized = []
        for paper in papers:
            summary = _summarize_abstract(paper["title"], paper["abstract"])
            summarized.append({
                "title": paper["title"],
                "summary": summary,
                "doi": paper["doi"],
            })
        digest[display_name] = summarized

    return digest


def render_weekly_digest():
    """Render the weekly digest in the Streamlit app."""
    today = datetime.today()
    week_start = (today - timedelta(days=7)).strftime("%b %d")
    week_end = today.strftime("%b %d, %Y")

    st.markdown("## 📰 Weekly Science Digest")
    st.caption(f"Top papers from the past week — {week_start} to {week_end}")
    st.markdown("---")

    with st.spinner("Loading this week's digest..."):
        digest = get_weekly_digest()

    if not any(digest.values()):
        st.info("No papers found for this week. Try again later.")
        return

    for journal, papers in digest.items():
        if not papers:
            continue

        st.markdown(f"### {journal}")

        for paper in papers:
            title = paper["title"]
            summary = paper["summary"]
            doi = paper["doi"]

            if doi:
                st.markdown(f"- **[{title}]({doi})** — {summary}")
            else:
                st.markdown(f"- **{title}** — {summary}")

        st.markdown("")

    st.caption("Summaries generated by Gemini Flash · Refreshed weekly · Powered by PubMed")
