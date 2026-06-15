import os
import time
from datetime import datetime, timedelta

from groq import Groq
import streamlit as st
from Bio import Entrez

# Journal display names mapped to their PubMed-searchable names
DIGEST_JOURNALS = {
    "Nature": "Nature",
    "Science": "Science",
    "Cell": "Cell",
}

PAPERS_PER_JOURNAL = 50  # Fetch up to 50; in practice covers all weekly publications


def _fetch_journal_papers(pubmed_name: str, start_date: str, end_date: str) -> list[dict]:
    """Fetch all papers published in a journal during the given date range."""
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
    """Generate a plain-English 3-bullet summary using Llama 3 via Groq."""
    try:
        try:
            api_key = st.secrets.get("GROQ_API_KEY")
        except Exception:
            api_key = os.getenv("GROQ_API_KEY")

        if not api_key:
            return "- Summary unavailable (API key not configured)."

        client = Groq(api_key=api_key)
        prompt = (
            "Summarize this scientific paper in exactly 3 short bullet points for a general audience with no scientific background. "
            "Each bullet should be one plain English sentence. "
            "Cover: (1) what problem or question the study addresses, (2) what they found, (3) why it matters. "
            "Use simple everyday language. No jargon. No markdown — just start each bullet with a dash (-).\n\n"
            f"Title: {title}\n\nAbstract: {abstract}"
        )
        response = client.chat.completions.create(
            messages=[{"role": "user", "content": prompt}],
            model="llama-3.1-8b-instant",
        )
        return response.choices[0].message.content.strip()

    except Exception as e:
        print(f"⚠️ Summarization failed: {e}")
        return f"- Summary unavailable ({str(e)})."


@st.cache_data(ttl=604800, show_spinner=False)
def get_weekly_digest() -> dict:
    """
    Fetch and summarize all papers from the past week in Nature, Science, and Cell.
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
    st.caption(f"All papers from the past week — {week_start} to {week_end}")
    st.markdown("---")

    with st.spinner("Loading this week's digest..."):
        digest = get_weekly_digest()

    if not any(digest.values()):
        st.info("No papers found for this week. Try again later.")
        return

    # Keyword filter
    keyword = st.text_input(
        "🔍 Filter by keyword",
        placeholder="e.g. cancer, CRISPR, climate, neuroscience...",
    )

    st.markdown("")

    keyword_lower = keyword.strip().lower()
    total_shown = 0

    for journal, papers in digest.items():
        if not papers:
            continue

        # Apply keyword filter across title and summary
        if keyword_lower:
            filtered = [
                p for p in papers
                if keyword_lower in p["title"].lower()
                or keyword_lower in p["summary"].lower()
            ]
        else:
            filtered = papers

        if not filtered:
            continue

        st.markdown(f"### {journal} ({len(filtered)} paper{'s' if len(filtered) != 1 else ''})")

        for paper in filtered:
            title = paper["title"]
            summary = paper["summary"]
            doi = paper["doi"]

            if doi:
                st.markdown(f"**[{title}]({doi})**")
            else:
                st.markdown(f"**{title}**")

            st.markdown(summary)
            st.markdown("")

        total_shown += len(filtered)

    if keyword_lower and total_shown == 0:
        st.info(f"No papers found matching **{keyword}**. Try a different keyword.")

    st.caption("Summaries generated by Llama 3 via Groq · Refreshed weekly · Powered by PubMed")
