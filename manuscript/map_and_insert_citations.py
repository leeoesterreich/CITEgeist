#!/usr/bin/env python3
"""
Map V3 Zotero citations to V4 references and create insertion script.

Strategy:
1. Extract author/title from V3 Zotero citations
2. Match against V4 bibliography
3. Create a mapping of V4 citation numbers to V3 Zotero field codes
4. Insert the field codes into V4.docx
"""
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}

def load_v3_citations():
    """Load the extracted V3 Zotero citations."""
    with open('zotero_full_extraction.json', 'r') as f:
        data = json.load(f)
    return data['citations']

def extract_v4_bibliography():
    """Extract bibliography from V4 markdown."""
    refs = {}
    with open('CITEgeist_Patterns_v4.md', 'r') as f:
        content = f.read()

    # Find bibliography section
    lines = content.split('\n')
    in_refs = False
    for line in lines:
        if line.strip().startswith('## References'):
            in_refs = True
            continue
        if in_refs:
            # Match numbered references like "1. Author Name..."
            match = re.match(r'^(\d+)\.\s+(.+)', line.strip())
            if match:
                num = match.group(1)
                ref_text = match.group(2)
                refs[num] = ref_text

    return refs

def match_citations(v3_citations, v4_refs):
    """Match V3 citations to V4 references by DOI, title keywords, or author+year."""
    mapping = {}  # v4_num -> v3_citation

    # Build lookup by DOI and title keywords
    for v3_cit in v3_citations:
        if 'json_data' not in v3_cit:
            continue

        items = v3_cit['json_data'].get('citationItems', [])
        for item in items:
            item_data = item.get('itemData', {})

            # Get identifiers
            doi = item_data.get('DOI', '')
            authors = item_data.get('author', [])
            first_author = authors[0].get('family', '') if authors else ''
            title = item_data.get('title', '')

            # Extract key title words (remove common words)
            skip_words = {'the', 'and', 'for', 'with', 'from', 'using', 'based', 'analysis'}
            title_keywords = [w.lower() for w in title.split() if len(w) > 4 and w.lower() not in skip_words][:5]

            # Get year from issued
            year = ''
            issued = item_data.get('issued', {})
            date_parts = issued.get('date-parts', [[]])
            if date_parts and date_parts[0]:
                year = str(date_parts[0][0])

            # Search V4 refs for match
            for v4_num, v4_text in v4_refs.items():
                if v4_num in mapping:
                    continue

                v4_lower = v4_text.lower()

                # Best match: DOI
                if doi and doi.lower() in v4_lower:
                    mapping[v4_num] = v3_cit
                    print(f"  V4[{v4_num}] <- V3[{v3_cit['formatted']}] (DOI): {first_author}")
                    continue

                # Good match: multiple title keywords + author
                keyword_matches = sum(1 for kw in title_keywords if kw in v4_lower)
                if keyword_matches >= 3 and first_author.lower() in v4_lower[:50]:
                    mapping[v4_num] = v3_cit
                    print(f"  V4[{v4_num}] <- V3[{v3_cit['formatted']}] (keywords): {first_author} - {title[:40]}...")
                    continue

                # Fallback: author + year + distinctive keyword
                if first_author.lower() in v4_lower[:40] and year in v4_text:
                    # Require at least one distinctive title word
                    if keyword_matches >= 1:
                        mapping[v4_num] = v3_cit
                        print(f"  V4[{v4_num}] <- V3[{v3_cit['formatted']}] (author+year): {first_author} {year}")

    return mapping

def create_v4_field_codes(mapping):
    """Create new Zotero field codes renumbered for V4."""
    new_citations = {}

    for v4_num, v3_cit in mapping.items():
        # Copy the V3 field code but update the formatted citation number
        old_field = v3_cit['full_field']

        # Update the formatted citation number
        new_field = re.sub(
            r'"formattedCitation":"[^"]*"',
            f'"formattedCitation":"[{v4_num}]"',
            old_field
        )
        new_field = re.sub(
            r'"plainCitation":"[^"]*"',
            f'"plainCitation":"[{v4_num}]"',
            new_field
        )

        new_citations[v4_num] = {
            'field_code': new_field,
            'original_v3': v3_cit['formatted']
        }

    return new_citations

def main():
    print("Loading V3 Zotero citations...")
    v3_citations = load_v3_citations()
    print(f"  {len(v3_citations)} citations loaded")

    print("\nExtracting V4 bibliography...")
    v4_refs = extract_v4_bibliography()
    print(f"  {len(v4_refs)} references found")

    print("\nMatching V3 citations to V4 references...")
    mapping = match_citations(v3_citations, v4_refs)
    print(f"  {len(mapping)} matches found")

    print("\nCreating V4 field codes...")
    new_citations = create_v4_field_codes(mapping)

    # Save the mapping
    output = {
        'v4_citation_count': len(v4_refs),
        'matched_count': len(mapping),
        'citations': new_citations
    }

    with open('v4_zotero_mapping.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nSaved mapping to v4_zotero_mapping.json")

    print("\n=== Summary ===")
    print(f"V4 references: {len(v4_refs)}")
    print(f"Matched to V3 Zotero: {len(mapping)}")
    print(f"Unmatched: {len(v4_refs) - len(mapping)}")

    unmatched = set(v4_refs.keys()) - set(mapping.keys())
    if unmatched:
        print(f"\nUnmatched V4 references: {sorted(unmatched, key=int)}")
        for num in sorted(unmatched, key=int):
            print(f"  [{num}] {v4_refs[num][:60]}...")

if __name__ == '__main__':
    main()
