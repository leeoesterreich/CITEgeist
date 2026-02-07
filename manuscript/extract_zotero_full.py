#!/usr/bin/env python3
"""
Extract all Zotero field codes from V3 docx with full JSON metadata.
"""
import zipfile
import re
import json
from pathlib import Path
from xml.etree import ElementTree as ET

# Word XML namespace
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}

def extract_all_zotero_data(docx_path):
    """Extract all Zotero citations and bibliography from docx."""

    with zipfile.ZipFile(docx_path, 'r') as zf:
        doc_xml = zf.read('word/document.xml').decode('utf-8')

    # Parse XML
    root = ET.fromstring(doc_xml)

    # Find all instrText elements
    instr_texts = root.findall('.//w:instrText', NS)
    print(f"Found {len(instr_texts)} instrText elements")

    citations = []
    bibliography = None

    for elem in instr_texts:
        text = elem.text or ""

        # Zotero citation
        if 'ZOTERO_ITEM CSL_CITATION' in text:
            # Extract the JSON part
            match = re.search(r'ADDIN ZOTERO_ITEM CSL_CITATION (\{.*)', text, re.DOTALL)
            if match:
                json_str = match.group(1).strip()
                try:
                    data = json.loads(json_str)
                    formatted = data.get('properties', {}).get('formattedCitation', 'unknown')
                    citations.append({
                        'formatted': formatted,
                        'citationID': data.get('citationID'),
                        'full_field': f'ADDIN ZOTERO_ITEM CSL_CITATION {json_str}',
                        'json_data': data
                    })
                    print(f"  Citation: {formatted} -> {data.get('citationID')}")
                except json.JSONDecodeError as e:
                    print(f"  Failed to parse JSON: {e}")
                    # Try to extract partial
                    citations.append({
                        'formatted': 'parse_error',
                        'raw_text': text[:200]
                    })

        # Zotero bibliography
        elif 'ZOTERO_BIBL' in text:
            match = re.search(r'ADDIN ZOTERO_BIBL (\{[^}]+\})', text)
            if match:
                bibliography = f'ADDIN ZOTERO_BIBL {match.group(1)}'
                print(f"  Bibliography field found")

    return citations, bibliography

def main():
    v3_path = Path('CITEgeistManuscript_v3.docx')

    if not v3_path.exists():
        print(f"ERROR: {v3_path} not found")
        return

    print(f"Extracting from {v3_path}...")
    citations, bib = extract_all_zotero_data(v3_path)

    print(f"\n=== Summary ===")
    print(f"Total citations found: {len(citations)}")
    if bib:
        print(f"Bibliography field: YES")

    # Save full extraction
    output = {
        'source': str(v3_path),
        'citation_count': len(citations),
        'citations': citations,
        'bibliography_field': bib
    }

    with open('zotero_full_extraction.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nSaved to zotero_full_extraction.json")

    # Also print citation mapping
    print("\n=== Citation Mapping ===")
    for i, c in enumerate(citations, 1):
        print(f"{i}. {c.get('formatted', '?')} -> ID: {c.get('citationID', '?')}")

if __name__ == '__main__':
    main()
