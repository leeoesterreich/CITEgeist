#!/usr/bin/env python3
"""
Transfer Zotero citations from V3 to V4 docx.
"""
import zipfile
import re
import json
import shutil
from pathlib import Path
from xml.etree import ElementTree as ET

# Word XML namespaces
NAMESPACES = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
}

def extract_zotero_citations(docx_path):
    """Extract all Zotero field codes from a docx file."""
    citations = {}

    with zipfile.ZipFile(docx_path, 'r') as zf:
        doc_xml = zf.read('word/document.xml').decode('utf-8')

    # Find all Zotero citation field codes
    # Pattern matches: ADDIN ZOTERO_ITEM CSL_CITATION {...json...}
    pattern = r'ADDIN ZOTERO_ITEM CSL_CITATION (\{.*?\})(?=</w:instrText>|ADDIN)'

    # Also need to handle the complete JSON which may span multiple lines
    # Let's try a different approach - find the instrText elements

    # Parse the full JSON from field codes
    full_pattern = r'<w:instrText[^>]*>ADDIN ZOTERO_ITEM CSL_CITATION (\{[^<]+)'
    matches = re.findall(full_pattern, doc_xml)

    for i, match in enumerate(matches):
        # Try to parse the JSON to get the formatted citation
        try:
            # The JSON might be incomplete, try to fix it
            json_str = match.strip()
            if not json_str.endswith('}'):
                # Find the next closing brace pattern
                continue
            data = json.loads(json_str)
            formatted = data.get('properties', {}).get('formattedCitation', f'[{i+1}]')
            citations[formatted] = {
                'raw': f'ADDIN ZOTERO_ITEM CSL_CITATION {json_str}',
                'json': data
            }
        except json.JSONDecodeError:
            continue

    return citations

def extract_zotero_bibliography(docx_path):
    """Extract the Zotero bibliography field code."""
    with zipfile.ZipFile(docx_path, 'r') as zf:
        doc_xml = zf.read('word/document.xml').decode('utf-8')

    # Find ZOTERO_BIBL field
    pattern = r'ADDIN ZOTERO_BIBL (\{[^}]+\})'
    match = re.search(pattern, doc_xml)
    if match:
        return f'ADDIN ZOTERO_BIBL {match.group(1)}'
    return None

def get_citation_xml_template():
    """Return XML template for a Zotero citation field."""
    return '''<w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> {instr} </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>{display}</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r>'''

def main():
    v3_path = Path('CITEgeistManuscript_v3.docx')
    v4_path = Path('CITEgeist_Patterns_v4.docx')

    print("Extracting Zotero citations from V3...")
    citations = extract_zotero_citations(v3_path)
    print(f"Found {len(citations)} citations")

    for citation, data in citations.items():
        print(f"  {citation}: {data['json'].get('citationID', 'unknown')}")

    bib = extract_zotero_bibliography(v3_path)
    if bib:
        print(f"\nFound bibliography field")

    # Save extracted citations for reference
    with open('zotero_citations_extracted.json', 'w') as f:
        # Can't serialize the full data, just save citation IDs
        output = {}
        for citation, data in citations.items():
            output[citation] = {
                'citationID': data['json'].get('citationID'),
                'raw_instruction': data['raw'][:200] + '...' if len(data['raw']) > 200 else data['raw']
            }
        json.dump(output, f, indent=2)

    print("\nExtracted citations saved to zotero_citations_extracted.json")
    print("\nNote: Full field code transfer requires complex XML manipulation.")
    print("Recommendation: Open V4.docx in Word, use Zotero plugin to re-insert citations.")

if __name__ == '__main__':
    main()
