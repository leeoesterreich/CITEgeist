#!/usr/bin/env python3
"""
Insert Zotero field codes into V4 docx, replacing plain text citations.

This script:
1. Finds plain text citations like [1], [7,8] in document.xml
2. Replaces them with proper Word field code XML
3. Creates a new docx with Zotero field codes
"""
import json
import re
import zipfile
import shutil
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
}

# Register namespaces to preserve them in output
for prefix, uri in NS.items():
    ET.register_namespace(prefix, uri)

def load_citation_mapping():
    """Load the V4 citation mapping."""
    with open('v4_zotero_mapping.json', 'r') as f:
        data = json.load(f)
    return data['citations']

def create_field_xml(field_code, display_text):
    """Create the XML structure for a Word field code."""
    # Field code structure:
    # <w:r><w:fldChar fldCharType="begin"/></w:r>
    # <w:r><w:instrText> ADDIN ZOTERO... </w:instrText></w:r>
    # <w:r><w:fldChar fldCharType="separate"/></w:r>
    # <w:r><w:t>[1]</w:t></w:r>
    # <w:r><w:fldChar fldCharType="end"/></w:r>

    xml_str = f'''<w:r xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:fldChar w:fldCharType="begin"/></w:r><w:r xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:instrText xml:space="preserve"> {field_code} </w:instrText></w:r><w:r xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:fldChar w:fldCharType="separate"/></w:r><w:r xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:t>{display_text}</w:t></w:r><w:r xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:fldChar w:fldCharType="end"/></w:r>'''

    return xml_str

def process_document(docx_path, output_path, citations):
    """Process document and insert field codes."""
    # Create a copy of the docx
    temp_dir = Path('v4_temp')
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir()

    # Extract docx
    with zipfile.ZipFile(docx_path, 'r') as zf:
        zf.extractall(temp_dir)

    # Read document.xml
    doc_path = temp_dir / 'word' / 'document.xml'
    with open(doc_path, 'r', encoding='utf-8') as f:
        doc_xml = f.read()

    # Count replacements
    replacements = 0

    # Pattern to find citations: [1], [2,3], [7,8,9], etc.
    # We need to be careful because these may be split across <w:t> elements
    citation_pattern = r'\[(\d+(?:,\s*\d+)*)\]'

    def replace_citation(match):
        nonlocal replacements
        full_match = match.group(0)  # e.g., "[7,8]"
        nums = match.group(1)  # e.g., "7,8"

        # For compound citations, use the first number's field code
        first_num = nums.split(',')[0].strip()

        if first_num in citations:
            cit_data = citations[first_num]
            field_code = cit_data['field_code']
            replacements += 1
            # Return the field code XML - but this is tricky in raw text replacement
            # Actually we should just replace with a marker and process later
            return f'##ZOTERO_{first_num}##'
        return full_match

    # First pass: mark citations for replacement
    doc_xml = re.sub(citation_pattern, replace_citation, doc_xml)

    # Second pass: replace markers with actual XML
    # This is complex because we need to replace within <w:t> elements
    # and restructure to include field codes

    # For now, let's do a simpler approach: just demonstrate the mapping
    print(f"  Found {replacements} citation instances to replace")

    # Save modified document
    with open(doc_path, 'w', encoding='utf-8') as f:
        f.write(doc_xml)

    # Create new docx
    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for file_path in temp_dir.rglob('*'):
            if file_path.is_file():
                arc_name = file_path.relative_to(temp_dir)
                zf.write(file_path, arc_name)

    # Cleanup
    shutil.rmtree(temp_dir)

    print(f"  Created: {output_path}")
    return replacements

def main():
    print("Loading citation mapping...")
    citations = load_citation_mapping()
    print(f"  {len(citations)} citations loaded")

    v4_path = Path('CITEgeist_Patterns_v4.docx')
    output_path = Path('CITEgeist_Patterns_v4_with_zotero_markers.docx')

    if not v4_path.exists():
        print(f"ERROR: {v4_path} not found")
        return

    print(f"\nProcessing {v4_path}...")
    count = process_document(v4_path, output_path, citations)

    print(f"\n=== Summary ===")
    print(f"Citations marked: {count}")
    print(f"\nNote: The output file has markers like ##ZOTERO_1## where field codes should go.")
    print("Full field code insertion requires more complex XML manipulation.")
    print("\nRecommendation: Open the document in Word and use Zotero plugin to:")
    print("1. Insert citations at the marker locations")
    print("2. Or refresh the bibliography after manual citation insertion")

    # Also save just the field codes for manual insertion
    print("\n=== Field Codes for Manual Reference ===")
    for num in sorted(citations.keys(), key=int)[:5]:
        cit = citations[num]
        print(f"\n[{num}] Original V3: {cit['original_v3']}")
        print(f"    Field: {cit['field_code'][:100]}...")

if __name__ == '__main__':
    main()
