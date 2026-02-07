#!/usr/bin/env python3
"""
Properly insert Zotero field codes into V4 docx.

This handles the XML structure correctly, replacing <w:t>[1]</w:t> with
the full field code structure.
"""
import json
import re
import zipfile
import shutil
from pathlib import Path
from lxml import etree

NS = {
    'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
}

def load_citation_mapping():
    """Load the V4 citation mapping."""
    with open('v4_zotero_mapping.json', 'r') as f:
        data = json.load(f)
    return data['citations']

def create_field_elements(field_code, display_text, nsmap):
    """Create the XML elements for a Word field code."""
    w = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'

    elements = []

    # Begin field
    r1 = etree.Element(f'{w}r', nsmap=nsmap)
    fld_begin = etree.SubElement(r1, f'{w}fldChar')
    fld_begin.set(f'{w}fldCharType', 'begin')
    elements.append(r1)

    # Instruction text
    r2 = etree.Element(f'{w}r', nsmap=nsmap)
    instr = etree.SubElement(r2, f'{w}instrText')
    instr.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    instr.text = f' {field_code} '
    elements.append(r2)

    # Separator
    r3 = etree.Element(f'{w}r', nsmap=nsmap)
    fld_sep = etree.SubElement(r3, f'{w}fldChar')
    fld_sep.set(f'{w}fldCharType', 'separate')
    elements.append(r3)

    # Display text
    r4 = etree.Element(f'{w}r', nsmap=nsmap)
    t = etree.SubElement(r4, f'{w}t')
    t.text = display_text
    elements.append(r4)

    # End field
    r5 = etree.Element(f'{w}r', nsmap=nsmap)
    fld_end = etree.SubElement(r5, f'{w}fldChar')
    fld_end.set(f'{w}fldCharType', 'end')
    elements.append(r5)

    return elements

def process_text_element(t_elem, citations, nsmap):
    """Process a <w:t> element, replacing citations with field codes."""
    if t_elem.text is None:
        return None

    text = t_elem.text
    citation_pattern = r'\[(\d+(?:,\s*\d+)*)\]'

    # Find all citations in this text
    matches = list(re.finditer(citation_pattern, text))
    if not matches:
        return None

    # Create replacement elements
    replacements = []
    last_end = 0

    for match in matches:
        # Text before the citation
        if match.start() > last_end:
            prefix = text[last_end:match.start()]
            if prefix:
                replacements.append(('text', prefix))

        # The citation itself
        full_citation = match.group(0)  # e.g., "[7,8]"
        nums = match.group(1)  # e.g., "7,8"
        first_num = nums.split(',')[0].strip()

        if first_num in citations:
            cit_data = citations[first_num]
            replacements.append(('field', cit_data['field_code'], full_citation))
        else:
            # Keep as plain text if not in mapping
            replacements.append(('text', full_citation))

        last_end = match.end()

    # Text after the last citation
    if last_end < len(text):
        suffix = text[last_end:]
        if suffix:
            replacements.append(('text', suffix))

    return replacements

def process_document(docx_path, output_path, citations):
    """Process document and insert field codes."""
    temp_dir = Path('v4_temp')
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir()

    # Extract docx
    with zipfile.ZipFile(docx_path, 'r') as zf:
        zf.extractall(temp_dir)

    # Parse document.xml with lxml
    doc_path = temp_dir / 'word' / 'document.xml'
    parser = etree.XMLParser(remove_blank_text=False)
    tree = etree.parse(str(doc_path), parser)
    root = tree.getroot()
    nsmap = root.nsmap

    w = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'

    # Find all <w:t> elements
    t_elements = root.findall(f'.//{w}t', nsmap)
    print(f"  Found {len(t_elements)} text elements")

    replacements_made = 0

    for t_elem in t_elements:
        replacements = process_text_element(t_elem, citations, nsmap)
        if replacements is None:
            continue

        # Get parent <w:r> element
        r_elem = t_elem.getparent()
        if r_elem is None or r_elem.tag != f'{w}r':
            continue

        # Get parent of <w:r> (e.g., <w:p>)
        parent = r_elem.getparent()
        if parent is None:
            continue

        # Get index of this <w:r> in parent
        idx = list(parent).index(r_elem)

        # Build new elements
        new_elements = []
        for repl in replacements:
            if repl[0] == 'text':
                # Create a new <w:r><w:t>text</w:t></w:r>
                new_r = etree.Element(f'{w}r', nsmap=nsmap)
                # Copy run properties if they exist
                rPr = r_elem.find(f'{w}rPr', nsmap)
                if rPr is not None:
                    new_r.append(rPr.__copy__())
                new_t = etree.SubElement(new_r, f'{w}t')
                new_t.text = repl[1]
                new_elements.append(new_r)
            elif repl[0] == 'field':
                field_code = repl[1]
                display = repl[2]
                field_elems = create_field_elements(field_code, display, nsmap)
                new_elements.extend(field_elems)
                replacements_made += 1

        # Replace the original <w:r> with new elements
        for i, new_elem in enumerate(new_elements):
            parent.insert(idx + i, new_elem)

        # Remove original <w:r>
        parent.remove(r_elem)

    print(f"  Made {replacements_made} field code insertions")

    # Write modified document
    tree.write(str(doc_path), xml_declaration=True, encoding='UTF-8', standalone=True)

    # Create new docx
    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for file_path in temp_dir.rglob('*'):
            if file_path.is_file():
                arc_name = file_path.relative_to(temp_dir)
                zf.write(file_path, arc_name)

    shutil.rmtree(temp_dir)
    print(f"  Created: {output_path}")
    return replacements_made

def main():
    print("Loading citation mapping...")
    citations = load_citation_mapping()
    print(f"  {len(citations)} citations loaded")

    v4_path = Path('CITEgeist_Patterns_v4.docx')
    output_path = Path('CITEgeist_Patterns_v4_zotero.docx')

    if not v4_path.exists():
        print(f"ERROR: {v4_path} not found")
        return

    print(f"\nProcessing {v4_path}...")
    count = process_document(v4_path, output_path, citations)

    print(f"\n=== Summary ===")
    print(f"Field codes inserted: {count}")
    print(f"Output: {output_path}")
    print("\nTo complete the Zotero integration in Word:")
    print("1. Open the document in Word")
    print("2. Open Zotero plugin")
    print("3. Click 'Refresh' to update citations")
    print("4. Add the bibliography using Zotero plugin")

if __name__ == '__main__':
    main()
