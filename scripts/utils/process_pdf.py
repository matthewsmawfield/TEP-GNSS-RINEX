#!/usr/bin/env python3
"""Generic PDF Processing Script
Compresses PDF and embeds metadata from CITATION.cff in one operation.

Reads title, abstract, keywords, version, date, DOI, URL, and authors
from the project's CITATION.cff file, then compresses the PDF with
Ghostscript and writes the metadata using exiftool.

Usage:
    python process_pdf.py <input_pdf> [--quality ebook|printer|prepress|default]

Example:
    python process_pdf.py site/public/docs/26-TEP-C0-v0.1-Athens.pdf --quality ebook
"""

import subprocess
import sys
import os
import re
from pathlib import Path
import argparse
import tempfile

from compress_pdf import compress_pdf as _compress_pdf

try:
    import yaml
except ImportError:
    yaml = None


def _extract_yaml_value(content, key):
    """Extract a scalar string value from CFF content."""
    # Match key: 'value' or key: value (non-greedy, single line)
    pattern = rf'^{re.escape(key)}:\s*[\'"]?(.+?)[\'"]?\s*$'
    match = re.search(pattern, content, re.MULTILINE)
    if match:
        return match.group(1).strip()
    return ''


def parse_citation_cff():
    """Parse CITATION.cff for PDF metadata."""
    base_dir = Path(__file__).parent.parent
    citation_file = base_dir / 'CITATION.cff'

    if not citation_file.exists():
        return None

    try:
        content = citation_file.read_text()

        if yaml:
            data = yaml.safe_load(content)
        else:
            # Manual parsing fallback
            data = {}
            data['title'] = _extract_yaml_value(content, 'title')
            data['version'] = _extract_yaml_value(content, 'version')
            data['date-released'] = _extract_yaml_value(content, 'date-released')
            data['doi'] = _extract_yaml_value(content, 'doi')
            data['url'] = _extract_yaml_value(content, 'url')
            data['license'] = _extract_yaml_value(content, 'license')

            # Abstract may span multiple lines in CFF
            abstract_match = re.search(
                r'^abstract:\s*[\'"]?((?:.*?\n)+?)(?:\n\w|\Z)',
                content, re.MULTILINE
            )
            if abstract_match:
                data['abstract'] = abstract_match.group(1).strip().replace("\n  ", " ")
            else:
                data['abstract'] = _extract_yaml_value(content, 'abstract')

            # Keywords list
            kw_match = re.search(r'^keywords:\n((?:- .+\n)+)', content, re.MULTILINE)
            if kw_match:
                kws = re.findall(r'- (.+)', kw_match.group(1))
                data['keywords'] = [k.strip() for k in kws]
            else:
                data['keywords'] = []

            # Authors
            authors = re.findall(r'family-names:\s*(.+)', content)
            givens = re.findall(r'given-names:\s*(.+)', content)
            data['authors'] = [
                {'family-names': f.strip(), 'given-names': g.strip()}
                for f, g in zip(authors, givens)
            ]

        return data
    except Exception:
        return None


def build_metadata(cff_data):
    """Build PDF metadata dict from parsed CFF data."""
    title = cff_data.get('title', '')

    # Author name
    authors = cff_data.get('authors', [])
    if authors and isinstance(authors, list) and len(authors) > 0:
        fam = authors[0].get('family-names', 'Smawfield')
        giv = authors[0].get('given-names', 'Matthew Lukin')
        author_name = f"{giv} {fam}"
    else:
        author_name = 'Matthew Lukin Smawfield'

    version = cff_data.get('version', 'v0.1')
    v_match = re.match(r'v?([\d.]+)(?:\s*\(([^)]+)\))?', str(version))
    if v_match:
        version_num = v_match.group(1)
        codename = v_match.group(2) or ''
    else:
        version_num = str(version).lstrip('v')
        codename = ''

    date = cff_data.get('date-released', '')
    if date:
        date_pdf = date.replace('-', ':')
    else:
        date_pdf = ''

    doi = cff_data.get('doi', '')
    url = cff_data.get('url', '')
    abstract = cff_data.get('abstract', '')
    license_str = cff_data.get('license', 'CC-BY-4.0')

    keywords_list = cff_data.get('keywords', [])
    if isinstance(keywords_list, list):
        keywords = '; '.join(str(k) for k in keywords_list)
    else:
        keywords = str(keywords_list)

    producer_label = f"TEP Research Project - Version {version_num}"
    if codename:
        producer_label += f" ({codename})"

    metadata = {
        'Title': title,
        'Author': author_name,
        'Creator': author_name,
        'Subject': abstract,
        'Keywords': keywords,
        'Producer': producer_label,
        'Copyright': f'Creative Commons Attribution 4.0 International License ({license_str})',
    }

    if date_pdf:
        metadata['CreationDate'] = f'{date_pdf} 00:00:00'
        metadata['ModifyDate'] = f'{date_pdf} 00:00:00'

    metadata['XMP-dc:Creator'] = author_name
    metadata['XMP-dc:Title'] = title
    metadata['XMP-dc:Description'] = abstract[:500] if abstract else ''
    metadata['XMP-dc:Rights'] = license_str
    metadata['XMP-dc:Publisher'] = 'Zenodo'
    metadata['XMP-dc:Type'] = 'Preprint'
    metadata['XMP-dc:Format'] = 'application/pdf'
    metadata['XMP-dc:Language'] = 'en'

    if doi:
        metadata['XMP-dc:Identifier'] = f'doi:{doi}'
        metadata['XMP-prism:DOI'] = doi
    if url:
        metadata['XMP-dc:Source'] = url
        metadata['XMP-prism:URL'] = url
    if date:
        metadata['XMP-dc:Date'] = date
    if version_num:
        metadata['XMP-prism:VersionIdentifier'] = version_num

    metadata['XMP-prism:PublicationName'] = 'TEP Research Series'
    metadata['XMP-pdfaid:Part'] = '1'
    metadata['XMP-pdfaid:Conformance'] = 'B'

    # Remove empty values
    metadata = {k: v for k, v in metadata.items() if v}
    return metadata


def compress_pdf(input_path, output_path, quality='ebook'):
    """Compress PDF using Ghostscript (wrapper with legacy return format)."""
    result = _compress_pdf(input_path, output_path, quality=quality)
    return {
        'original_mb': result['original_size'] / (1024 * 1024),
        'compressed_mb': result['compressed_size'] / (1024 * 1024),
        'reduction_pct': result['reduction_percent']
    }


def embed_metadata(pdf_path, metadata):
    """Embed metadata into PDF using exiftool."""
    cmd = ['exiftool']
    for key, value in metadata.items():
        cmd.extend([f'-{key}={value}'])
    cmd.extend(['-overwrite_original', pdf_path])
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode() if e.stderr else str(e)
        raise RuntimeError(f"Exiftool metadata embedding failed: {stderr}")


def verify_metadata(pdf_path, expected_fields):
    """Verify metadata was embedded correctly."""
    cmd = ['exiftool'] + [f'-{field}' for field in expected_fields] + [pdf_path]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Compress PDF and embed metadata in one operation'
    )
    parser.add_argument('input_pdf', help='Path to input PDF file')
    parser.add_argument(
        '--quality',
        choices=['screen', 'ebook', 'printer', 'prepress', 'default'],
        default='ebook',
        help='Compression quality (default: ebook)'
    )

    args = parser.parse_args()

    input_path = Path(args.input_pdf).resolve()

    if not input_path.exists():
        print(f"Error: File not found: {input_path}")
        sys.exit(1)

    print(f"Processing: {input_path}")
    print(f"Quality: {args.quality}")
    print()

    # Load metadata from CITATION.cff
    cff_data = parse_citation_cff()
    if cff_data:
        metadata = build_metadata(cff_data)
        print(f"Loaded metadata from CITATION.cff")
        print(f"  Title: {metadata.get('Title', 'N/A')[:60]}...")
    else:
        print("Warning: Could not load CITATION.cff, using minimal metadata")
        metadata = {
            'Title': 'TEP Manuscript',
            'Author': 'Matthew Lukin Smawfield',
            'Creator': 'Matthew Lukin Smawfield',
        }
    print()

    # Step 1: Compress PDF
    print("Step 1: Compressing PDF...")
    with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        stats = compress_pdf(str(input_path), tmp_path, args.quality)
        os.replace(tmp_path, str(input_path))
        print(f"  Original:    {stats['original_mb']:.2f} MB")
        print(f"  Compressed:  {stats['compressed_mb']:.2f} MB")
        print(f"  Reduction:   {stats['reduction_pct']:.1f}%")
        print()

    except Exception as e:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        print(f"Error during compression: {e}")
        sys.exit(1)

    # Step 2: Embed metadata
    print("Step 2: Embedding metadata...")
    try:
        embed_metadata(str(input_path), metadata)
        print("  Metadata embedded successfully")
        print()
    except Exception as e:
        print(f"Error during metadata embedding: {e}")
        sys.exit(1)

    # Step 3: Verify
    print("Step 3: Verifying metadata...")
    verification = verify_metadata(
        str(input_path),
        ['Title', 'Author', 'Subject', 'Keywords', 'Creator', 'Copyright']
    )

    if verification:
        print("  ✓ Metadata verified")
        print()
        print("Verification output:")
        print(verification)
    else:
        print("  ⚠ Could not verify metadata")

    print()
    print(f"✓ Processing complete: {input_path}")
    print(f"  Final size: {os.path.getsize(input_path) / (1024 * 1024):.2f} MB")


if __name__ == '__main__':
    main()
