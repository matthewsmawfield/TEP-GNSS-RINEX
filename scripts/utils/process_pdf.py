#!/usr/bin/env python3
"""Unified PDF Processing Script
Compresses PDF and embeds project-specific metadata from CITATION.cff.

This script processes TEP manuscript PDFs by compressing them for web distribution
and embedding complete academic metadata for proper indexing and citation.
Metadata is auto-detected from the project's CITATION.cff file.

Usage:
    python process_pdf.py <input_pdf> [--quality ebook|printer|prepress|default]
    
Example:
    python process_pdf.py site/public/docs/Smawfield_2026_TEP-J0437_v0.1_Sintra.pdf --quality ebook
"""

import subprocess
import sys
import os
import re
from pathlib import Path
import argparse
import tempfile

try:
    import yaml
except ImportError:
    yaml = None


def compress_pdf(input_path, output_path, quality='ebook'):
    """Compress PDF using Ghostscript."""
    # Validate inputs are safe paths (prevent directory traversal)
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    
    # Ensure input exists and is a file
    if not input_path.is_file():
        raise ValueError(f"Input path must be an existing file: {input_path}")
    
    # Ensure output parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    quality_settings = {
        'screen': '/screen',      # 72 dpi
        'ebook': '/ebook',        # 150 dpi
        'printer': '/printer',    # 300 dpi
        'prepress': '/prepress',  # 300 dpi, color preserving
        'default': '/default'
    }
    
    if quality not in quality_settings:
        raise ValueError(f"Quality must be one of: {', '.join(quality_settings.keys())}")
    
    gs_quality = quality_settings[quality]
    
    # Get original size
    original_size = input_path.stat().st_size
    
    # Compress using Ghostscript
    cmd = [
        'gs',
        '-sDEVICE=pdfwrite',
        '-dCompatibilityLevel=1.4',
        f'-dPDFSETTINGS={gs_quality}',
        '-dNOPAUSE',
        '-dQUIET',
        '-dBATCH',
        f'-sOutputFile={output_path}',
        str(input_path)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        compressed_size = output_path.stat().st_size
        reduction = ((original_size - compressed_size) / original_size) * 100
        
        return {
            'original_mb': original_size / (1024 * 1024),
            'compressed_mb': compressed_size / (1024 * 1024),
            'reduction_pct': reduction
        }
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Ghostscript compression failed: {e.stderr.decode()}")


def load_citation_metadata(pdf_path):
    """Load metadata from CITATION.cff in the same project as the PDF."""
    pdf_path = Path(pdf_path).resolve()
    
    # Find project root (look for CITATION.cff)
    project_dir = pdf_path.parent
    while project_dir != project_dir.parent:
        citation_file = project_dir / 'CITATION.cff'
        if citation_file.exists():
            break
        project_dir = project_dir.parent
    else:
        # Fallback: assume 3 levels up from docs folder
        project_dir = pdf_path.parent.parent.parent.parent
        citation_file = project_dir / 'CITATION.cff'
    
    if not citation_file.exists():
        print(f"⚠️  CITATION.cff not found, using minimal defaults")
        return {
            'title': 'TEP Manuscript',
            'author': 'Matthew Lukin Smawfield',
            'version': 'v0.1',
            'codename': 'Unknown',
            'doi': '',
            'abstract': '',
            'keywords': [],
            'date': '2026-01-01',
            'url': '',
            'project_short': 'TEP'
        }
    
    try:
        if yaml:
            with open(citation_file, 'r') as f:
                data = yaml.safe_load(f)
        else:
            # Manual parsing fallback
            with open(citation_file, 'r') as f:
                content = f.read()
            data = {}
            # Extract key fields with regex
            title_match = re.search(r'title:\s*"([^"]+)"', content)
            data['title'] = title_match.group(1) if title_match else 'TEP Manuscript'
            
            version_match = re.search(r'version:\s*"?([^"\n]+)"?', content)
            data['version'] = version_match.group(1).strip() if version_match else 'v0.1'
            
            doi_match = re.search(r'doi:\s*"?([^"\n]+)"?', content)
            data['doi'] = doi_match.group(1).strip() if doi_match else ''
            
            abstract_match = re.search(r'abstract:\s*>?\s*"?([^"\n]+(?:\n[^"]+)*)"?', content, re.DOTALL)
            data['abstract'] = abstract_match.group(1).strip() if abstract_match else ''
            
            url_match = re.search(r'url:\s*"([^"]+)"', content)
            data['url'] = url_match.group(1) if url_match else ''
            
            date_match = re.search(r'date-released:\s*"?([^"\n]+)"?', content)
            data['date-released'] = date_match.group(1).strip() if date_match else '2026-01-01'
            
            # Parse authors
            author_matches = re.findall(r'family-names:\s*"([^"]+)"\s*\n\s*given-names:\s*"([^"]+)"', content)
            if author_matches:
                data['authors'] = [{'family-names': fam, 'given-names': giv} for fam, giv in author_matches]
            else:
                data['authors'] = [{'family-names': 'Smawfield', 'given-names': 'Matthew Lukin'}]
            
            # Parse keywords
            kw_section = re.search(r'keywords:\s*\n((?:\s+-\s*[^\n]+\n?)+)', content)
            if kw_section:
                keywords = re.findall(r'-\s*([^\n]+)', kw_section.group(1))
                data['keywords'] = [k.strip() for k in keywords]
            else:
                data['keywords'] = []
        
        # Extract version and codename
        version_str = data.get('version', 'v0.1')
        pattern = r'^(v?[\d.]+)(?:\s*\(([^)]+)\))?$'
        match = re.match(pattern, version_str.strip())
        
        if match:
            version = match.group(1).lstrip('v')
            codename = match.group(2) or 'Unknown'
        else:
            version = version_str.lstrip('v')
            codename = 'Unknown'
        
        # Get author
        authors = data.get('authors', [])
        author_name = "Matthew Lukin Smawfield"  # default
        if authors:
            first_author = authors[0]
            family = first_author.get('family-names', '')
            given = first_author.get('given-names', '')
            author_name = f"{given} {family}".strip()
        
        # Get date
        date = data.get('date-released', '2026-01-01')
        if hasattr(date, 'strftime'):
            date = date.strftime('%Y-%m-%d')
        else:
            date = str(date)
        
        # Get project short name from directory
        project_short = project_dir.name.replace('TEP-', '') if project_dir.name.startswith('TEP-') else 'TEP'
        
        # Get URL from data or construct from project name
        url = data.get('url', '')
        if not url and 'repository-code' in data:
            url = data['repository-code']
        
        return {
            'title': data.get('title', 'TEP Manuscript'),
            'author': author_name,
            'version': version,
            'codename': codename,
            'doi': data.get('doi', ''),
            'abstract': data.get('abstract', ''),
            'keywords': data.get('keywords', []),
            'date': date,
            'url': url,
            'project_short': project_short,
            'project_dir': project_dir
        }
        
    except Exception as e:
        print(f"⚠️  Error parsing CITATION.cff: {e}, using defaults")
        return {
            'title': 'TEP Manuscript',
            'author': 'Matthew Lukin Smawfield',
            'version': 'v0.1',
            'codename': 'Unknown',
            'doi': '',
            'abstract': '',
            'keywords': [],
            'date': '2026-01-01',
            'url': '',
            'project_short': 'TEP'
        }


def embed_metadata(pdf_path, metadata):
    """Embed metadata into PDF using exiftool."""
    cmd = ['exiftool']
    
    # Add all metadata fields
    for key, value in metadata.items():
        cmd.extend([f'-{key}={value}'])
    
    # Overwrite original
    cmd.extend(['-overwrite_original', pdf_path])
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Exiftool metadata embedding failed: {e.stderr.decode()}")


def build_metadata_dict(meta):
    """Build complete metadata dict from CITATION.cff data."""
    
    # Format keywords as semicolon-separated string
    keywords_str = '; '.join(meta['keywords']) if meta['keywords'] else 'TEP; Temporal Equivalence Principle'
    
    # Format creation date for PDF
    date_parts = meta['date'].split('-')
    creation_date = f"{date_parts[0]}:{date_parts[1]}:{date_parts[2]} 00:00:00" if len(date_parts) == 3 else meta['date']
    
    # Build subject/abstract
    subject = meta['abstract'][:2000] if meta['abstract'] else f"TEP manuscript: {meta['title']}"
    
    # Project identifier
    project_id = f"TEP-{meta['project_short']} v{meta['version']} ({meta['codename']})"
    
    # DOI with prefix
    doi_full = f"doi:{meta['doi']}" if meta['doi'] else ''
    
    metadata = {
        # Core identification
        'Title': meta['title'],
        'Author': meta['author'],
        'Creator': meta['author'],
        
        # Abstract/Subject
        'Subject': subject,
        
        # Keywords
        'Keywords': keywords_str,
        
        # Production metadata
        'Producer': project_id,
        
        # Rights
        'Copyright': 'Creative Commons Attribution 4.0 International License (CC BY 4.0)',
        
        # Dates
        'CreationDate': creation_date,
        'ModifyDate': creation_date,
    }
    
    # XMP Dublin Core
    metadata['XMP-dc:Creator'] = meta['author']
    metadata['XMP-dc:Title'] = meta['title']
    metadata['XMP-dc:Description'] = meta['abstract'][:500] if meta['abstract'] else meta['title']
    metadata['XMP-dc:Rights'] = 'CC BY 4.0'
    if doi_full:
        metadata['XMP-dc:Identifier'] = doi_full
    if meta['url']:
        metadata['XMP-dc:Source'] = meta['url']
    metadata['XMP-dc:Publisher'] = 'Zenodo'
    metadata['XMP-dc:Date'] = meta['date']
    metadata['XMP-dc:Type'] = 'Preprint'
    metadata['XMP-dc:Format'] = 'application/pdf'
    metadata['XMP-dc:Language'] = 'en'
    
    # PRISM metadata
    if meta['doi']:
        metadata['XMP-prism:DOI'] = meta['doi']
    if meta['url']:
        metadata['XMP-prism:URL'] = meta['url']
    metadata['XMP-prism:VersionIdentifier'] = meta['version']
    metadata['XMP-prism:PublicationName'] = 'TEP Research Series'
    
    # PDF/A
    metadata['XMP-pdfaid:Part'] = '1'
    metadata['XMP-pdfaid:Conformance'] = 'B'
    
    return metadata


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
    
    # Step 1: Compress PDF
    print("Step 1: Compressing PDF...")
    with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        stats = compress_pdf(str(input_path), tmp_path, args.quality)
        
        # Replace original with compressed version
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
    print("Step 2: Loading metadata from CITATION.cff...")
    
    # Load project-specific metadata from CITATION.cff
    citation_meta = load_citation_metadata(str(input_path))
    print(f"  Project: TEP-{citation_meta['project_short']}")
    print(f"  Title: {citation_meta['title'][:50]}...")
    print(f"  Version: v{citation_meta['version']} ({citation_meta['codename']})")
    print()
    
    print("Step 3: Embedding metadata...")
    
    # Build complete metadata dict
    metadata = build_metadata_dict(citation_meta)
    
    try:
        embed_metadata(str(input_path), metadata)
        print("  Metadata embedded successfully")
        print()
        
    except Exception as e:
        print(f"Error during metadata embedding: {e}")
        sys.exit(1)
    
    # Step 4: Verify
    print("Step 4: Verifying metadata...")
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
