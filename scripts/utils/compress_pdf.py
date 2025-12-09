#!/usr/bin/env python3
"""
PDF Compression Utility
Uses Ghostscript for high-quality compression.

Usage:
    python compress_pdf.py input.pdf                    # Creates input_compressed.pdf
    python compress_pdf.py input.pdf output.pdf         # Custom output name
    python compress_pdf.py input.pdf -q screen          # Lower quality, smaller size
    python compress_pdf.py input.pdf -q prepress        # Higher quality, larger size

Quality presets (from smallest to largest):
    screen   - 72 dpi, lowest quality, smallest size
    ebook    - 150 dpi, good for viewing (default)
    printer  - 300 dpi, high quality
    prepress - 300 dpi, highest quality, largest size
"""

import argparse
import subprocess
import sys
from pathlib import Path


def compress_pdf(input_path: str, output_path: str = None, quality: str = "ebook") -> dict:
    """
    Compress a PDF file using Ghostscript.
    
    Args:
        input_path: Path to input PDF
        output_path: Path to output PDF (default: input_compressed.pdf)
        quality: Compression preset (screen, ebook, printer, prepress)
    
    Returns:
        dict with original_size, compressed_size, reduction_percent
    """
    input_file = Path(input_path)
    
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    if not input_file.suffix.lower() == '.pdf':
        raise ValueError(f"Input must be a PDF file: {input_path}")
    
    # Default output name
    if output_path is None:
        output_file = input_file.parent / f"{input_file.stem}_compressed.pdf"
    else:
        output_file = Path(output_path)
    
    # Validate quality preset
    valid_presets = ["screen", "ebook", "printer", "prepress"]
    if quality not in valid_presets:
        raise ValueError(f"Quality must be one of: {valid_presets}")
    
    # Run Ghostscript
    cmd = [
        "gs",
        "-sDEVICE=pdfwrite",
        "-dCompatibilityLevel=1.4",
        f"-dPDFSETTINGS=/{quality}",
        "-dNOPAUSE",
        "-dQUIET",
        "-dBATCH",
        f"-sOutputFile={output_file}",
        str(input_file)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        raise RuntimeError(
            "Ghostscript not found. Install with:\n"
            "  macOS:  brew install ghostscript\n"
            "  Ubuntu: sudo apt install ghostscript\n"
            "  Windows: https://ghostscript.com/releases/gsdnld.html"
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Ghostscript error: {e.stderr}")
    
    # Calculate sizes
    orig_size = input_file.stat().st_size
    new_size = output_file.stat().st_size
    reduction = (1 - new_size / orig_size) * 100
    
    return {
        "input": str(input_file),
        "output": str(output_file),
        "original_size": orig_size,
        "compressed_size": new_size,
        "reduction_percent": reduction
    }


def format_size(bytes_size: int) -> str:
    """Format bytes as human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_size < 1024:
            return f"{bytes_size:.2f} {unit}"
        bytes_size /= 1024
    return f"{bytes_size:.2f} TB"


def main():
    parser = argparse.ArgumentParser(
        description="Compress PDF files using Ghostscript",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("input", help="Input PDF file")
    parser.add_argument("output", nargs="?", help="Output PDF file (optional)")
    parser.add_argument(
        "-q", "--quality",
        choices=["screen", "ebook", "printer", "prepress"],
        default="ebook",
        help="Compression quality preset (default: ebook)"
    )
    
    args = parser.parse_args()
    
    try:
        result = compress_pdf(args.input, args.output, args.quality)
        
        print(f"Input:       {result['input']}")
        print(f"Output:      {result['output']}")
        print(f"Original:    {format_size(result['original_size'])}")
        print(f"Compressed:  {format_size(result['compressed_size'])}")
        print(f"Reduction:   {result['reduction_percent']:.1f}%")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
