#!/bin/bash
"""
Setup HTML to PDF Converter
============================

This script installs the required Playwright browser for HTML to PDF conversion.
"""

echo "🔧 Setting up HTML to PDF converter..."

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is required but not installed"
    exit 1
fi

# Check if we're in a virtual environment (recommended)
if [[ "$VIRTUAL_ENV" != "" ]]; then
    echo "✅ Virtual environment detected: $VIRTUAL_ENV"
else
    echo "⚠️  Warning: Not in a virtual environment. Consider using one."
fi

# Install Python dependencies
echo "📦 Installing Python dependencies..."
pip install playwright==1.48.0

# Install Playwright browser
echo "🌐 Installing Playwright Chromium browser..."
playwright install chromium

# Install system dependencies (if needed)
echo "🔧 Installing system dependencies..."
playwright install-deps chromium

echo "✅ Setup complete!"
echo ""
echo "Usage examples:"
echo "  # Convert any HTML file to PDF"
echo "  python scripts/utils/html_to_pdf.py input.html output.pdf"
echo ""  
echo "  # Convert with high quality preset"
echo "  python scripts/utils/html_to_pdf.py input.html output.pdf --preset high_quality"
echo ""
echo "  # Generate PDF from the TEP-GNSS site"
echo "  python scripts/generate_site_pdf.py"
echo ""
echo "  # Custom options"
echo "  python scripts/utils/html_to_pdf.py input.html output.pdf --scale 1.5 --format A3 --wait-time 5"
