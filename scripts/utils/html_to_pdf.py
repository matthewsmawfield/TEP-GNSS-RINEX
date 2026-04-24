#!/usr/bin/env python3
"""
Advanced HTML to PDF Converter
=============================

High-quality HTML to PDF conversion with comprehensive control over:
- Image resolution and quality
- Page breaks and layout
- Zoom/scale factors
- Print media CSS
- Custom page dimensions
- Headers and footers
- Margins and styling

Dependencies:
    pip install playwright
    playwright install chromium

Usage:
    python html_to_pdf.py input.html output.pdf --scale 1.5 --dpi 300
"""

import argparse
import asyncio
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Union
from urllib.parse import urljoin

from playwright.async_api import async_playwright

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class HTMLToPDFConverter:
    """Advanced HTML to PDF converter using Playwright."""
    
    def __init__(self):
        self.browser = None
        self.context = None
        
    async def __aenter__(self):
        """Async context manager entry."""
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(
            headless=True,
            args=[
                '--no-sandbox',
                '--disable-dev-shm-usage',
                '--disable-gpu',
                '--disable-extensions',
                '--font-render-hinting=none',  # Better font rendering
                '--enable-font-antialiasing',
                '--disable-background-timer-throttling',
            ]
        )
        return self
        
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        if self.browser:
            await self.browser.close()
        if self.playwright:
            await self.playwright.stop()
    
    def _get_css_for_pdf(self, options: Dict) -> str:
        """Generate CSS optimizations for PDF output."""
        css = """
        <style>
        @media print {
            /* Ensure high-quality images */
            img {
                max-width: 100% !important;
                height: auto !important;
                image-rendering: -webkit-optimize-contrast !important;
                image-rendering: crisp-edges !important;
                print-color-adjust: exact !important;
                -webkit-print-color-adjust: exact !important;
            }
            
            /* Better page breaks */
            .page-break-before { page-break-before: always !important; }
            .page-break-after { page-break-after: always !important; }
            .page-break-avoid { page-break-inside: avoid !important; }
            
            /* Avoid breaking these elements */
            h1, h2, h3, h4, h5, h6 { 
                page-break-after: avoid !important;
                page-break-inside: avoid !important;
            }
            
            figure, table, .figure, .table {
                page-break-inside: avoid !important;
                margin: 1em 0 !important;
            }
            
            /* Improve text rendering */
            body {
                -webkit-font-smoothing: subpixel-antialiased !important;
                font-smooth: always !important;
                text-rendering: optimizeLegibility !important;
            }
            
            /* Ensure backgrounds print */
            * {
                print-color-adjust: exact !important;
                -webkit-print-color-adjust: exact !important;
            }
        }
        </style>
        """
        
        # Add custom CSS if provided
        if options.get('custom_css'):
            css += f"<style>{options['custom_css']}</style>"
            
        return css
    
    async def convert_file(
        self, 
        input_path: Union[str, Path], 
        output_path: Union[str, Path],
        options: Optional[Dict] = None
    ) -> bool:
        """
        Convert HTML file to PDF.
        
        Args:
            input_path: Path to HTML file
            output_path: Path for output PDF
            options: Conversion options
            
        Returns:
            True if successful, False otherwise
        """
        if options is None:
            options = {}
            
        input_path = Path(input_path)
        output_path = Path(output_path)
        
        if not input_path.exists():
            logger.error(f"Input file not found: {input_path}")
            return False
            
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # Create browser context with device emulation if specified
            context_options = {
                'viewport': options.get('viewport', {'width': 1920, 'height': 1080}),
                'device_scale_factor': options.get('device_scale_factor', 1.0),
                'has_touch': False,
                'is_mobile': False,
            }
            
            self.context = await self.browser.new_context(**context_options)
            page = await self.context.new_page()
            
            # Handle network issues more gracefully
            page.on('requestfailed', lambda request: 
                logger.warning(f"Request failed: {request.url} - {request.failure}")
            )
            
            # Navigate to the HTML file
            file_url = input_path.resolve().as_uri()
            logger.info(f"Loading HTML file: {file_url}")
            
            await page.goto(file_url, wait_until='networkidle', timeout=60000)
            
            # Inject CSS for better PDF rendering
            css_content = self._get_css_for_pdf(options)
            await page.add_style_tag(content=css_content)
            
            # Wait for any dynamic content
            if options.get('wait_for_selector'):
                logger.info(f"Waiting for selector: {options['wait_for_selector']}")
                await page.wait_for_selector(options['wait_for_selector'], timeout=60000)
                # Additional wait to ensure all content is loaded
                await page.wait_for_timeout(2000)
            elif options.get('wait_time'):
                logger.info(f"Waiting {options['wait_time']} seconds for content to load")
                await page.wait_for_timeout(options['wait_time'] * 1000)
            else:
                # Default wait for images to load
                await page.wait_for_load_state('networkidle')
            
            # Additional check for dynamic content completion
            try:
                # Wait for the loading indicator to disappear (if present)
                loading_hidden = await page.wait_for_function(
                    "() => { const loading = document.getElementById('loading'); return !loading || loading.style.display === 'none' || loading.offsetParent === null; }",
                    timeout=30000
                )
                logger.info("Loading indicator hidden, content should be ready")
            except:
                logger.warning("Could not detect loading indicator state, proceeding anyway")
            
            # Final check for content presence
            try:
                content_ready = await page.wait_for_function(
                    "() => { const content = document.getElementById('content'); return content && content.children.length > 0; }",
                    timeout=10000
                )
                logger.info("Content container has content loaded")
            except:
                logger.warning("Could not verify content loading, proceeding anyway")
            
            # Wait specifically for images to load
            try:
                images_loaded = await page.wait_for_function(
                    """() => {
                        const images = document.querySelectorAll('img');
                        if (images.length === 0) return true;
                        return Array.from(images).every(img => 
                            img.complete && img.naturalHeight !== 0
                        );
                    }""",
                    timeout=15000
                )
                logger.info(f"All images loaded successfully")
            except:
                logger.warning("Could not verify all images loaded, proceeding anyway")
                # Log how many images we found
                try:
                    image_count = await page.evaluate("document.querySelectorAll('img').length")
                    logger.info(f"Found {image_count} images in document")
                except:
                    pass
            
            # Configure PDF options
            pdf_options = self._build_pdf_options(options)
            
            logger.info(f"Generating PDF with options: {json.dumps(pdf_options, indent=2)}")
            
            # Generate PDF
            await page.pdf(path=str(output_path), **pdf_options)
            
            logger.info(f"PDF generated successfully: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error converting HTML to PDF: {e}")
            return False
        finally:
            if self.context:
                await self.context.close()
    
    def _build_pdf_options(self, options: Dict) -> Dict:
        """Build PDF generation options."""
        pdf_options = {
            'format': options.get('format', 'A4'),
            'print_background': True,
            'prefer_css_page_size': options.get('prefer_css_page_size', False),
            'margin': {
                'top': options.get('margin_top', '1cm'),
                'right': options.get('margin_right', '1cm'), 
                'bottom': options.get('margin_bottom', '1cm'),
                'left': options.get('margin_left', '1cm'),
            }
        }
        
        # Scale/zoom
        if options.get('scale'):
            pdf_options['scale'] = options['scale']
            
        # Custom page dimensions
        if options.get('width') and options.get('height'):
            pdf_options['width'] = options['width']
            pdf_options['height'] = options['height']
            
        # Headers and footers
        if options.get('display_header_footer'):
            pdf_options['display_header_footer'] = True
            pdf_options['header_template'] = options.get(
                'header_template', 
                '<div style="font-size:10px; text-align:center; width:100%;">'
                '<span class="title"></span></div>'
            )
            pdf_options['footer_template'] = options.get(
                'footer_template',
                '<div style="font-size:10px; text-align:center; width:100%;">'
                'Page <span class="pageNumber"></span> of <span class="totalPages"></span></div>'
            )
            
        # Page ranges
        if options.get('page_ranges'):
            pdf_options['page_ranges'] = options['page_ranges']
            
        return pdf_options


def create_preset_configs() -> Dict[str, Dict]:
    """Create preset configurations for common use cases."""
    return {
        'high_quality': {
            'scale': 1.0,
            'device_scale_factor': 2.0,
            'format': 'A4',
            'margin_top': '1.5cm',
            'margin_bottom': '1.5cm',
            'margin_left': '2cm',
            'margin_right': '2cm',
            'prefer_css_page_size': True,
            'wait_time': 3,
            'custom_css': '''
                img { image-rendering: -webkit-optimize-contrast !important; }
                * { font-smooth: always !important; }
            '''
        },
        'print_ready': {
            'scale': 1.0,
            'format': 'A4',
            'display_header_footer': True,
            'header_template': '<div></div>',  # Empty header
            'footer_template': '''
                <div style="font-size:9px; text-align:center; width:100%; 
                           font-family: system-ui, -apple-system, sans-serif;">
                    Page <span class="pageNumber"></span> of <span class="totalPages"></span>
                </div>
            ''',
            'margin_top': '2cm',
            'margin_bottom': '2cm',
            'margin_left': '2cm',
            'margin_right': '2cm',
        },
        'web_optimized': {
            'scale': 0.8,
            'format': 'A4',
            'margin_top': '1cm',
            'margin_bottom': '1cm',
            'margin_left': '1cm',
            'margin_right': '1cm',
            'wait_time': 2,
        },
        'large_format': {
            'format': 'A3',
            'scale': 1.2,
            'device_scale_factor': 1.5,
            'margin_top': '2cm',
            'margin_bottom': '2cm',
            'margin_left': '2.5cm',
            'margin_right': '2.5cm',
        }
    }


async def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description='Advanced HTML to PDF converter',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion
  python html_to_pdf.py input.html output.pdf
  
  # High quality with custom scale
  python html_to_pdf.py input.html output.pdf --preset high_quality --scale 1.5
  
  # Custom dimensions and margins
  python html_to_pdf.py input.html output.pdf --format A3 --margin-top 2cm --scale 1.2
  
  # With headers and footers
  python html_to_pdf.py input.html output.pdf --header-footer --wait-time 5
        """
    )
    
    # Required arguments
    parser.add_argument('input', help='Input HTML file path')
    parser.add_argument('output', help='Output PDF file path')
    
    # Preset configurations
    parser.add_argument('--preset', choices=create_preset_configs().keys(),
                       help='Use preset configuration')
    
    # Page format and dimensions
    parser.add_argument('--format', default='A4',
                       help='Page format (A4, A3, Letter, Legal, Tabloid)')
    parser.add_argument('--width', help='Custom page width (e.g., 21cm, 8.5in)')
    parser.add_argument('--height', help='Custom page height (e.g., 29.7cm, 11in)')
    
    # Scaling and quality
    parser.add_argument('--scale', type=float, default=1.0,
                       help='Scale factor for the page (0.1 to 2.0)')
    parser.add_argument('--dpi', type=int, default=96,
                       help='Device scale factor equivalent (96, 150, 192, 300)')
    
    # Margins
    parser.add_argument('--margin-top', default='1cm', help='Top margin')
    parser.add_argument('--margin-right', default='1cm', help='Right margin')
    parser.add_argument('--margin-bottom', default='1cm', help='Bottom margin')
    parser.add_argument('--margin-left', default='1cm', help='Left margin')
    
    # Headers and footers
    parser.add_argument('--header-footer', action='store_true',
                       help='Enable headers and footers')
    parser.add_argument('--header-template', help='Custom header HTML template')
    parser.add_argument('--footer-template', help='Custom footer HTML template')
    
    # Wait options
    parser.add_argument('--wait-time', type=float, help='Wait time in seconds before PDF generation')
    parser.add_argument('--wait-for-selector', help='CSS selector to wait for before PDF generation')
    
    # Advanced options
    parser.add_argument('--page-ranges', help='Page ranges to print (e.g., "1-3, 5, 8-")')
    parser.add_argument('--css-page-size', action='store_true',
                       help='Prefer CSS page size over format')
    parser.add_argument('--custom-css-file', help='Path to custom CSS file for PDF generation')
    
    # Viewport
    parser.add_argument('--viewport-width', type=int, default=1920,
                       help='Browser viewport width')
    parser.add_argument('--viewport-height', type=int, default=1080,
                       help='Browser viewport height')
    
    args = parser.parse_args()
    
    # Build options from arguments
    options = {}
    
    # Apply preset if specified
    if args.preset:
        presets = create_preset_configs()
        options.update(presets[args.preset])
        logger.info(f"Applied preset: {args.preset}")
    
    # Override with command line arguments
    options.update({
        'format': args.format,
        'scale': args.scale,
        'device_scale_factor': max(1.0, args.dpi / 96.0),
        'margin_top': args.margin_top,
        'margin_right': args.margin_right,
        'margin_bottom': args.margin_bottom,
        'margin_left': args.margin_left,
        'prefer_css_page_size': args.css_page_size,
        'viewport': {
            'width': args.viewport_width,
            'height': args.viewport_height
        }
    })
    
    if args.width and args.height:
        options.update({'width': args.width, 'height': args.height})
    
    if args.header_footer:
        options['display_header_footer'] = True
        if args.header_template:
            options['header_template'] = args.header_template
        if args.footer_template:
            options['footer_template'] = args.footer_template
    
    if args.wait_time:
        options['wait_time'] = args.wait_time
    
    if args.wait_for_selector:
        options['wait_for_selector'] = args.wait_for_selector
        
    if args.page_ranges:
        options['page_ranges'] = args.page_ranges
    
    if args.custom_css_file:
        css_path = Path(args.custom_css_file)
        if css_path.exists():
            options['custom_css'] = css_path.read_text()
        else:
            logger.warning(f"Custom CSS file not found: {css_path}")
    
    # Convert HTML to PDF
    async with HTMLToPDFConverter() as converter:
        success = await converter.convert_file(args.input, args.output, options)
        
        if success:
            output_size = Path(args.output).stat().st_size / (1024 * 1024)
            logger.info(f"✅ Conversion completed successfully")
            logger.info(f"📄 Output: {args.output} ({output_size:.1f} MB)")
            return 0
        else:
            logger.error("❌ Conversion failed")
            return 1


if __name__ == '__main__':
    try:
        exit_code = asyncio.run(main())
        sys.exit(exit_code)
    except KeyboardInterrupt:
        logger.info("Conversion cancelled by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)
