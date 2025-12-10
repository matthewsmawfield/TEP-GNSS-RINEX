#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

/**
 * HTML to Markdown Converter for TEP-GNSS-RINEX Site
 * Converts the built static HTML site into a clean markdown document
 */

class HTMLToMarkdownConverter {
    constructor() {
        this.output = '';
        this.currentSection = '';
    }

    /**
     * Convert HTML string to markdown with proper academic formatting
     */
    htmlToMarkdown(html) {
        // Remove script tags and their content
        html = html.replace(/<script\b[^<]*(?:(?!<\/script>)<[^<]*)*<\/script>/gi, '');
        
        // Remove style tags and their content
        html = html.replace(/<style\b[^<]*(?:(?!<\/style>)<[^<]*)*<\/style>/gi, '');
        
        // Remove comments
        html = html.replace(/<!--[\s\S]*?-->/g, '');
        
        // Preserve MathJax expressions before processing
        const mathExpressions = [];
        html = html.replace(/<span[^>]*class=["'][^"']*MathJax[^"']*["'][^>]*>.*?<\/span>/gi, (match) => {
            mathExpressions.push(match);
            return `__MATH_EXPRESSION_${mathExpressions.length - 1}__`;
        });
        
        // Convert manuscript sections to proper markdown structure FIRST
        html = html.replace(/<div[^>]*class=["'][^"']*manuscript-section[^"']*["'][^>]*data-section=["']([^"']*)["'][^>]*>/gi, '\n\n## $1\n\n');
        
        // Convert headers
        html = html.replace(/<h1[^>]*>(.*?)<\/h1>/gi, '\n# $1\n\n');
        html = html.replace(/<h2[^>]*>(.*?)<\/h2>/gi, '\n## $1\n\n');
        html = html.replace(/<h3[^>]*>(.*?)<\/h3>/gi, '\n### $1\n\n');
        html = html.replace(/<h4[^>]*>(.*?)<\/h4>/gi, '\n#### $1\n\n');
        
        // Convert paragraphs
        html = html.replace(/<p[^>]*>(.*?)<\/p>/gi, '$1\n\n');
        
        // Convert strong/bold
        html = html.replace(/<(strong|b)[^>]*>(.*?)<\/(strong|b)>/gi, '**$2**');
        
        // Convert emphasis/italic
        html = html.replace(/<(em|i)[^>]*>(.*?)<\/(em|i)>/gi, '*$2*');
        
        // Convert links
        html = html.replace(/<a[^>]*href=["']([^"']*)["'][^>]*>(.*?)<\/a>/gi, '[$2]($1)');
        
        // Convert lists
        html = html.replace(/<ul[^>]*>/gi, '\n');
        html = html.replace(/<\/ul>/gi, '\n');
        html = html.replace(/<ol[^>]*>/gi, '\n');
        html = html.replace(/<\/ol>/gi, '\n');
        html = html.replace(/<li[^>]*>(.*?)<\/li>/gi, '- $1\n');
        
        // Convert blockquotes
        html = html.replace(/<blockquote[^>]*>(.*?)<\/blockquote>/gi, '\n> $1\n\n');
        
        // Convert code blocks
        html = html.replace(/<pre[^>]*><code[^>]*>(.*?)<\/code><\/pre>/gi, '\n```\n$1\n```\n\n');
        html = html.replace(/<code[^>]*>(.*?)<\/code>/gi, '`$1`');
        
        // Convert line breaks
        html = html.replace(/<br\s*\/?>/gi, '\n');
        
        // Convert horizontal rules
        html = html.replace(/<hr\s*\/?>/gi, '\n---\n\n');
        
        // Convert divs with special classes to markdown equivalents
        html = html.replace(/<div[^>]*class=["'][^"']*abstract[^"']*["'][^>]*>/gi, '');
        html = html.replace(/<div[^>]*class=["'][^"']*theorem[^"']*["'][^>]*>/gi, '\n**Theorem:**\n');
        html = html.replace(/<div[^>]*class=["'][^"']*principle[^"']*["'][^>]*>/gi, '\n**Principle:**\n');
        html = html.replace(/<div[^>]*class=["'][^"']*proof[^"']*["'][^>]*>/gi, '\n*Proof:*\n');
        html = html.replace(/<div[^>]*class=["'][^"']*experimental-section[^"']*["'][^>]*>/gi, '\n**Experimental Section:**\n');
        html = html.replace(/<div[^>]*class=["'][^"']*critical-analysis[^"']*["'][^>]*>/gi, '\n**Critical Analysis:**\n');
        html = html.replace(/<div[^>]*class=["'][^"']*significance[^"']*["'][^>]*>/gi, '\n**Significance:**\n');
        html = html.replace(/<div[^>]*class=["'][^"']*key-finding[^"']*["'][^>]*>/gi, '\n**Key Finding:**\n');
        html = html.replace(/<div[^>]*class=["'][^"']*validation-box[^"']*["'][^>]*>/gi, '\n**Validation:**\n');
        
        // Handle abstract section specially
        html = html.replace(/## Abstract\s*\n\s*<h2>Abstract<\/h2>/gi, '## Abstract\n\n');
        
        // Convert tables
        html = html.replace(/<table[^>]*>(.*?)<\/table>/gis, (match) => {
            return this.convertTable(match);
        });
        
        // Remove remaining HTML tags
        html = html.replace(/<[^>]+>/g, '');
        
        // Restore MathJax expressions
        mathExpressions.forEach((expr, index) => {
            html = html.replace(`__MATH_EXPRESSION_${index}__`, expr);
        });
        
        // Decode HTML entities
        html = html.replace(/&amp;/g, '&');
        html = html.replace(/&lt;/g, '<');
        html = html.replace(/&gt;/g, '>');
        html = html.replace(/&quot;/g, '"');
        html = html.replace(/&#39;/g, "'");
        html = html.replace(/&nbsp;/g, ' ');
        html = html.replace(/&times;/g, '×');
        html = html.replace(/&minus;/g, '−');
        html = html.replace(/&plusmn;/g, '±');
        html = html.replace(/&sup2;/g, '²');
        html = html.replace(/&sup3;/g, '³');
        html = html.replace(/&sup1;/g, '¹');
        html = html.replace(/&deg;/g, '°');
        html = html.replace(/&lambda;/g, 'λ');
        html = html.replace(/&mu;/g, 'μ');
        html = html.replace(/&sigma;/g, 'σ');
        
        // Clean up whitespace
        html = html.replace(/\n\s*\n\s*\n/g, '\n\n');
        html = html.replace(/^\s+|\s+$/g, '');
        
        // Remove duplicate headers (same header appearing consecutively)
        html = html.replace(/(##\s+[^\n]+)\n+\1/g, '$1');
        
        // Clean up any remaining formatting issues
        html = html.replace(/\n{3,}/g, '\n\n');
        
        return html;
    }

    /**
     * Convert HTML table to markdown table
     */
    convertTable(tableHtml) {
        const rows = [];
        const rowMatches = tableHtml.match(/<tr[^>]*>(.*?)<\/tr>/gis);
        
        if (!rowMatches) return '';
        
        rowMatches.forEach((row, index) => {
            const cells = row.match(/<t[dh][^>]*>(.*?)<\/t[dh]>/gis);
            if (cells) {
                const cellTexts = cells.map(cell => 
                    cell.replace(/<[^>]+>/g, '').trim()
                );
                rows.push(cellTexts);
            }
        });
        
        if (rows.length === 0) return '';
        
        // Create markdown table
        let markdown = '\n';
        rows.forEach((row, index) => {
            markdown += '| ' + row.join(' | ') + ' |\n';
            if (index === 0) {
                // Add separator row
                markdown += '|' + row.map(() => ' --- ').join('|') + '|\n';
            }
        });
        markdown += '\n';
        
        return markdown;
    }

    /**
     * Extract title and metadata from HTML
     */
    extractMetadata(html) {
        const titleMatch = html.match(/<title[^>]*>(.*?)<\/title>/i);
        const title = titleMatch ? titleMatch[1] : 'Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks';
        
        const authorMatch = html.match(/<meta[^>]*name=["']author["'][^>]*content=["']([^"']*)["']/i);
        const author = authorMatch ? authorMatch[1] : 'Matthew Lukin Smawfield';
        
        const versionMatch = html.match(/<div[^>]*class=["'][^"']*version[^"']*["'][^>]*>(.*?)<\/div>/i);
        const version = versionMatch ? versionMatch[1]
            .replace(/<[^>]+>/g, '')
            .replace(/^Version:\s*/i, '')
            .trim() : 'v0.1 (Istanbul)';
        
        const dateMatch = html.match(/<div[^>]*class=["'][^"']*date[^"']*["'][^>]*>(.*?)<\/div>/i);
        const date = dateMatch ? dateMatch[1].replace(/<[^>]+>/g, '').trim() : 'First published: 29 November 2025';
        
        const doiMatch = html.match(/DOI:\s*<a[^>]*href=["']([^"']*)["'][^>]*>(.*?)<\/a>/i);
        const doi = doiMatch ? doiMatch[2] : '10.5281/zenodo.17860166';
        
        return { title, author, version, date, doi };
    }

    /**
     * Extract main content from HTML
     */
    extractMainContent(html) {
        // Find the manuscript-content div and extract everything until the closing main tag
        const startMatch = html.match(/<div[^>]*id=["']manuscript-content["'][^>]*>/i);
        if (!startMatch) {
            throw new Error('Could not find manuscript-content div');
        }
        
        const startIndex = startMatch.index + startMatch[0].length;
        
        // Find the closing main tag
        const endMatch = html.match(/<\/main>/i);
        if (!endMatch) {
            throw new Error('Could not find closing main tag');
        }
        
        const endIndex = endMatch.index;
        
        // Extract the content between these points
        const content = html.substring(startIndex, endIndex);
        
        return content;
    }

    /**
     * Convert the built HTML site to markdown
     */
    async convertSiteToMarkdown() {
        console.log('🔄 Converting HTML site to markdown...');
        
        try {
            // Read the built HTML file
            const htmlPath = path.join(__dirname, 'dist', 'index.html');
            if (!fs.existsSync(htmlPath)) {
                throw new Error('Built HTML file not found. Please run "npm run build" first.');
            }
            
            const html = fs.readFileSync(htmlPath, 'utf8');
            
            // Extract metadata
            const metadata = this.extractMetadata(html);
            
            // Extract main content
            const mainContent = this.extractMainContent(html);
            
            // Convert to markdown
            const markdownContent = this.htmlToMarkdown(mainContent);
            
            // Build the complete markdown document
            const markdown = this.buildMarkdownDocument(metadata, markdownContent);
            
            // Write to file
            const outputPath = path.join(__dirname, '..', 'manuscript-rinex.md');
            fs.writeFileSync(outputPath, markdown, 'utf8');
            
            console.log('✅ Markdown conversion complete!');
            console.log(`📄 Output: ${outputPath}`);
            console.log(`📊 Document: ${metadata.title}`);
            console.log(`👤 Author: ${metadata.author}`);
            console.log(`📅 Version: ${metadata.version}`);
            
            return outputPath;
            
        } catch (error) {
            console.error('❌ Markdown conversion failed:', error.message);
            process.exit(1);
        }
    }

    /**
     * Build the complete markdown document with metadata
     */
    buildMarkdownDocument(metadata, content) {
        const timestamp = new Date().toISOString().split('T')[0];
        
        // Clean up the title to remove the author part
        const cleanTitle = metadata.title.replace(' | Matthew Lukin Smawfield', '');
        
        return `# ${cleanTitle}

**Author:** ${metadata.author}  
**Version:** ${metadata.version}  
**Date:** ${metadata.date}  
**DOI:** ${metadata.doi}  
**Generated:** ${timestamp}  
**Paper Series:** TEP-GNSS Paper 3 (Raw RINEX Validation)

---

${content}

---

*This document was automatically generated from the TEP-GNSS-RINEX research site. For the interactive version with figures and enhanced formatting, visit: https://matthewsmawfield.github.io/TEP-GNSS-RINEX/*

*Related Papers:*
- *Paper 1 (Multi-Center Validation): https://matthewsmawfield.github.io/TEP-GNSS/*
- *Paper 2 (25-Year CODE Analysis): https://matthewsmawfield.github.io/TEP-GNSS/code-longspan/*

*Source code and data available at: https://github.com/matthewsmawfield/TEP-GNSS-RINEX*
`;
    }
}

// Main execution
async function main() {
    const converter = new HTMLToMarkdownConverter();
    await converter.convertSiteToMarkdown();
}

// Run if called directly
if (require.main === module) {
    main();
}

module.exports = { HTMLToMarkdownConverter };
