#!/usr/bin/env node

const chokidar = require('chokidar');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { buildStaticSite } = require('./build.js');
const { HTMLToMarkdownConverter } = require('./html-to-markdown.js');

class DevServer {
    constructor() {
        this.isBuilding = false;
        this.buildQueue = false;
        this.liveServerProcess = null;
        this.watcherReady = false;
        this.watcherRestarting = false;
    }

    async startLiveServer() {
        console.log('🚀 Starting live server...');
        
        // Kill any existing live server process
        if (this.liveServerProcess) {
            this.liveServerProcess.kill();
        }

        // Start live-server serving the dist directory
        this.liveServerProcess = spawn('npx', [
            'live-server',
            'dist',
            '--port=55503',
            '--host=localhost',
            '--no-browser',
            '--wait=500'
        ], {
            stdio: 'pipe',
            cwd: __dirname
        });

        this.liveServerProcess.stdout.on('data', (data) => {
            const output = data.toString();
            if (output.includes('Serving') || output.includes('Ready')) {
                console.log('📡 Live server:', output.trim());
            }
        });

        this.liveServerProcess.stderr.on('data', (data) => {
            const error = data.toString();
            if (!error.includes('GET /') && !error.includes('reload')) {
                console.error('🚨 Live server error:', error.trim());
            }
        });
    }

    async build() {
        if (this.isBuilding) {
            this.buildQueue = true;
            return;
        }

        this.isBuilding = true;
        console.log('\n🔄 Rebuilding site...');
        
        try {
            await buildStaticSite();
            console.log('✅ Build complete! Changes will be reflected in the browser.');
            
            // Auto-generate markdown after successful build
            try {
                console.log('📝 Generating markdown...');
                const converter = new HTMLToMarkdownConverter();
                await converter.convertSiteToMarkdown();
                console.log('✅ Markdown generated successfully!\n');
            } catch (markdownError) {
                console.warn('⚠️  Markdown generation failed:', markdownError.message);
                console.log('✅ Build complete! Changes will be reflected in the browser.\n');
            }
            
            // If there was a queued build, run it now
            if (this.buildQueue) {
                this.buildQueue = false;
                setTimeout(() => {
                    this.isBuilding = false;
                    this.build();
                }, 100);
                return;
            }
        } catch (error) {
            console.error('❌ Build failed:', error.message);
        }
        
        this.isBuilding = false;
    }

    async start() {
        console.log('🎯 TEP-GNSS-RINEX Development Server');
        console.log('=====================================\n');

        // Ensure dist directory exists and do initial build
        const distDir = path.join(__dirname, 'dist');
        if (!fs.existsSync(distDir)) {
            fs.mkdirSync(distDir, { recursive: true });
        }

        console.log('🔨 Initial build...');
        await this.build();

        // Start live server
        await this.startLiveServer();

        // Watch for file changes with absolute paths
        const watchPaths = [
            path.join(__dirname, 'components'),
            path.join(__dirname, 'index.html'),
            path.join(__dirname, 'manifest.json'),
            path.join(__dirname, 'figures'),
            path.join(__dirname, 'data'),
            path.join(__dirname, 'public')
        ];
        
        console.log('👁️  Watching paths:', watchPaths);
        
        const watcher = chokidar.watch(watchPaths, {
            ignored: [
                'dist/**',
                'node_modules/**',
                '**/.git/**',
                '**/.DS_Store'
            ],
            persistent: true,
            ignoreInitial: true,
            usePolling: true,
            interval: 1000,
            binaryInterval: 1000,
            atomic: false,
            alwaysStat: true,
            depth: 10,
            awaitWriteFinish: {
                stabilityThreshold: 100,
                pollInterval: 100
            }
        });

        watcher.on('ready', () => {
            if (!this.watcherReady) {
                console.log('👁️  File watcher is ready and monitoring for changes...');
                console.log('👁️  Watched files:', watcher.getWatched());
                this.watcherReady = true;
            }
        });

        watcher.on('change', (filepath) => {
            console.log(`📝 File changed: ${filepath}`);
            this.build();
        });

        watcher.on('add', (filepath) => {
            console.log(`➕ File added: ${filepath}`);
            this.build();
        });

        watcher.on('unlink', (filepath) => {
            console.log(`🗑️  File removed: ${filepath}`);
            this.build();
        });

        watcher.on('error', (error) => {
            console.error('🚨 File watcher error:', error.message);
            if (!this.watcherRestarting) {
                this.watcherRestarting = true;
                console.log('🔄 Attempting to restart file watcher...');
                setTimeout(() => {
                    this.watcherRestarting = false;
                    this.watcherReady = false;
                }, 2000);
            }
        });

        console.log('👀 Watching for changes in:');
        console.log('   • components/*.html');
        console.log('   • index.html');
        console.log('   • manifest.json');
        console.log('   • figures/*.png');
        console.log('   • data/*.json');
        console.log('   • public/*');
        console.log('\n🌐 Server running at: http://localhost:55503');
        console.log('📱 The page will auto-reload when you make changes!');
        console.log('📝 Markdown will be auto-generated after each build!');
        console.log('💡 If auto-reload doesn\'t work, run: npm run build');
        console.log('');

        // Graceful shutdown
        process.on('SIGINT', () => {
            console.log('\n🛑 Shutting down development server...');
            if (watcher) {
                watcher.close();
            }
            if (this.liveServerProcess) {
                this.liveServerProcess.kill();
            }
            process.exit(0);
        });

        process.on('SIGTERM', () => {
            console.log('\n🛑 Shutting down development server...');
            if (watcher) {
                watcher.close();
            }
            if (this.liveServerProcess) {
                this.liveServerProcess.kill();
            }
            process.exit(0);
        });
    }
}

// Run if called directly
if (require.main === module) {
    const devServer = new DevServer();
    devServer.start().catch(console.error);
}

module.exports = DevServer;
