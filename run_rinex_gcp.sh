#!/bin/bash
# =============================================================================
# TEP-GNSS-RINEX Google Cloud Platform (GCP) Analysis
# OPTIMIZED 100-STATION VERSION for 5-Year TEP Validation
# =============================================================================
#
# Downloads RINEX files for 100 curated stations (2020-2024),
# processes with RTKLIB SPP, and runs TEP orbital coupling analysis.
#
# Station Selection: 100 optimal stations (50 N / 50 S hemisphere)
#   - Data-driven stability selection
#   - Regional clustering for reliable decay fits
#   - Expected significance: 3.2σ with 5 years of data
#
# Author: Matthew Lukin Smawfield
# Date: 9 December 2025
# License: CC-BY-4.0
#
# =============================================================================
# PREREQUISITES: Create GCP Instance
# =============================================================================
#
# Create the compute instance with 200GB data disk:
#    gcloud compute instances create tep-optimal \
#        --zone=us-central1-a \
#        --machine-type=n2d-highcpu-16 \
#        --provisioning-model=SPOT \
#        --instance-termination-action=STOP \
#        --boot-disk-size=50GB \
#        --boot-disk-type=pd-ssd \
#        --image-family=ubuntu-2204-lts \
#        --image-project=ubuntu-os-cloud
#
# Why these settings:
#   - us-central1-a: Good connectivity, low preemption rate
#   - n2d-highcpu-16: 16 AMD vCPUs, 16GB RAM - optimal for parallel processing
#   - SPOT: ~65% cheaper (~$0.12/hr vs $0.34/hr)
#   - 200GB SSD data disk: Sufficient for 100 stations × 5 years
#
# Estimated costs:
#   - Full 5-year run: ~3 days, ~$10-15 (SPOT pricing)
#
# Update .env.local with instance details:
#    GCP_PROJECT_ID=your-project-id
#    GCP_ZONE=us-central1-a
#    GCP_INSTANCE_NAME=tep-optimal
#
# 4. Run this script:
#    ./run_rinex_gcp.sh
#
# When complete, download results:
#    gcloud compute scp --recurse tep-optimal:~/tep-gnss-rinex/results/ ./gcp_results/ --zone=us-central1-a
#    gcloud compute scp --recurse tep-optimal:~/tep-gnss-rinex/data/processed/ ./gcp_processed/ --zone=us-central1-a
#
# Clean up (delete instance):
#    gcloud compute instances delete tep-optimal --zone=us-central1-a --quiet
#
# =============================================================================

set -e

# GCP Configuration - Load from environment variables
if [ -f ".env.local" ]; then
    export $(grep -v '^#' .env.local | xargs)
elif [ -f ".env" ]; then
    export $(grep -v '^#' .env | xargs)
fi

# Get configuration from environment variables
PROJECT_ID="${GCP_PROJECT_ID:-}"
ZONE="${GCP_ZONE:-}"
INSTANCE_NAME="${GCP_INSTANCE_NAME:-}"
PACKAGE_NAME="tep-gnss-rinex-gcp.tar.gz"

# CDDIS credentials (required for RINEX downloads)
CDDIS_USER="${CDDIS_USER:-}"
CDDIS_PASS="${CDDIS_PASS:-}"

echo "🚀 TEP-GNSS-RINEX Google Cloud Platform Analysis"
echo "================================================="

# Check if gcloud CLI is available
if ! command -v gcloud &> /dev/null; then
    echo "❌ gcloud CLI not found. Please install it first:"
    echo "   https://cloud.google.com/sdk/docs/install"
    exit 1
fi

# Check GCP configuration
if [ -z "$PROJECT_ID" ] || [ -z "$ZONE" ] || [ -z "$INSTANCE_NAME" ]; then
    echo "❌ GCP configuration missing!"
    echo ""
    echo "Please set the following environment variables:"
    echo "   GCP_PROJECT_ID: Your GCP project ID"
    echo "   GCP_ZONE: Your GCP zone (e.g., us-central1-a)"
    echo "   GCP_INSTANCE_NAME: Your GCP instance name"
    echo ""
    echo "Options:"
    echo "1. Create .env.local file:"
    echo "   GCP_PROJECT_ID=your-project-id"
    echo "   GCP_ZONE=us-central1-a"
    echo "   GCP_INSTANCE_NAME=your-instance-name"
    echo "   CDDIS_USER=your-cddis-username"
    echo "   CDDIS_PASS=your-cddis-password"
    echo ""
    echo "2. Set environment variables directly:"
    echo "   export GCP_PROJECT_ID=your-project-id"
    echo "   export GCP_ZONE=us-central1-a"
    echo "   export GCP_INSTANCE_NAME=your-instance-name"
    echo "   export CDDIS_USER=your-cddis-username"
    echo "   export CDDIS_PASS=your-cddis-password"
    exit 1
fi

# Check CDDIS credentials
if [ -z "$CDDIS_USER" ] || [ -z "$CDDIS_PASS" ]; then
    echo "❌ CDDIS credentials missing!"
    echo ""
    echo "RINEX downloads require CDDIS authentication."
    echo "Register at: https://urs.earthdata.nasa.gov/"
    echo ""
    echo "Set credentials:"
    echo "   export CDDIS_USER=your-username"
    echo "   export CDDIS_PASS=your-password"
    exit 1
fi

echo "✅ GCP Configuration:"
echo "   Project: $PROJECT_ID"
echo "   Zone: $ZONE"
echo "   Instance: $INSTANCE_NAME"
echo "   CDDIS User: $CDDIS_USER"

# Set the project
gcloud config set project $PROJECT_ID

# Get instance external IP
EXTERNAL_IP=$(gcloud compute instances describe $INSTANCE_NAME \
    --zone=$ZONE \
    --format='get(networkInterfaces[0].accessConfigs[0].natIP)' 2>/dev/null || echo "")

if [ -z "$EXTERNAL_IP" ]; then
    echo "⚠️  Instance may not be running or doesn't exist yet."
    echo "   Create instance with (SPOT instance recommended):"
    echo "   gcloud compute instances create $INSTANCE_NAME \\"
    echo "       --zone=$ZONE \\"
    echo "       --machine-type=t2d-standard-8 \\"
    echo "       --provisioning-model=SPOT \\"
    echo "       --boot-disk-size=50GB \\"
    echo "       --boot-disk-type=pd-balanced \\"
    echo "       --image-family=ubuntu-2204-lts \\"
    echo "       --image-project=ubuntu-os-cloud"
    exit 1
fi

echo "   External IP: $EXTERNAL_IP"

# Create analysis package
echo "📦 Creating TEP-GNSS-RINEX analysis package..."
rm -f $PACKAGE_NAME

# Package the RINEX project (exclude large data files and compiled binaries)
# RTKLIB will be recompiled on Linux - only need source code
# Include optimal_100_metadata.json for station selection
tar --exclude='*.pyc' --exclude='__pycache__' --exclude='.git' \
    --exclude='*.log' --exclude='*.pid' --exclude='*.npz' \
    --exclude='data/rinex/*' --exclude='data/nav/*' --exclude='data/sp3/*' \
    --exclude='results/*' --exclude='logs/*' \
    --exclude='venv' --exclude='*.tar.gz' --exclude='*.zip' \
    --exclude='*.o' --exclude='*.a' --exclude='*.so' \
    --exclude='RTKLIB/.git' \
    --exclude='RTKLIB/test' \
    --exclude='RTKLIB/data' \
    --exclude='RTKLIB/doc' \
    --exclude='RTKLIB/util' \
    --exclude='RTKLIB/bin' \
    --exclude='RTKLIB/lib/iers' \
    --exclude='RTKLIB/app/winapp' \
    --exclude='RTKLIB/app/qtapp' \
    -czf $PACKAGE_NAME \
    scripts/ RTKLIB/src/ RTKLIB/app/consapp/rnx2rtkp/ \
    run_full_analysis.py \
    data/processed/optimal_100_metadata.json \
    data/processed/station_coordinates.json 2>/dev/null || \
tar -czf $PACKAGE_NAME scripts/ RTKLIB/src/ RTKLIB/app/consapp/rnx2rtkp/ \
    run_full_analysis.py \
    data/processed/optimal_100_metadata.json data/processed/station_coordinates.json

PACKAGE_SIZE=$(du -h $PACKAGE_NAME | cut -f1)
echo "   Package size: $PACKAGE_SIZE"

# Test SSH connection
echo "🔗 Testing SSH connection..."
if ! gcloud compute ssh $INSTANCE_NAME --zone=$ZONE --command="echo 'SSH OK'" >/dev/null 2>&1; then
    echo "❌ SSH connection failed!"
    echo "   Please check:"
    echo "   1. Instance is running"
    echo "   2. Firewall rules allow SSH (port 22)"
    echo "   3. You have the correct permissions"
    exit 1
fi

echo "✅ SSH connection established"

# Transfer package
echo "📤 Transferring analysis package to GCP..."
gcloud compute scp $PACKAGE_NAME $INSTANCE_NAME:~/ --zone=$ZONE

# Create the setup script that will run on GCP
echo "🚀 Creating GCP setup and execution script..."
cat > /tmp/gcp_run_rinex.sh << REMOTE_SCRIPT_EOF
#!/bin/bash
set -e

echo "🔧 Setting up GCP environment for TEP-GNSS-RINEX analysis..."

# CDDIS credentials (passed from local machine)
export CDDIS_USER="$CDDIS_USER"
export CDDIS_PASS="$CDDIS_PASS"

# Create .netrc for CDDIS authentication
echo "🔐 Setting up CDDIS authentication..."
cat > ~/.netrc << NETRC_EOF
machine urs.earthdata.nasa.gov
    login $CDDIS_USER
    password $CDDIS_PASS
NETRC_EOF
chmod 600 ~/.netrc
echo "   ✅ .netrc configured"

# Update system and install dependencies
echo "📦 Installing system dependencies..."
sudo apt update >/dev/null 2>&1
sudo apt install -y python3 python3-pip python3-venv python3-dev \
    build-essential gfortran gcc g++ make cmake \
    libhdf5-dev libnetcdf-dev htop iotop wget curl \
    >/dev/null 2>&1

# Get system specs
CORES=\$(nproc)
MEMORY_GB=\$(free -g | awk "/^Mem:/{print int(\\\$2*0.8)}")

echo ""
echo "⚡ GCP System Configuration:"
echo "  CPU cores: \$CORES"
echo "  Memory: \${MEMORY_GB}GB"
echo "  Analysis: OPTIMIZED 100 stations × 5 years (2020-2024)"
echo "  Expected: 182,600 station-days → 3.2σ TEP signal"

# Setup working directory
WORK_DIR=\$HOME/tep-gnss-rinex
mkdir -p \$WORK_DIR

# MOUNT ATTACHED SSD (Critical for storage)
# Check multiple possible device paths for attached disks
echo "💾 Looking for attached data disk..."
DISK_DEVICE=""
for dev in /dev/sdb /dev/nvme0n2 /dev/disk/by-id/google-data-disk; do
    if [ -b "\$dev" ] || [ -L "\$dev" ]; then
        DISK_DEVICE="\$dev"
        echo "   Found disk at: \$DISK_DEVICE"
        break
    fi
done

if [ -n "\$DISK_DEVICE" ]; then
    # Check if already mounted
    if mount | grep -q "\$WORK_DIR"; then
        echo "   ✅ Disk already mounted to \$WORK_DIR"
    else
        # Check if formatted
        if ! sudo blkid "\$DISK_DEVICE" >/dev/null 2>&1; then
            echo "   Formatting disk..."
            sudo mkfs.ext4 -m 0 -E lazy_itable_init=0,lazy_journal_init=0,discard "\$DISK_DEVICE"
        fi
        
        # Mount it
        echo "   Mounting disk to \$WORK_DIR..."
        sudo mount -o discard,defaults "\$DISK_DEVICE" "\$WORK_DIR"
        sudo chown \$USER:\$USER "\$WORK_DIR"
        echo "   ✅ 200GB SSD mounted successfully"
    fi
else
    echo "⚠️  No secondary disk found. Using boot disk (watch space!)."
fi

cd \$WORK_DIR
echo "  Working directory: \$(pwd)"

# Extract package
echo "📦 Extracting analysis package..."
cp ~/tep-gnss-rinex-gcp.tar.gz ./
tar -xzf tep-gnss-rinex-gcp.tar.gz

# Create required directories
mkdir -p data/{rinex,nav,sp3,processed}
mkdir -p logs results/{outputs,figures}

# Build RTKLIB - ALWAYS recompile for Linux (Mac binaries won't work)
echo "🔧 Building RTKLIB for Linux..."
cd RTKLIB/app/consapp/rnx2rtkp/gcc

# Check if binary works on this platform
if ./rnx2rtkp --help >/dev/null 2>&1; then
    echo "  ✅ RTKLIB already built for this platform"
else
    echo "  Compiling RTKLIB from source..."
    make clean 2>/dev/null || true
    make
    echo "  ✅ RTKLIB compiled successfully"
fi

cd \$WORK_DIR

# Verify RTKLIB
if [ -x "RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp" ]; then
    echo "  ✅ RTKLIB executable ready"
else
    echo "  ❌ RTKLIB build failed!"
    exit 1
fi

# Setup Python virtual environment
echo "🐍 Setting up Python environment..."
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip setuptools wheel >/dev/null 2>&1
pip install numpy scipy requests tqdm pandas matplotlib hatanaka astropy jplephem unlzw3 >/dev/null 2>&1
echo "  ✅ Python environment ready"

# Set environment variables for processing
export TEP_WORKERS=\$CORES
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=\$CORES

echo ""
echo "🚀 Starting TEP-GNSS-RINEX OPTIMIZED analysis..."
echo "   Date Range: 2020-01-01 to 2024-12-31 (5 years)"
echo "   Stations: 100 optimal (50 N / 50 S balanced)"
echo "   Workers: \$CORES parallel workers"
echo "   CDDIS User: \$CDDIS_USER"
echo ""
echo "   Pipeline:"
echo "     1. Download RINEX for 100 curated stations"
echo "     2. Download broadcast navigation files"
echo "     3. Process with RTKLIB SPP"
echo "     4. Extract clock bias and position jitter"
echo "     5. Run Step 2.0 coherence analysis"
echo "     6. Run Step 2.2 orbital coupling analysis"
echo ""
echo "   Start time: \$(date)"
echo ""

# Run the full analysis pipeline
echo "📡 Starting TEP-GNSS-RINEX full pipeline..."
source venv/bin/activate

# Create a wrapper script that runs the full TEP pipeline using run_full_analysis.py
cat > run_full_pipeline.sh << 'PIPELINE_EOF'
#!/bin/bash
set -e
cd \$HOME/tep-gnss-rinex
source venv/bin/activate

echo "[\$(date)] =========================================="
echo "TEP-GNSS-RINEX FULL PIPELINE (All Steps)"
echo "=========================================="
echo ""

# Run the complete pipeline orchestrator with all filters
python run_full_analysis.py --filters optimal_100 all_stations

echo ""
echo "[\$(date)] ===== PIPELINE COMPLETE ====="
echo ""
echo "Results saved to:"
echo "  - results/outputs/step_*.json"
echo "  - results/figures/step_*.png"
echo "  - logs/*.log"
echo ""
echo "All 11 analysis steps completed:"
echo "  1.0  Data Acquisition"
echo "  1.1  Dynamic 50 Metadata"
echo "  2.0  Raw SPP Analysis"
echo "  2.1  Control Tests + Elevation"
echo "  2.2  Anisotropy Analysis"
echo "  2.3  Temporal + Kp Stratification"
echo "  2.4  Diurnal + Null Tests"
echo "  2.5  Orbital Coupling"
echo "  2.6  Planetary Events"
echo "  2.7  CMB Frame Analysis"
PIPELINE_EOF
chmod +x run_full_pipeline.sh

# Run the pipeline in background
nohup ./run_full_pipeline.sh > rinex_processing.log 2>&1 &
PROC_PID=\$!
echo "  ✅ Full pipeline started with PID: \$PROC_PID"
echo "  Monitor: tail -f \$WORK_DIR/rinex_processing.log"

# Save PID for later
echo \$PROC_PID > processing.pid

echo ""
echo "✅ TEP-GNSS-RINEX full pipeline is now running on GCP!"
echo "   Step 1: Download & Process RINEX files"
echo "   Step 2: Run correlation analysis (automatically after Step 1)"
echo "   You can disconnect safely - the process will keep running."
echo ""
echo "📊 Useful commands:"
echo "   Monitor progress: tail -f \$WORK_DIR/rinex_processing.log"
echo "   Check status: ps aux | grep \$PROC_PID"
echo "   Count processed: ls \$WORK_DIR/data/processed/*.npz 2>/dev/null | wc -l"
REMOTE_SCRIPT_EOF

# Transfer the script to GCP
echo "📤 Transferring setup script to GCP..."
gcloud compute scp /tmp/gcp_run_rinex.sh $INSTANCE_NAME:~/ --zone=$ZONE

# Execute the script on GCP - RUN INTERACTIVELY for real-time output
echo "▶️  Executing setup and analysis script on GCP..."
echo "   (You will see real-time output from the remote machine)"
echo ""
gcloud compute ssh $INSTANCE_NAME --zone=$ZONE --command="chmod +x ~/gcp_run_rinex.sh && ~/gcp_run_rinex.sh"

echo ""
echo "✅ GCP RINEX analysis setup initiated!"
echo ""
echo "📊 Monitor setup progress:"
echo "gcloud compute ssh $INSTANCE_NAME --zone=$ZONE --command='tail -f ~/gcp_setup.log'"
echo ""
echo "📊 Monitor processing progress (after setup completes):"
echo "gcloud compute ssh $INSTANCE_NAME --zone=$ZONE --command='tail -f ~/tep-gnss-rinex/rinex_processing.log'"
echo ""
echo "📊 Connect to instance:"
echo "gcloud compute ssh $INSTANCE_NAME --zone=$ZONE"
echo ""
echo "📊 Check processed file count:"
echo "gcloud compute ssh $INSTANCE_NAME --zone=$ZONE --command='ls ~/tep-gnss-rinex/data/processed/*.npz 2>/dev/null | wc -l'"
echo ""
echo "📥 Download results when complete:"
echo "gcloud compute scp --recurse $INSTANCE_NAME:~/tep-gnss-rinex/data/processed/ ./gcp_processed/ --zone=$ZONE"
echo "gcloud compute scp --recurse $INSTANCE_NAME:~/tep-gnss-rinex/results/ ./gcp_results/ --zone=$ZONE"
echo ""
echo "💰 To stop instance (saves money):"
echo "gcloud compute instances stop $INSTANCE_NAME --zone=$ZONE"
echo ""
echo "🛑 To delete instance:"
echo "gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE"
echo ""
echo "Instance: $INSTANCE_NAME"
echo "External IP: $EXTERNAL_IP"
echo ""
echo "💡 Optimized Instance for 100-Station TEP Analysis:"
echo "   - n2d-highcpu-16: 16 AMD vCPUs, 16 GB RAM (RECOMMENDED)"
echo "   - SPOT pricing: ~\$0.12/hr (65% savings)"
echo "   - Estimated runtime: ~3 days for 5-year dataset"
echo "   - Estimated total cost: ~\$10-15"
echo ""
echo "📊 Expected Results:"
echo "   - 182,600 processed station-days"
echo "   - Orbital correlation: ~3.2σ significance"
echo "   - E-W/N-S anisotropy with orbital tracking"
