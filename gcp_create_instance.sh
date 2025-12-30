#!/bin/bash
# =============================================================================
# TEP-GNSS-RINEX: Create Optimized GCP Instance
# =============================================================================
#
# Creates a cost-effective, high-performance GCP instance optimized for
# CDDIS data downloads and RINEX processing.
#
# Key optimizations:
#   - us-east4: Closest region to NASA Goddard/CDDIS (~40 miles)
#   - SPOT: 60-80% cheaper than on-demand
#   - SSD: Fast I/O for RINEX processing
#   - Premium Network: Uses Google's backbone for faster downloads
#
# Estimated cost: ~$0.48/hour (SPOT)
#
# Author: Matthew Lukin Smawfield
# Date: 9 December 2025
# License: CC-BY-4.0
#
# =============================================================================

set -e

# Configuration
INSTANCE_NAME="${1:-tep-rinex-fast}"
ZONE="us-east4-c"
MACHINE_TYPE="n1-highcpu-64"  # 64 vCPU, 57.6 GB RAM

echo "🚀 Creating optimized GCP instance for TEP-GNSS-RINEX..."
echo "   Instance: $INSTANCE_NAME"
echo "   Zone: $ZONE (closest to CDDIS/NASA Goddard)"
echo "   Machine: $MACHINE_TYPE"
echo "   Provisioning: SPOT (~80% cheaper)"
echo ""

gcloud compute instances create "$INSTANCE_NAME" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --provisioning-model=SPOT \
    --boot-disk-size=20GB \
    --boot-disk-type=pd-ssd \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --network-tier=PREMIUM \
    --create-disk=name="${INSTANCE_NAME}-data",size=200GB,type=pd-ssd,auto-delete=no

echo ""
echo "✅ Instance created successfully!"
echo ""
echo "📝 Update your .env.local:"
echo "   GCP_ZONE=$ZONE"
echo "   GCP_INSTANCE_NAME=$INSTANCE_NAME"
echo ""
echo "🚀 Then run: ./run_rinex_gcp.sh"
echo ""
echo "💰 Cost: ~\$0.48/hour (SPOT pricing)"
echo ""
echo "⚠️  Data disk has auto-delete=NO - will persist after instance deletion"
echo ""
echo "🛑 When done, stop or delete to save money:"
echo "   gcloud compute instances stop $INSTANCE_NAME --zone=$ZONE"
echo "   gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE"
