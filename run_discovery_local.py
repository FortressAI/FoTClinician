#!/usr/bin/env python3
"""
FoTChemistry Local Discovery Runner

Runs autonomous discovery campaigns locally without GitHub Actions.
Connects to existing Neo4j instance and exports results for Streamlit.
"""

import argparse
import logging
import time
import glob
import yaml
from pathlib import Path
from orchestrator import AutonomousOrchestrator
from akg.client import AKG

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Run FoTChemistry discovery campaigns locally')
    parser.add_argument('--campaigns', nargs='+', 
                       default=['configs/campaigns/*.yaml'],
                       help='Campaign configuration files or patterns')
    parser.add_argument('--budget', default='1h',
                       choices=['30m', '1h', '2h', '4h', '8h', 'infinite'],
                       help='Time budget for discovery run (use "infinite" for perpetual)')
    parser.add_argument('--continuous', action='store_true',
                       help='Run continuously (restart after budget expires)')
    parser.add_argument('--export-only', action='store_true',
                       help='Only export existing data to Streamlit (no discovery)')
    parser.add_argument('--test-connection', action='store_true',
                       help='Test Neo4j connection and exit')
    
    args = parser.parse_args()
    
    # Convert budget to seconds
    budget_map = {"30m": 1800, "1h": 3600, "2h": 7200, "4h": 14400, "8h": 28800, "infinite": float('inf')}
    budget_seconds = budget_map[args.budget]
    
    logger.info(f"🚀 Starting FoTChemistry Discovery Runner")
    budget_display = "♾️ INFINITE" if args.budget == "infinite" else f"{budget_seconds}s"
    logger.info(f"📅 Budget: {args.budget} ({budget_display})")
    logger.info(f"🎯 Campaigns: {args.campaigns}")
    
    # Initialize AKG
    akg = AKG()
    
    # Test connection
    health = akg.health_check()
    logger.info(f"🏥 AKG Health: Neo4j={health['neo4j']}")
    
    if args.test_connection:
        if health['neo4j']:
            logger.info("✅ Neo4j connection successful!")
            stats = akg.query_discovery_statistics()
            logger.info(f"📊 Current discoveries: {stats['total_discoveries']}")
        else:
            logger.error("❌ Neo4j connection failed!")
        akg.close()
        return
    
    if not health['neo4j']:
        logger.error("❌ Neo4j connection required for discovery. Check credentials.")
        akg.close()
        return
        
    # Export existing data for Streamlit
    logger.info("📁 Exporting current discoveries for Streamlit...")
    export_data = akg.export_for_streamlit()
    logger.info(f"✅ Exported {export_data['total_discoveries']} discoveries")
    
    if args.export_only:
        logger.info("📊 Export-only mode complete")
        akg.close()
        return
    
    # Load campaign configurations
    campaign_files = []
    for pattern in args.campaigns:
        campaign_files.extend(glob.glob(pattern))
    
    if not campaign_files:
        logger.error(f"❌ No campaign files found for patterns: {args.campaigns}")
        akg.close()
        return
        
    campaigns = []
    for file_path in campaign_files:
        try:
            with open(file_path, 'r') as f:
                campaign = yaml.safe_load(f)
                campaigns.append(campaign)
                logger.info(f"📋 Loaded campaign: {campaign.get('name', file_path)}")
        except Exception as e:
            logger.error(f"❌ Failed to load campaign {file_path}: {e}")
    
    if not campaigns:
        logger.error("❌ No valid campaigns loaded")
        akg.close()
        return
    
    # Run discovery
    try:
        if args.continuous:
            logger.info("🔄 Running in continuous mode...")
            cycle = 1
            while True:
                logger.info(f"🔄 Starting discovery cycle {cycle}")
                orchestrator = AutonomousOrchestrator(campaign_files, budget_seconds)
                orchestrator.run_autonomous_discovery()
                
                # Export results after each cycle
                export_data = akg.export_for_streamlit()
                logger.info(f"📁 Cycle {cycle} complete: {export_data['total_discoveries']} total discoveries")
                
                cycle += 1
                logger.info(f"⏱️  Waiting 60s before next cycle...")
                time.sleep(60)
        else:
            logger.info("🎯 Running single discovery session...")
            orchestrator = AutonomousOrchestrator(campaign_files, budget_seconds)
            orchestrator.run_autonomous_discovery()
            
            # Final export
            export_data = akg.export_for_streamlit()
            logger.info(f"✅ Discovery complete: {export_data['total_discoveries']} total discoveries")
            
    except KeyboardInterrupt:
        logger.info("⚠️ Discovery interrupted by user")
    except Exception as e:
        logger.error(f"❌ Discovery failed: {e}")
    finally:
        akg.close()
        logger.info("🏁 Discovery session ended")

if __name__ == "__main__":
    main()
