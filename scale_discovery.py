#!/usr/bin/env python3
"""
FoTChemistry Discovery Scaling System

This script demonstrates how to scale molecular discovery to find thousands
of new molecules using quantum-guided search and parallel processing.
"""

import argparse
import multiprocessing as mp
import logging
from typing import List, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import os
from datetime import datetime

from continuous_chemistry_discovery import ContinuousChemistryDiscoveryEngine, ChemistryDiscoveryConfig

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('scaling_discovery.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ScaledDiscoveryConfig:
    """Configuration for scaled molecular discovery"""
    
    def __init__(self):
        # Scaling parameters
        self.total_target_molecules = 10000  # Target: 10K new molecules
        self.parallel_processes = min(8, mp.cpu_count() - 1)  # Leave 1 CPU free
        self.batch_size_per_process = 20  # Molecules per batch per process
        self.max_concurrent_batches = 4  # Limit concurrent batches to avoid memory issues
        
        # Discovery objectives
        self.objectives = [
            "drug_discovery",
            "green_chemistry", 
            "materials_science",
            "catalysis",
            "energy_storage"
        ]
        
        # Seed molecule libraries
        self.seed_libraries = {
            "pharmaceuticals": [
                "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
                "CC(=O)NC1=CC=C(C=C1)O",  # Paracetamol
                "CN(C)CCN(C1=CC=CC=C1)C2=CC=CC=C2",  # Diphenhydramine
            ],
            "natural_products": [
                "CC1=C(C(=O)C2=C(C1=O)C=CC=C2O)C",  # Vitamin K
                "CC(=O)OC1=CC=CC=C1C(=O)O",  # Salicylic acid derivative
                "C1=CC(=C(C=C1C=CC(=O)O)O)O",  # Caffeic acid
                "COC1=CC(=CC(=C1O)OC)C=CC(=O)O",  # Ferulic acid
            ],
            "materials": [
                "C1=CC=C2C(=C1)C=CC=C2",  # Naphthalene
                "C1=CC=C(C=C1)C2=CC=CC=C2",  # Biphenyl
                "C1=CC=C(C=C1)C=CC2=CC=CC=C2",  # Stilbene
                "C1=CC=C2C(=C1)SC=C2",  # Benzothiophene
            ],
            "green_solvents": [
                "CCO",  # Ethanol
                "CC(C)O",  # Isopropanol
                "CCOCCO",  # Ethylene glycol monomethyl ether
                "CC(=O)C",  # Acetone (for comparison)
            ]
        }
        
        # Quality thresholds
        self.min_combined_score = 0.3  # Lower threshold for scaled discovery
        self.target_drug_likeness = 0.6
        self.target_safety_score = 0.7

class ParallelDiscoveryOrchestrator:
    """Orchestrates parallel molecular discovery across multiple processes"""
    
    def __init__(self, config: ScaledDiscoveryConfig):
        self.config = config
        self.discovered_molecules = set()  # Track unique SMILES
        self.total_discoveries = 0
        self.start_time = None
        
    def run_scaled_discovery(self) -> Dict[str, Any]:
        """Run scaled molecular discovery campaign"""
        self.start_time = time.time()
        logger.info(f"🚀 Starting scaled discovery campaign")
        logger.info(f"🎯 Target: {self.config.total_target_molecules} unique molecules")
        logger.info(f"⚡ Parallel processes: {self.config.parallel_processes}")
        logger.info(f"📦 Batch size per process: {self.config.batch_size_per_process}")
        
        all_discoveries = []
        
        try:
            with ProcessPoolExecutor(max_workers=self.config.parallel_processes) as executor:
                # Submit discovery batches
                future_to_config = {}
                
                batches_submitted = 0
                while self.total_discoveries < self.config.total_target_molecules:
                    # Limit concurrent batches to avoid memory issues
                    if len(future_to_config) >= self.config.max_concurrent_batches:
                        # Wait for at least one batch to complete
                        completed = next(as_completed(future_to_config))
                        batch_result = completed.result()
                        all_discoveries.extend(batch_result['discoveries'])
                        self._update_progress(batch_result)
                        del future_to_config[completed]
                    
                    # Submit new batch
                    batch_config = self._create_batch_config(batches_submitted)
                    future = executor.submit(self._run_discovery_batch, batch_config)
                    future_to_config[future] = batch_config
                    batches_submitted += 1
                    
                    # Progress update
                    if batches_submitted % 10 == 0:
                        self._log_progress()
                
                # Wait for remaining batches
                for future in as_completed(future_to_config):
                    batch_result = future.result()
                    all_discoveries.extend(batch_result['discoveries'])
                    self._update_progress(batch_result)
                    
        except KeyboardInterrupt:
            logger.info("🛑 Discovery campaign interrupted by user")
        except Exception as e:
            logger.error(f"❌ Discovery campaign failed: {e}")
            
        # Final statistics
        end_time = time.time()
        duration = end_time - self.start_time
        
        results = {
            'total_discoveries': len(all_discoveries),
            'unique_molecules': len(self.discovered_molecules),
            'duration_seconds': duration,
            'discoveries_per_hour': len(all_discoveries) / (duration / 3600),
            'discoveries': all_discoveries,
            'campaign_timestamp': datetime.now().isoformat()
        }
        
        self._log_final_results(results)
        return results
    
    def _create_batch_config(self, batch_number: int) -> ChemistryDiscoveryConfig:
        """Create configuration for a discovery batch"""
        # Rotate through objectives and seed libraries
        objective = self.config.objectives[batch_number % len(self.config.objectives)]
        
        # Select seed library based on objective
        if objective == "drug_discovery":
            seeds = self.config.seed_libraries["pharmaceuticals"]
        elif objective == "green_chemistry":
            seeds = self.config.seed_libraries["green_solvents"]
        elif objective in ["materials_science", "energy_storage"]:
            seeds = self.config.seed_libraries["materials"]
        else:
            seeds = self.config.seed_libraries["natural_products"]
        
        return ChemistryDiscoveryConfig(
            batch_size=self.config.batch_size_per_process,
            min_combined_score=self.config.min_combined_score,
            max_memory_usage_gb=4.0,  # Reduced memory per process
            max_cpu_usage_percent=50.0  # Reduced CPU per process
        )
    
    def _run_discovery_batch(self, config: ChemistryDiscoveryConfig) -> Dict[str, Any]:
        """Run a single discovery batch (executed in separate process)"""
        try:
            # Initialize discovery engine
            engine = ContinuousChemistryDiscoveryEngine(config)
            
            # Use the private method directly with correct parameter
            batch_result = engine._run_discovery_batch("drug_discovery")
            
            # Extract discoveries from result
            discoveries = batch_result.discoveries if batch_result else []
            
            # Generate batch ID
            batch_id = f"scaled_batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
            
            return {
                'batch_id': batch_id,
                'process_id': os.getpid(),
                'discoveries': discoveries,
                'success': True
            }
            
        except Exception as e:
            logger.error(f"❌ Batch discovery failed in process {os.getpid()}: {e}")
            return {
                'batch_id': f"failed_{os.getpid()}",
                'process_id': os.getpid(),
                'discoveries': [],
                'success': False,
                'error': str(e)
            }
    
    def _update_progress(self, batch_result: Dict[str, Any]):
        """Update progress tracking"""
        if batch_result['success']:
            new_discoveries = batch_result['discoveries']
            for discovery in new_discoveries:
                smiles = discovery.smiles
                if smiles not in self.discovered_molecules:
                    self.discovered_molecules.add(smiles)
                    self.total_discoveries += 1
    
    def _log_progress(self):
        """Log current progress"""
        elapsed = time.time() - self.start_time
        rate = self.total_discoveries / (elapsed / 3600) if elapsed > 0 else 0
        
        logger.info(f"📊 Progress: {self.total_discoveries}/{self.config.total_target_molecules} "
                   f"unique molecules ({rate:.1f}/hour)")
    
    def _log_final_results(self, results: Dict[str, Any]):
        """Log final campaign results"""
        logger.info("🎉 Scaled Discovery Campaign Complete!")
        logger.info(f"📊 Total discoveries: {results['total_discoveries']}")
        logger.info(f"🧬 Unique molecules: {results['unique_molecules']}")
        logger.info(f"⏱️ Duration: {results['duration_seconds']/3600:.2f} hours")
        logger.info(f"🚀 Rate: {results['discoveries_per_hour']:.1f} discoveries/hour")
        
        # Save detailed results
        results_file = f"scaled_discovery_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        import json
        with open(results_file, 'w') as f:
            # Convert ChemicalDiscovery objects to dicts for JSON serialization
            serializable_discoveries = []
            for d in results['discoveries']:
                if hasattr(d, '__dict__'):
                    serializable_discoveries.append(d.__dict__)
                else:
                    serializable_discoveries.append(d)
            
            results_copy = results.copy()
            results_copy['discoveries'] = serializable_discoveries
            json.dump(results_copy, f, indent=2, default=str)
        
        logger.info(f"💾 Results saved to: {results_file}")

def main():
    """Main entry point for scaled discovery"""
    parser = argparse.ArgumentParser(description="Scale FoTChemistry molecular discovery")
    parser.add_argument("--target", type=int, default=1000, 
                       help="Target number of molecules to discover (default: 1000)")
    parser.add_argument("--processes", type=int, default=None,
                       help="Number of parallel processes (default: auto)")
    parser.add_argument("--batch-size", type=int, default=20,
                       help="Molecules per batch per process (default: 20)")
    parser.add_argument("--test", action="store_true",
                       help="Run quick test with 100 molecules")
    
    args = parser.parse_args()
    
    # Configure scaling
    config = ScaledDiscoveryConfig()
    
    if args.test:
        config.total_target_molecules = 100
        config.parallel_processes = 2
        config.max_concurrent_batches = 2
    else:
        config.total_target_molecules = args.target
        if args.processes:
            config.parallel_processes = args.processes
        config.batch_size_per_process = args.batch_size
    
    # Run scaled discovery
    orchestrator = ParallelDiscoveryOrchestrator(config)
    results = orchestrator.run_scaled_discovery()
    
    print(f"\n🎉 Discovery campaign completed!")
    print(f"📊 Found {results['unique_molecules']} unique molecules")
    print(f"⏱️ Completed in {results['duration_seconds']/3600:.2f} hours")
    print(f"🚀 Rate: {results['discoveries_per_hour']:.1f} discoveries/hour")

if __name__ == "__main__":
    main()
