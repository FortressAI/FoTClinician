#!/usr/bin/env python3
"""
Neo4j Problem-Solution Loader

Loads FoTChem problem-solution claims into the existing Neo4j database
while preserving existing protein and discovery data.

Uses FoTChem_ namespace to avoid conflicts.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Any
from neo4j import GraphDatabase
import uuid
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FoTChemNeo4jLoader:
    """Loads problem-solution data into Neo4j with FoTChem namespace"""
    
    def __init__(self, uri="bolt://localhost:7687", user="neo4j", password="fotquantum"):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
        
    def close(self):
        """Close the database connection"""
        self.driver.close()
    
    def create_schema(self):
        """Create indexes and constraints for FoTChem nodes"""
        with self.driver.session() as session:
            # Create constraints for unique IDs
            constraints = [
                "CREATE CONSTRAINT fotchem_compound_id IF NOT EXISTS FOR (c:FoTChem_Compound) REQUIRE c.id IS UNIQUE",
                "CREATE CONSTRAINT fotchem_problem_id IF NOT EXISTS FOR (p:FoTChem_Problem) REQUIRE p.id IS UNIQUE", 
                "CREATE CONSTRAINT fotchem_claim_id IF NOT EXISTS FOR (cl:FoTChem_Claim) REQUIRE cl.id IS UNIQUE",
                "CREATE CONSTRAINT fotchem_measurement_id IF NOT EXISTS FOR (m:FoTChem_Measurement) REQUIRE m.id IS UNIQUE"
            ]
            
            # Create indexes for performance
            indexes = [
                "CREATE INDEX fotchem_compound_smiles IF NOT EXISTS FOR (c:FoTChem_Compound) ON (c.smiles)",
                "CREATE INDEX fotchem_problem_name IF NOT EXISTS FOR (p:FoTChem_Problem) ON (p.name)",
                "CREATE INDEX fotchem_claim_verdict IF NOT EXISTS FOR (cl:FoTChem_Claim) ON (cl.verdict_ok)",
                "CREATE INDEX fotchem_measurement_value IF NOT EXISTS FOR (m:FoTChem_Measurement) ON (m.value)"
            ]
            
            for constraint in constraints:
                try:
                    session.run(constraint)
                    logger.info(f"✅ Created constraint: {constraint.split('FOR')[1].split('REQUIRE')[0]}")
                except Exception as e:
                    logger.warning(f"Constraint might already exist: {e}")
                    
            for index in indexes:
                try:
                    session.run(index)
                    logger.info(f"✅ Created index: {index.split('FOR')[1].split('ON')[0]}")
                except Exception as e:
                    logger.warning(f"Index might already exist: {e}")
    
    def load_problems(self):
        """Load the predefined chemistry problems"""
        problems = [
            {
                "id": "fct:PFASRemovalProblem",
                "name": "PFAS Removal", 
                "description": "Remove PFAS contaminants to below threshold",
                "threshold_value": 25.0,
                "threshold_unit": "ng/L",
                "higher_is_better": False
            },
            {
                "id": "fct:CO2ElectroReductionProblem",
                "name": "CO2 Electrocatalysis",
                "description": "Efficient CO2 electroreduction to CO",
                "threshold_value": 0.65,
                "threshold_unit": "fraction", 
                "higher_is_better": True
            },
            {
                "id": "fct:PKaLogPConsistencyProblem",
                "name": "Thermodynamic Consistency",
                "description": "Thermodynamic cycle closure within error",
                "threshold_value": 0.3,
                "threshold_unit": "kcal/mol",
                "higher_is_better": False
            },
            {
                "id": "fct:GreenSynthesisProblem", 
                "name": "Green Synthesis",
                "description": "Sustainable synthesis improvements",
                "threshold_value": 0.5,
                "threshold_unit": "fraction",
                "higher_is_better": True
            }
        ]
        
        with self.driver.session() as session:
            for problem in problems:
                result = session.run("""
                    MERGE (p:FoTChem_Problem {id: $id})
                    SET p.name = $name,
                        p.description = $description,
                        p.threshold_value = $threshold_value,
                        p.threshold_unit = $threshold_unit,
                        p.higher_is_better = $higher_is_better,
                        p.created_at = datetime()
                    RETURN p.id as problem_id
                """, **problem)
                
                record = result.single()
                logger.info(f"✅ Loaded problem: {record['problem_id']}")
    
    def load_compounds_from_claims(self, claims_file: str):
        """Load compounds referenced in claims"""
        logger.info(f"📦 Loading compounds from {claims_file}")
        
        compounds = {}
        
        # Extract unique compounds from claims
        with open(claims_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                try:
                    claim = json.loads(line.strip())
                    compound_ref = claim.get('aboutCompound', '')
                    compound_id = compound_ref.replace('compound:', '')
                    
                    if compound_id and compound_id not in compounds:
                        # Try to get SMILES from the original discovery data
                        compounds[compound_id] = {
                            "id": compound_id,
                            "reference": compound_ref,
                            "smiles": None  # Will try to fetch from main dataset
                        }
                        
                except Exception as e:
                    logger.warning(f"Error parsing claim on line {line_num}: {e}")
        
        logger.info(f"📊 Found {len(compounds)} unique compounds in claims")
        
        # Try to get SMILES from main discovery dataset
        try:
            with open("results/chemistry_discoveries.json", 'r') as f:
                discoveries = json.load(f)
                
            for discovery in discoveries.get('discoveries', []):
                discovery_id = discovery.get('discovery_id', '')
                smiles = discovery.get('smiles', '')
                
                if discovery_id in compounds and smiles:
                    compounds[discovery_id]['smiles'] = smiles
                    
            logger.info(f"🧬 Matched SMILES for compounds from discovery dataset")
            
        except Exception as e:
            logger.warning(f"Could not load main discovery dataset: {e}")
        
        # Load compounds into Neo4j
        with self.driver.session() as session:
            for compound_id, compound_data in compounds.items():
                result = session.run("""
                    MERGE (c:FoTChem_Compound {id: $id})
                    SET c.smiles = $smiles,
                        c.reference = $reference,
                        c.created_at = datetime()
                    RETURN c.id as compound_id
                """, 
                id=compound_id,
                smiles=compound_data['smiles'],
                reference=compound_data['reference']
                )
                
                if (len(compounds) <= 10) or (len([k for k in compounds.keys() if compounds[k] == compound_data]) % 1000 == 0):
                    logger.info(f"✅ Loaded compound: {compound_id}")
        
        logger.info(f"📦 Loaded {len(compounds)} compounds into Neo4j")
        return len(compounds)
    
    def load_claims(self, claims_file: str):
        """Load problem-solution claims from JSONL file"""
        logger.info(f"📋 Loading claims from {claims_file}")
        
        claims_loaded = 0
        
        with self.driver.session() as session:
            with open(claims_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    try:
                        claim = json.loads(line.strip())
                        
                        # Extract claim data
                        claim_id = claim.get('id', f'claim_{uuid.uuid4()}')
                        compound_ref = claim.get('aboutCompound', '').replace('compound:', '')
                        problem_ref = claim.get('addressesProblem', '')
                        
                        # Extract measurement data
                        measurement = claim.get('hasMeasurement', {})
                        measurement_id = measurement.get('id', f'measurement_{uuid.uuid4()}')
                        metric_name = measurement.get('hasMetric', '').replace('fct:', '')
                        value = measurement.get('value', 0.0)
                        unit = measurement.get('unit', '')
                        uncertainty = measurement.get('uncertainty', 0.0)
                        
                        # Extract virtue vector
                        virtue_vector = claim.get('hasVirtueVector', {})
                        honesty = virtue_vector.get('honesty', 0.0)
                        prudence = virtue_vector.get('prudence', 0.0)
                        temperance = virtue_vector.get('temperance', 0.0)
                        beneficence = virtue_vector.get('beneficence', 0.0)
                        virtue_sum = honesty + prudence + temperance + beneficence
                        
                        # Extract collapse policy
                        collapse_policy = claim.get('hasCollapsePolicy', {})
                        required_replications = collapse_policy.get('replicationsCount', 1)
                        max_uncertainty = collapse_policy.get('uncertainty', 1.0)
                        
                        # Extract verdict
                        verdict = claim.get('resultsInVerdict', {})
                        verdict_ok = verdict.get('ok', False)
                        
                        # Create claim with relationships
                        result = session.run("""
                            // Find or create measurement
                            MERGE (m:FoTChem_Measurement {id: $measurement_id})
                            SET m.metric_name = $metric_name,
                                m.value = $value,
                                m.unit = $unit,
                                m.uncertainty = $uncertainty,
                                m.created_at = datetime()
                            
                            // Find compound and problem
                            MATCH (c:FoTChem_Compound {id: $compound_id})
                            MATCH (p:FoTChem_Problem {id: $problem_id})
                            
                            // Create claim
                            MERGE (cl:FoTChem_Claim {id: $claim_id})
                            SET cl.honesty = $honesty,
                                cl.prudence = $prudence,
                                cl.temperance = $temperance,
                                cl.beneficence = $beneficence,
                                cl.virtue_sum = $virtue_sum,
                                cl.required_replications = $required_replications,
                                cl.max_uncertainty = $max_uncertainty,
                                cl.verdict_ok = $verdict_ok,
                                cl.created_at = datetime()
                            
                            // Create relationships
                            MERGE (cl)-[:ABOUT_COMPOUND]->(c)
                            MERGE (cl)-[:ADDRESSES_PROBLEM]->(p) 
                            MERGE (cl)-[:HAS_MEASUREMENT]->(m)
                            
                            RETURN cl.id as claim_id
                        """,
                        measurement_id=measurement_id,
                        metric_name=metric_name,
                        value=float(value),
                        unit=unit,
                        uncertainty=float(uncertainty),
                        compound_id=compound_ref,
                        problem_id=problem_ref,
                        claim_id=claim_id,
                        honesty=float(honesty),
                        prudence=float(prudence),
                        temperance=float(temperance),
                        beneficence=float(beneficence),
                        virtue_sum=float(virtue_sum),
                        required_replications=int(required_replications),
                        max_uncertainty=float(max_uncertainty),
                        verdict_ok=bool(verdict_ok)
                        )
                        
                        claims_loaded += 1
                        
                        if claims_loaded % 1000 == 0:
                            logger.info(f"   Loaded {claims_loaded} claims...")
                            
                    except Exception as e:
                        logger.warning(f"Error loading claim on line {line_num}: {e}")
                        continue
        
        logger.info(f"📋 Loaded {claims_loaded} claims into Neo4j")
        return claims_loaded
    
    def create_summary_views(self):
        """Create summary data for easy querying"""
        with self.driver.session() as session:
            
            # Create problem summary statistics
            session.run("""
                MATCH (p:FoTChem_Problem)<-[:ADDRESSES_PROBLEM]-(cl:FoTChem_Claim)
                WITH p, 
                     count(cl) as total_claims,
                     sum(CASE WHEN cl.verdict_ok THEN 1 ELSE 0 END) as valid_claims
                SET p.total_claims = total_claims,
                    p.valid_claims = valid_claims,
                    p.success_rate = toFloat(valid_claims) / total_claims
            """)
            
            # Create compound problem-solving capabilities
            session.run("""
                MATCH (c:FoTChem_Compound)<-[:ABOUT_COMPOUND]-(cl:FoTChem_Claim)
                WHERE cl.verdict_ok = true
                WITH c, collect(DISTINCT cl) as valid_claims
                SET c.problems_solved = size(valid_claims),
                    c.avg_virtue_sum = reduce(sum = 0.0, cl IN valid_claims | sum + cl.virtue_sum) / size(valid_claims)
            """)
            
            logger.info("✅ Created summary statistics")
    
    def run_test_queries(self):
        """Run test queries to verify data loading"""
        with self.driver.session() as session:
            
            # Test 1: Count nodes by type
            result = session.run("""
                RETURN 
                    size((:FoTChem_Problem)) as problems,
                    size((:FoTChem_Compound)) as compounds, 
                    size((:FoTChem_Claim)) as claims,
                    size((:FoTChem_Measurement)) as measurements
            """)
            record = result.single()
            logger.info(f"📊 Node counts: {dict(record)}")
            
            # Test 2: Top PFAS solutions
            result = session.run("""
                MATCH (p:FoTChem_Problem {name: 'PFAS Removal'})<-[:ADDRESSES_PROBLEM]-(cl:FoTChem_Claim)-[:ABOUT_COMPOUND]->(c:FoTChem_Compound)
                MATCH (cl)-[:HAS_MEASUREMENT]->(m:FoTChem_Measurement)
                WHERE cl.verdict_ok = true
                RETURN c.id, c.smiles, m.value as pfas_residual, cl.virtue_sum
                ORDER BY m.value ASC
                LIMIT 5
            """)
            
            logger.info("🏆 Top 5 PFAS removal solutions:")
            for record in result:
                logger.info(f"   {record['c.id']}: {record['pfas_residual']:.1f} ng/L (virtue: {record['cl.virtue_sum']:.2f})")
            
            # Test 3: Problem success rates
            result = session.run("""
                MATCH (p:FoTChem_Problem)
                RETURN p.name, p.total_claims, p.valid_claims, p.success_rate
                ORDER BY p.success_rate DESC
            """)
            
            logger.info("📈 Problem success rates:")
            for record in result:
                logger.info(f"   {record['p.name']}: {record['p.valid_claims']}/{record['p.total_claims']} ({record['p.success_rate']:.1%})")

def main():
    """Main loader function"""
    print("🔗 FoTChem Neo4j Problem-Solution Loader")
    print("=" * 50)
    
    loader = FoTChemNeo4jLoader()
    
    try:
        # Create schema
        logger.info("🏗️ Creating Neo4j schema...")
        loader.create_schema()
        
        # Load problems
        logger.info("🧪 Loading chemistry problems...")
        loader.load_problems()
        
        # Load compounds
        logger.info("🧬 Loading compounds...")
        compounds_loaded = loader.load_compounds_from_claims("problem_solution_analysis/problem_solution_claims.jsonl")
        
        # Load claims
        logger.info("📋 Loading problem-solution claims...")
        claims_loaded = loader.load_claims("problem_solution_analysis/problem_solution_claims.jsonl")
        
        # Create summary views
        logger.info("📊 Creating summary statistics...")
        loader.create_summary_views()
        
        # Run test queries
        logger.info("🧪 Running validation queries...")
        loader.run_test_queries()
        
        print("\n🎉 SUCCESS: FoTChem Data Loaded!")
        print("=" * 50)
        print(f"📦 Compounds loaded: {compounds_loaded}")
        print(f"📋 Claims loaded: {claims_loaded}")
        print(f"🧪 Problems: 4 chemistry challenges")
        print(f"🔗 Relationships: Claims → Compounds → Problems")
        print()
        print("🔍 Query examples:")
        print("   // Find PFAS solutions")
        print("   MATCH (p:FoTChem_Problem {name: 'PFAS Removal'})<-[:ADDRESSES_PROBLEM]-(cl:FoTChem_Claim)")
        print("   WHERE cl.verdict_ok = true")
        print("   RETURN cl, p")
        print()
        print("   // Find multi-problem solvers")
        print("   MATCH (c:FoTChem_Compound)")
        print("   WHERE c.problems_solved > 1")
        print("   RETURN c ORDER BY c.problems_solved DESC")
        
    except Exception as e:
        logger.error(f"❌ Error loading data: {e}")
        raise
    finally:
        loader.close()

if __name__ == "__main__":
    main()
