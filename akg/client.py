#!/usr/bin/env python3
"""
FoTChemistry Agentic Knowledge Graph (AKG) Client

Provides unified interface to the multi-modal knowledge graph:
- Neo4j: Property graph for molecules, reactions, measurements
- GraphDB: RDF triplestore for ontological reasoning
- Fuseki: SPARQL endpoint for complex queries

All discoveries are recorded with full provenance and FoT validation.
"""

import json
import logging
from typing import Dict, List, Any, Optional
from datetime import datetime
from dataclasses import asdict
import uuid

from neo4j import GraphDatabase
import requests

logger = logging.getLogger(__name__)


class AKG:
    """Agentic Knowledge Graph client for FoTChemistry."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize AKG connections."""
        self.config = config or {
            'neo4j': {
                'uri': 'bolt://localhost:7687',
                'user': 'neo4j', 
                'password': 'fotquantum'
            },
            'graphdb': {
                'uri': 'http://localhost:7200',
                'repository': 'fotchemistry'
            },
            'fuseki': {
                'uri': 'http://localhost:3030',
                'dataset': 'fotchemistry'
            }
        }
        
        # Initialize Neo4j connection (connect to EXISTING instance)
        try:
            self.neo4j_driver = GraphDatabase.driver(
                self.config['neo4j']['uri'],
                auth=(self.config['neo4j']['user'], self.config['neo4j']['password'])
            )
            logger.info("✅ Connected to EXISTING Neo4j instance")
            # Ensure chemistry schema with safe namespacing
            self._ensure_safe_chemistry_schema()
        except Exception as e:
            logger.warning(f"⚠️ Neo4j connection failed: {e}")
            self.neo4j_driver = None
        
        # Initialize HTTP session for SPARQL endpoints
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'FoTChemistry-AKG/1.0',
            'Accept': 'application/json'
        })
    
    def _ensure_safe_chemistry_schema(self):
        """Create chemistry schema with FoTChem namespace (SAFE for existing apps)"""
        if not self.neo4j_driver:
            return
            
        with self.neo4j_driver.session() as session:
            # Create constraints for chemistry entities (isolated namespace)
            # Uses FoTChem_ prefix to avoid conflicts with protein/fluid schemas
            constraints = [
                "CREATE CONSTRAINT FoTChem_molecule_id IF NOT EXISTS FOR (m:FoTChem_Molecule) REQUIRE m.id IS UNIQUE",
                "CREATE CONSTRAINT FoTChem_claim_id IF NOT EXISTS FOR (c:FoTChem_Claim) REQUIRE c.id IS UNIQUE", 
                "CREATE CONSTRAINT FoTChem_verdict_id IF NOT EXISTS FOR (v:FoTChem_Verdict) REQUIRE v.id IS UNIQUE",
                "CREATE CONSTRAINT FoTChem_discovery_id IF NOT EXISTS FOR (d:FoTChem_Discovery) REQUIRE d.id IS UNIQUE",
                "CREATE CONSTRAINT FoTChem_campaign_id IF NOT EXISTS FOR (camp:FoTChem_Campaign) REQUIRE camp.id IS UNIQUE"
            ]
            
            for constraint in constraints:
                try:
                    session.run(constraint)
                    logger.debug(f"✅ Created constraint: {constraint}")
                except Exception as e:
                    # Constraint might already exist, that's fine
                    logger.debug(f"Constraint exists or failed: {e}")
    
    def health_check(self) -> Dict[str, bool]:
        """Check health of all AKG components."""
        health = {}
        
        # Neo4j health
        health['neo4j'] = self.neo4j_health()
        
        # GraphDB health  
        health['graphdb'] = self.graphdb_health()
        
        # Fuseki health
        health['fuseki'] = self.fuseki_health()
        
        return health
    
    def neo4j_health(self) -> bool:
        """Check Neo4j connectivity."""
        if not self.neo4j_driver:
            return False
        
        try:
            with self.neo4j_driver.session() as session:
                result = session.run("RETURN 1 as health")
                return result.single()['health'] == 1
        except Exception as e:
            logger.warning(f"Neo4j health check failed: {e}")
            return False
    
    def graphdb_health(self) -> bool:
        """Check GraphDB connectivity."""
        try:
            response = self.session.get(
                f"{self.config['graphdb']['uri']}/rest/repositories",
                timeout=5
            )
            return response.status_code == 200
        except Exception as e:
            logger.warning(f"GraphDB health check failed: {e}")
            return False
    
    def fuseki_health(self) -> bool:
        """Check Fuseki connectivity."""
        try:
            response = self.session.get(
                f"{self.config['fuseki']['uri']}/$/ping",
                timeout=5
            )
            return response.status_code == 200
        except Exception as e:
            logger.warning(f"Fuseki health check failed: {e}")
            return False
    
    def record_discovery_verdict(self, verdict: Dict[str, Any]) -> str:
        """Record a discovery verdict with full provenance."""
        verdict_id = str(uuid.uuid4())
        timestamp = datetime.now().isoformat()
        
        try:
            # Store in Neo4j property graph
            self._store_verdict_neo4j(verdict_id, verdict, timestamp)
            
            # Store in GraphDB as RDF
            self._store_verdict_rdf(verdict_id, verdict, timestamp)
            
            logger.info(f"📝 Recorded discovery verdict: {verdict_id}")
            return verdict_id
            
        except Exception as e:
            logger.error(f"❌ Failed to record verdict: {e}")
            raise
    
    def _store_verdict_neo4j(self, verdict_id: str, verdict: Dict[str, Any], timestamp: str):
        """Store verdict in Neo4j property graph (SAFE namespace)."""
        if not self.neo4j_driver:
            logger.warning("Neo4j not available, skipping storage")
            return
        
        with self.neo4j_driver.session() as session:
            # Create verdict node with FoTChem namespace (SAFE for existing apps)
            session.run("""
                CREATE (v:FoTChem_Verdict {
                    id: $verdict_id,
                    status: $status,
                    confidence: $confidence,
                    evidence_strength: $evidence_strength,
                    virtue_score: $virtue_score,
                    replication_count: $replication_count,
                    reasoning: $reasoning,
                    recommendation: $recommendation,
                    timestamp: $timestamp,
                    raw_data: $raw_data
                })
            """, 
                verdict_id=verdict_id,
                status=verdict.get('status', 'unknown'),
                confidence=verdict.get('confidence', 0.0),
                evidence_strength=verdict.get('evidence_strength', 0.0),
                virtue_score=verdict.get('virtue_score', 0.0),
                replication_count=verdict.get('replication_count', 0),
                reasoning=verdict.get('reasoning', ''),
                recommendation=verdict.get('recommendation', ''),
                timestamp=timestamp,
                raw_data=json.dumps(verdict)
            )
            
            logger.debug(f"✅ Stored verdict in Neo4j (FoTChem namespace): {verdict_id}")
    
    def _store_verdict_rdf(self, verdict_id: str, verdict: Dict[str, Any], timestamp: str):
        """Store verdict in GraphDB as RDF triples."""
        try:
            # Construct RDF triples using FoT Chemistry ontology
            triples = f"""
            @prefix fot: <http://fieldoftruth.org/ontology/chemistry#> .
            @prefix prov: <http://www.w3.org/ns/prov#> .
            @prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
            
            fot:verdict_{verdict_id} a fot:Verdict ;
                fot:hasStatus "{verdict.get('status', 'unknown')}" ;
                fot:hasConfidence "{verdict.get('confidence', 0.0)}"^^xsd:double ;
                fot:hasEvidenceStrength "{verdict.get('evidence_strength', 0.0)}"^^xsd:double ;
                fot:hasVirtueScore "{verdict.get('virtue_score', 0.0)}"^^xsd:double ;
                fot:hasReplicationCount "{verdict.get('replication_count', 0)}"^^xsd:integer ;
                prov:generatedAtTime "{timestamp}"^^xsd:dateTime ;
                fot:hasReasoning "{verdict.get('reasoning', '')}" ;
                fot:hasRecommendation "{verdict.get('recommendation', '')}" .
            """
            
            # Store in GraphDB
            response = self.session.post(
                f"{self.config['graphdb']['uri']}/repositories/{self.config['graphdb']['repository']}/statements",
                data=triples,
                headers={'Content-Type': 'text/turtle'},
                timeout=10
            )
            
            if response.status_code in [200, 204]:
                logger.debug(f"✅ Stored verdict in GraphDB: {verdict_id}")
            else:
                logger.warning(f"⚠️ GraphDB storage failed: {response.status_code}")
                
        except Exception as e:
            logger.warning(f"⚠️ GraphDB storage failed: {e}")
    
    def query_new_molecules(self, since: datetime) -> List[Dict[str, Any]]:
        """Query for new chemistry molecules since timestamp (SAFE namespace)."""
        if not self.neo4j_driver:
            return []
        
        try:
            with self.neo4j_driver.session() as session:
                result = session.run("""
                    MATCH (m:FoTChem_Molecule)
                    WHERE m.created >= $since
                    RETURN m.id as id, m.smiles as smiles, m.inchi as inchi, 
                           m.created as created
                    ORDER BY m.created DESC
                    LIMIT 100
                """, since=since.isoformat())
                
                molecules = []
                for record in result:
                    molecules.append({
                        'id': record['id'],
                        'smiles': record['smiles'],
                        'inchi': record['inchi'],
                        'created': record['created']
                    })
                
                return molecules
                
        except Exception as e:
            logger.warning(f"⚠️ Query new molecules failed: {e}")
            return []
    
    def query_new_reactions(self, since: datetime) -> List[Dict[str, Any]]:
        """Query for new reactions since timestamp."""
        if not self.neo4j_driver:
            return []
        
        try:
            with self.neo4j_driver.session() as session:
                result = session.run("""
                    MATCH (r:Reaction)
                    WHERE r.created >= $since
                    RETURN r.id as id, r.reaction_smiles as reaction_smiles,
                           r.yield as yield, r.created as created
                    ORDER BY r.created DESC
                    LIMIT 100
                """, since=since.isoformat())
                
                reactions = []
                for record in result:
                    reactions.append({
                        'id': record['id'],
                        'reaction_smiles': record['reaction_smiles'],
                        'yield': record['yield'],
                        'created': record['created']
                    })
                
                return reactions
                
        except Exception as e:
            logger.warning(f"⚠️ Query new reactions failed: {e}")
            return []
    
    def query_new_measurements(self, since: datetime) -> List[Dict[str, Any]]:
        """Query for new measurements since timestamp."""
        if not self.neo4j_driver:
            return []
        
        try:
            with self.neo4j_driver.session() as session:
                result = session.run("""
                    MATCH (m:Measurement)
                    WHERE m.created >= $since
                    RETURN m.id as id, m.property as property, m.value as value,
                           m.uncertainty as uncertainty, m.created as created
                    ORDER BY m.created DESC
                    LIMIT 100
                """, since=since.isoformat())
                
                measurements = []
                for record in result:
                    measurements.append({
                        'id': record['id'],
                        'property': record['property'],
                        'value': record['value'],
                        'uncertainty': record['uncertainty'],
                        'created': record['created']
                    })
                
                return measurements
                
        except Exception as e:
            logger.warning(f"⚠️ Query new measurements failed: {e}")
            return []
    
    def query_discovery_statistics(self) -> Dict[str, Any]:
        """Query overall discovery statistics."""
        stats = {
            'total_discoveries': 0,
            'truth_collapsed': 0,
            'refuted_claims': 0,
            'needs_evidence': 0,
            'by_campaign': {},
            'recent_discoveries': []
        }
        
        if not self.neo4j_driver:
            return stats
        
        try:
            with self.neo4j_driver.session() as session:
                # Total verdicts by status (SAFE namespace)
                result = session.run("""
                    MATCH (v:FoTChem_Verdict)
                    RETURN v.status as status, count(v) as count
                """)
                
                for record in result:
                    status = record['status']
                    count = record['count']
                    
                    if status == 'truth':
                        stats['truth_collapsed'] = count
                    elif status == 'refuted':
                        stats['refuted_claims'] = count
                    elif status == 'needs_evidence':
                        stats['needs_evidence'] = count
                    
                    stats['total_discoveries'] += count
                
                # Recent discoveries (SAFE namespace)
                result = session.run("""
                    MATCH (v:FoTChem_Verdict)
                    WHERE v.status = 'truth'
                    RETURN v.id as id, v.reasoning as reasoning, v.timestamp as timestamp
                    ORDER BY v.timestamp DESC
                    LIMIT 10
                """)
                
                for record in result:
                    stats['recent_discoveries'].append({
                        'id': record['id'],
                        'reasoning': record['reasoning'],
                        'timestamp': record['timestamp']
                    })
                
                return stats
                
        except Exception as e:
            logger.warning(f"⚠️ Query statistics failed: {e}")
            return stats
    
    def sparql_query(self, query: str) -> List[Dict[str, Any]]:
        """Execute SPARQL query against the triplestore."""
        try:
            response = self.session.post(
                f"{self.config['fuseki']['uri']}/{self.config['fuseki']['dataset']}/query",
                data={'query': query},
                headers={'Accept': 'application/json'},
                timeout=30
            )
            
            if response.status_code == 200:
                results = response.json()
                return results.get('results', {}).get('bindings', [])
            else:
                logger.warning(f"SPARQL query failed: {response.status_code}")
                return []
                
        except Exception as e:
            logger.warning(f"SPARQL query failed: {e}")
            return []
    
    def store_claim(self, claim: Dict[str, Any]) -> str:
        """Store a new chemistry claim (SAFE namespace)"""
        claim_id = claim.get('id', str(uuid.uuid4()))
        
        if not self.neo4j_driver:
            logger.warning("Neo4j not available, skipping claim storage")
            return claim_id
            
        with self.neo4j_driver.session() as session:
            result = session.run("""
                CREATE (c:FoTChem_Claim {
                    id: $claim_id,
                    objective: $objective,
                    status: 'active',
                    virtue_weighting: $virtue_weighting,
                    collapse_rules: $collapse_rules,
                    created_at: datetime(),
                    campaign: $campaign
                })
                RETURN c.id as claim_id
            """, 
                claim_id=claim_id,
                objective=claim.get('objective', ''),
                virtue_weighting=json.dumps(claim.get('virtue_weighting', {})),
                collapse_rules=json.dumps(claim.get('collapse_rules', {})),
                campaign=claim.get('campaign', 'unknown')
            )
            
            logger.info(f"📝 Stored chemistry claim: {claim_id}")
            return result.single()['claim_id']
            
    def get_claims(self, status: Optional[str] = None, campaign: Optional[str] = None) -> List[Dict]:
        """Retrieve chemistry claims (SAFE namespace)"""
        if not self.neo4j_driver:
            return []
            
        with self.neo4j_driver.session() as session:
            where_clause = "WHERE 1=1"
            params = {}
            
            if status:
                where_clause += " AND c.status = $status"
                params['status'] = status
                
            if campaign:
                where_clause += " AND c.campaign = $campaign" 
                params['campaign'] = campaign
                
            result = session.run(f"""
                MATCH (c:FoTChem_Claim)
                {where_clause}
                RETURN c.id as id, c.objective as objective, c.status as status,
                       c.virtue_weighting as virtue_weighting, c.collapse_rules as collapse_rules,
                       c.created_at as created_at, c.campaign as campaign
                ORDER BY c.created_at DESC
            """, **params)
            
            claims = []
            for record in result:
                claim = dict(record)
                # Parse JSON fields
                claim['virtue_weighting'] = json.loads(claim['virtue_weighting'] or '{}')
                claim['collapse_rules'] = json.loads(claim['collapse_rules'] or '{}')
                claims.append(claim)
                
            return claims
            
    def store_evidence(self, claim_id: str, evidence: Dict[str, Any]) -> str:
        """Store evidence for a chemistry claim (SAFE namespace)"""
        evidence_id = str(uuid.uuid4())
        
        if not self.neo4j_driver:
            logger.warning("Neo4j not available, skipping evidence storage")
            return evidence_id
        
        with self.neo4j_driver.session() as session:
            session.run("""
                MATCH (c:FoTChem_Claim {id: $claim_id})
                CREATE (e:FoTChem_Evidence {
                    id: $evidence_id,
                    metrics: $metrics,
                    uncertainty: $uncertainty,
                    virtue_vector: $virtue_vector,
                    generated_at: datetime(),
                    agent_type: $agent_type
                })
                CREATE (c)-[:HAS_EVIDENCE]->(e)
            """,
                claim_id=claim_id,
                evidence_id=evidence_id,
                metrics=json.dumps(evidence.get('metrics', {})),
                uncertainty=evidence.get('uncertainty', 1.0),
                virtue_vector=json.dumps(evidence.get('virtue_vector', [])),
                agent_type=evidence.get('agent_type', 'unknown')
            )
            
        logger.info(f"📊 Stored evidence: {evidence_id} for claim: {claim_id}")
        return evidence_id
        
    def collapse_claim(self, claim_id: str, verdict: str, virtues: List[float], evidence: Dict) -> bool:
        """Collapse a claim to truth/refute/needs-evidence (SAFE namespace)"""
        if not self.neo4j_driver:
            logger.warning("Neo4j not available, skipping claim collapse")
            return False
            
        with self.neo4j_driver.session() as session:
            session.run("""
                MATCH (c:FoTChem_Claim {id: $claim_id})
                SET c.status = $verdict,
                    c.final_virtues = $virtues,
                    c.collapsed_at = datetime(),
                    c.final_evidence = $evidence
                CREATE (d:FoTChem_Discovery {
                    id: $discovery_id,
                    claim_id: $claim_id,
                    verdict: $verdict,
                    virtues: $virtues,
                    discovered_at: datetime()
                })
                CREATE (c)-[:COLLAPSED_TO]->(d)
            """,
                claim_id=claim_id,
                verdict=verdict,
                virtues=virtues,
                evidence=json.dumps(evidence),
                discovery_id=str(uuid.uuid4())
            )
            
        logger.info(f"🎯 Collapsed claim {claim_id} to {verdict}")
        return True
        
    def export_for_streamlit(self, output_file: str = "results/chemistry_discoveries.json"):
        """Export chemistry data for Streamlit dashboard (Git-trackable)"""
        import os
        os.makedirs("results", exist_ok=True)
        
        if not self.neo4j_driver:
            logger.warning("Neo4j not available, creating empty export")
            empty_data = {
                "export_timestamp": datetime.now().isoformat(),
                "total_discoveries": 0,
                "discoveries": []
            }
            with open(output_file, 'w') as f:
                json.dump(empty_data, f, indent=2)
            return empty_data
        
        with self.neo4j_driver.session() as session:
            # Get all chemistry discoveries with claims and evidence (SAFE namespace)
            result = session.run("""
                MATCH (d:FoTChem_Discovery)
                OPTIONAL MATCH (d)<-[:COLLAPSED_TO]-(c:FoTChem_Claim)
                OPTIONAL MATCH (c)-[:HAS_EVIDENCE]->(e:FoTChem_Evidence)
                RETURN d, c, collect(e) as evidence
                ORDER BY d.discovered_at DESC
            """)
            
            export_data = {
                "export_timestamp": datetime.now().isoformat(),
                "total_discoveries": 0,
                "discoveries": []
            }
            
            for record in result:
                discovery = dict(record['d']) if record['d'] else {}
                claim = dict(record['c']) if record['c'] else {}
                evidence_list = [dict(e) for e in record['evidence']] if record['evidence'] else []
                
                export_data["discoveries"].append({
                    "discovery": discovery,
                    "claim": claim,
                    "evidence": evidence_list
                })
                
            export_data["total_discoveries"] = len(export_data["discoveries"])
            
            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)
                
            logger.info(f"📁 Exported {export_data['total_discoveries']} discoveries to {output_file}")
            return export_data
    
    def close(self):
        """Close all connections."""
        if self.neo4j_driver:
            self.neo4j_driver.close()
        
        self.session.close()


# Example usage and testing
if __name__ == "__main__":
    # Test AKG client
    akg = AKG()
    
    # Check health
    health = akg.health_check()
    print(f"🏥 AKG Health: {health}")
    
    # Test recording a discovery verdict
    test_verdict = {
        'status': 'truth',
        'confidence': 0.87,
        'evidence_strength': 0.92,
        'virtue_score': 0.85,
        'replication_count': 3,
        'reasoning': 'All criteria met with high confidence',
        'recommendation': 'Collapse to truth and publish'
    }
    
    try:
        verdict_id = akg.record_discovery_verdict(test_verdict)
        print(f"✅ Recorded test verdict: {verdict_id}")
    except Exception as e:
        print(f"❌ Failed to record verdict: {e}")
    
    # Query statistics
    stats = akg.query_discovery_statistics()
    print(f"📊 Discovery statistics: {stats}")
    
    # Close connections
    akg.close()
