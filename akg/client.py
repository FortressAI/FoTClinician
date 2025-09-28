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
                'password': 'fotchemistry123'
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
        
        # Initialize Neo4j connection
        try:
            self.neo4j_driver = GraphDatabase.driver(
                self.config['neo4j']['uri'],
                auth=(self.config['neo4j']['user'], self.config['neo4j']['password'])
            )
            logger.info("✅ Connected to Neo4j")
        except Exception as e:
            logger.warning(f"⚠️ Neo4j connection failed: {e}")
            self.neo4j_driver = None
        
        # Initialize HTTP session for SPARQL endpoints
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'FoTChemistry-AKG/1.0',
            'Accept': 'application/json'
        })
    
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
        """Store verdict in Neo4j property graph."""
        if not self.neo4j_driver:
            logger.warning("Neo4j not available, skipping storage")
            return
        
        with self.neo4j_driver.session() as session:
            # Create verdict node
            session.run("""
                CREATE (v:Verdict {
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
            
            logger.debug(f"✅ Stored verdict in Neo4j: {verdict_id}")
    
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
        """Query for new molecules since timestamp."""
        if not self.neo4j_driver:
            return []
        
        try:
            with self.neo4j_driver.session() as session:
                result = session.run("""
                    MATCH (m:Molecule)
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
                # Total verdicts by status
                result = session.run("""
                    MATCH (v:Verdict)
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
                
                # Recent discoveries
                result = session.run("""
                    MATCH (v:Verdict)
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
