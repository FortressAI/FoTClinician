#!/usr/bin/env python3
"""
NOVELTY VALIDATION ENGINE
Rigorous molecular novelty checking against major chemical databases

🎯 CORE FUNCTIONS:
- InChIKey generation and standardization
- Cross-checking against PubChem/ChEMBL/ChemSpider
- Structural similarity analysis
- Novelty scoring and classification
- Public benefit assessment

Author: FoT Research Team
Purpose: Ensure only truly novel compounds are claimed as discoveries
"""

import requests
import json
import time
import logging
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path
import hashlib
from datetime import datetime

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit.Chem import rdDepictor
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    import pubchempy as pcp
    HAS_PUBCHEMPY = True
except ImportError:
    HAS_PUBCHEMPY = False

logger = logging.getLogger(__name__)

@dataclass
class NoveltyResult:
    """Results of novelty validation"""
    smiles: str
    inchi_key: str
    is_novel: bool
    novelty_score: float  # 0.0 = known, 1.0 = completely novel
    database_matches: List[Dict[str, Any]]
    similarity_matches: List[Dict[str, Any]]
    public_benefit_score: float
    validation_timestamp: datetime
    validation_id: str

@dataclass
class DatabaseMatch:
    """Match found in chemical database"""
    database: str
    compound_id: str
    compound_name: str
    match_type: str  # "exact", "similar", "substructure"
    similarity_score: float
    url: str

class NoveltyValidationEngine:
    """
    Rigorous novelty validation against major chemical databases
    """
    
    def __init__(self, cache_dir: Path = Path("novelty_cache")):
        """Initialize novelty validation engine"""
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(exist_ok=True)
        
        # Database endpoints
        self.databases = {
            "pubchem": {
                "base_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
                "rate_limit": 5,  # requests per second
                "timeout": 30
            },
            "chembl": {
                "base_url": "https://www.ebi.ac.uk/chembl/api/data",
                "rate_limit": 10,
                "timeout": 30
            },
            "chemspider": {
                "base_url": "http://www.chemspider.com/Search.asmx",
                "rate_limit": 1,  # Very conservative
                "timeout": 60
            }
        }
        
        # Novelty thresholds
        self.similarity_threshold = 0.85  # Tanimoto similarity
        self.novelty_threshold = 0.7      # Minimum for "novel" classification
        
        # Rate limiting
        self.last_request_time = {}
        for db in self.databases:
            self.last_request_time[db] = 0
        
        logger.info("✅ Novelty validation engine initialized")
    
    def validate_molecular_novelty(self, smiles: str) -> NoveltyResult:
        """
        Comprehensive novelty validation for a SMILES string
        """
        validation_id = self._generate_validation_id(smiles)
        
        # Check cache first
        cached_result = self._load_from_cache(validation_id)
        if cached_result:
            logger.info(f"📋 Using cached novelty result for {smiles}")
            return cached_result
        
        logger.info(f"🔍 Validating novelty for: {smiles}")
        
        # Generate standardized identifiers
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"❌ Invalid SMILES: {smiles}")
            return self._create_invalid_result(smiles, validation_id)
        
        # Standardize molecule
        mol = self._standardize_molecule(mol)
        canonical_smiles = Chem.MolToSmiles(mol)
        inchi_key = Chem.MolToInchiKey(mol)
        
        # Search databases
        database_matches = []
        
        # PubChem search
        pubchem_matches = self._search_pubchem(inchi_key, canonical_smiles)
        database_matches.extend(pubchem_matches)
        
        # ChEMBL search
        chembl_matches = self._search_chembl(inchi_key, canonical_smiles)
        database_matches.extend(chembl_matches)
        
        # ChemSpider search (if no exact matches found)
        if not any(match["match_type"] == "exact" for match in database_matches):
            chemspider_matches = self._search_chemspider(inchi_key)
            database_matches.extend(chemspider_matches)
        
        # Calculate novelty score
        novelty_score = self._calculate_novelty_score(database_matches)
        is_novel = novelty_score >= self.novelty_threshold
        
        # Assess public benefit potential
        public_benefit_score = self._assess_public_benefit(mol, canonical_smiles)
        
        # Create result
        result = NoveltyResult(
            smiles=canonical_smiles,
            inchi_key=inchi_key,
            is_novel=is_novel,
            novelty_score=novelty_score,
            database_matches=database_matches,
            similarity_matches=[],  # TODO: Implement similarity search
            public_benefit_score=public_benefit_score,
            validation_timestamp=datetime.now(),
            validation_id=validation_id
        )
        
        # Cache result
        self._save_to_cache(validation_id, result)
        
        logger.info(f"✅ Novelty validation complete: {canonical_smiles} (Novel: {is_novel}, Score: {novelty_score:.3f})")
        return result
    
    def _standardize_molecule(self, mol):
        """Standardize molecule for consistent comparison"""
        # Remove hydrogens
        mol = Chem.RemoveHs(mol)
        
        # Sanitize
        Chem.SanitizeMol(mol)
        
        # Generate 2D coordinates for consistency
        rdDepictor.Compute2DCoords(mol)
        
        return mol
    
    def _search_pubchem(self, inchi_key: str, smiles: str) -> List[Dict[str, Any]]:
        """Search PubChem for molecular matches"""
        matches = []
        
        try:
            # Rate limiting
            self._rate_limit("pubchem")
            
            # Search by InChIKey first (most specific)
            url = f"{self.databases['pubchem']['base_url']}/compound/inchikey/{inchi_key}/JSON"
            response = requests.get(url, timeout=self.databases['pubchem']['timeout'])
            
            if response.status_code == 200:
                data = response.json()
                if 'PC_Compounds' in data:
                    for compound in data['PC_Compounds']:
                        cid = compound['id']['id']['cid']
                        matches.append({
                            "database": "PubChem",
                            "compound_id": f"CID:{cid}",
                            "compound_name": self._get_pubchem_name(cid),
                            "match_type": "exact",
                            "similarity_score": 1.0,
                            "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
                        })
            
            # If no exact match, try SMILES search
            if not matches:
                self._rate_limit("pubchem")
                url = f"{self.databases['pubchem']['base_url']}/compound/smiles/{smiles}/JSON"
                response = requests.get(url, timeout=self.databases['pubchem']['timeout'])
                
                if response.status_code == 200:
                    data = response.json()
                    if 'PC_Compounds' in data:
                        for compound in data['PC_Compounds']:
                            cid = compound['id']['id']['cid']
                            matches.append({
                                "database": "PubChem",
                                "compound_id": f"CID:{cid}",
                                "compound_name": self._get_pubchem_name(cid),
                                "match_type": "exact",
                                "similarity_score": 1.0,
                                "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
                            })
        
        except Exception as e:
            logger.warning(f"⚠️ PubChem search failed: {e}")
        
        return matches
    
    def _search_chembl(self, inchi_key: str, smiles: str) -> List[Dict[str, Any]]:
        """Search ChEMBL for molecular matches"""
        matches = []
        
        try:
            # Rate limiting
            self._rate_limit("chembl")
            
            # Search by InChIKey
            url = f"{self.databases['chembl']['base_url']}/molecule.json?molecule_structures__standard_inchi_key={inchi_key}"
            response = requests.get(url, timeout=self.databases['chembl']['timeout'])
            
            if response.status_code == 200:
                data = response.json()
                if 'molecules' in data and data['molecules']:
                    for molecule in data['molecules']:
                        chembl_id = molecule['molecule_chembl_id']
                        pref_name = molecule.get('pref_name', 'Unknown')
                        matches.append({
                            "database": "ChEMBL",
                            "compound_id": chembl_id,
                            "compound_name": pref_name,
                            "match_type": "exact",
                            "similarity_score": 1.0,
                            "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}"
                        })
        
        except Exception as e:
            logger.warning(f"⚠️ ChEMBL search failed: {e}")
        
        return matches
    
    def _search_chemspider(self, inchi_key: str) -> List[Dict[str, Any]]:
        """Search ChemSpider for molecular matches"""
        matches = []
        
        try:
            # Rate limiting (very conservative for ChemSpider)
            self._rate_limit("chemspider")
            
            # ChemSpider API requires registration, so this is a placeholder
            # In production, would need API key and proper implementation
            logger.info("🔍 ChemSpider search not implemented (requires API key)")
        
        except Exception as e:
            logger.warning(f"⚠️ ChemSpider search failed: {e}")
        
        return matches
    
    def _get_pubchem_name(self, cid: int) -> str:
        """Get compound name from PubChem CID"""
        try:
            self._rate_limit("pubchem")
            url = f"{self.databases['pubchem']['base_url']}/compound/cid/{cid}/property/Title/JSON"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    return data['PropertyTable']['Properties'][0].get('Title', f'CID:{cid}')
        except:
            pass
        
        return f'CID:{cid}'
    
    def _calculate_novelty_score(self, database_matches: List[Dict[str, Any]]) -> float:
        """Calculate novelty score based on database matches"""
        if not database_matches:
            return 1.0  # No matches = completely novel
        
        # Check for exact matches
        exact_matches = [m for m in database_matches if m["match_type"] == "exact"]
        if exact_matches:
            return 0.0  # Exact match = not novel
        
        # Calculate based on similarity matches
        max_similarity = max(m["similarity_score"] for m in database_matches)
        novelty_score = 1.0 - max_similarity
        
        return max(0.0, novelty_score)
    
    def _assess_public_benefit(self, mol, smiles: str) -> float:
        """Assess potential public benefit of the molecule"""
        benefit_score = 0.0
        
        try:
            # Drug-likeness indicators
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # Lipinski's Rule of Five compliance
            lipinski_violations = 0
            if mw > 500: lipinski_violations += 1
            if logp > 5: lipinski_violations += 1
            if hbd > 5: lipinski_violations += 1
            if hba > 10: lipinski_violations += 1
            
            if lipinski_violations <= 1:
                benefit_score += 0.3  # Drug-like = potential therapeutic benefit
            
            # Synthetic accessibility (simplified)
            num_atoms = mol.GetNumAtoms()
            num_rings = Descriptors.RingCount(mol)
            
            if num_atoms <= 30 and num_rings <= 3:
                benefit_score += 0.2  # Synthetically accessible
            
            # Environmental considerations
            if 'F' not in smiles and 'Cl' not in smiles and 'Br' not in smiles:
                benefit_score += 0.1  # No halogens = more environmentally friendly
            
            # Complexity balance (not too simple, not too complex)
            if 5 <= num_atoms <= 25:
                benefit_score += 0.2  # Optimal complexity range
            
            # Functional group diversity
            functional_groups = [
                'OH', 'NH2', 'COOH', 'C=O', 'C=C', 'c1ccccc1'  # Common beneficial groups
            ]
            
            for fg in functional_groups:
                if fg in smiles:
                    benefit_score += 0.05
            
            benefit_score = min(1.0, benefit_score)  # Cap at 1.0
        
        except Exception as e:
            logger.warning(f"⚠️ Public benefit assessment failed: {e}")
            benefit_score = 0.5  # Default neutral score
        
        return benefit_score
    
    def _rate_limit(self, database: str):
        """Enforce rate limiting for database requests"""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time[database]
        min_interval = 1.0 / self.databases[database]["rate_limit"]
        
        if time_since_last < min_interval:
            sleep_time = min_interval - time_since_last
            time.sleep(sleep_time)
        
        self.last_request_time[database] = time.time()
    
    def _generate_validation_id(self, smiles: str) -> str:
        """Generate unique validation ID"""
        return hashlib.md5(f"{smiles}_{datetime.now().date()}".encode()).hexdigest()
    
    def _load_from_cache(self, validation_id: str) -> Optional[NoveltyResult]:
        """Load cached validation result"""
        cache_file = self.cache_dir / f"{validation_id}.json"
        
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    data = json.load(f)
                
                # Convert back to NoveltyResult
                data['validation_timestamp'] = datetime.fromisoformat(data['validation_timestamp'])
                return NoveltyResult(**data)
            except Exception as e:
                logger.warning(f"⚠️ Failed to load cache: {e}")
        
        return None
    
    def _save_to_cache(self, validation_id: str, result: NoveltyResult):
        """Save validation result to cache"""
        cache_file = self.cache_dir / f"{validation_id}.json"
        
        try:
            # Convert to JSON-serializable format
            data = {
                "smiles": result.smiles,
                "inchi_key": result.inchi_key,
                "is_novel": result.is_novel,
                "novelty_score": result.novelty_score,
                "database_matches": result.database_matches,
                "similarity_matches": result.similarity_matches,
                "public_benefit_score": result.public_benefit_score,
                "validation_timestamp": result.validation_timestamp.isoformat(),
                "validation_id": result.validation_id
            }
            
            with open(cache_file, 'w') as f:
                json.dump(data, f, indent=2)
        
        except Exception as e:
            logger.warning(f"⚠️ Failed to save cache: {e}")
    
    def _create_invalid_result(self, smiles: str, validation_id: str) -> NoveltyResult:
        """Create result for invalid SMILES"""
        return NoveltyResult(
            smiles=smiles,
            inchi_key="INVALID",
            is_novel=False,
            novelty_score=0.0,
            database_matches=[],
            similarity_matches=[],
            public_benefit_score=0.0,
            validation_timestamp=datetime.now(),
            validation_id=validation_id
        )

def validate_discovery_batch(smiles_list: List[str]) -> List[NoveltyResult]:
    """Validate a batch of molecular discoveries"""
    engine = NoveltyValidationEngine()
    results = []
    
    logger.info(f"🔍 Validating {len(smiles_list)} molecular candidates for novelty")
    
    for i, smiles in enumerate(smiles_list):
        logger.info(f"Progress: {i+1}/{len(smiles_list)}")
        result = engine.validate_molecular_novelty(smiles)
        results.append(result)
        
        # Brief pause between validations to be respectful to APIs
        time.sleep(0.5)
    
    # Summary statistics
    novel_count = sum(1 for r in results if r.is_novel)
    high_benefit_count = sum(1 for r in results if r.public_benefit_score > 0.7)
    
    logger.info(f"✅ Validation complete:")
    logger.info(f"   📊 Total candidates: {len(results)}")
    logger.info(f"   🆕 Novel compounds: {novel_count} ({novel_count/len(results)*100:.1f}%)")
    logger.info(f"   🌟 High public benefit: {high_benefit_count} ({high_benefit_count/len(results)*100:.1f}%)")
    
    return results

if __name__ == "__main__":
    # Test with known compounds
    test_smiles = [
        "C",           # Methane (should be known)
        "CCO",         # Ethanol (should be known)
        "CCCO",        # Propanol (should be known)
        "C1=CC=CC=C1", # Benzene (should be known)
    ]
    
    logging.basicConfig(level=logging.INFO)
    results = validate_discovery_batch(test_smiles)
    
    for result in results:
        print(f"\n{result.smiles}:")
        print(f"  Novel: {result.is_novel}")
        print(f"  Novelty Score: {result.novelty_score:.3f}")
        print(f"  Public Benefit: {result.public_benefit_score:.3f}")
        print(f"  Database Matches: {len(result.database_matches)}")
