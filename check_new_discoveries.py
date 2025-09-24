#!/usr/bin/env python3
"""
Check for New Protein Discoveries Since Last Export

Quickly checks if we have new discoveries in Neo4j that need to be exported
to update GitHub wiki and Streamlit data.
"""

import sys
import os
from datetime import datetime
from neo4j_discovery_engine import Neo4jDiscoveryEngine
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def check_for_new_discoveries():
    """Check for new discoveries since last export timestamp"""
    
    print("🔍 CHECKING FOR NEW PROTEIN DISCOVERIES")
    print("=" * 50)
    
    # Last export timestamp from our data files
    last_export_time = "2025-09-18T08:47:28.829671"
    
    try:
        # Connect to Neo4j
        print("🔗 Connecting to Neo4j discovery database...")
        engine = Neo4jDiscoveryEngine()
        
        with engine.driver.session() as session:
            # Check for new discoveries
            new_discoveries_query = """
            MATCH (d:Discovery)
            WHERE d.timestamp > datetime($last_export)
            RETURN count(d) as new_count
            """
            
            result = session.run(new_discoveries_query, {'last_export': last_export_time})
            new_count = result.single()['new_count']
            
            print(f"📅 Last export: {last_export_time}")
            print(f"🔢 New discoveries since export: {new_count}")
            
            if new_count > 0:
                print(f"🆕 FOUND {new_count} NEW DISCOVERIES!")
                print("📊 Getting discovery details...")
                
                # Get details of new discoveries
                details_query = """
                MATCH (d:Discovery)-[:HAS_SEQUENCE]->(s:Sequence)
                WHERE d.timestamp > datetime($last_export)
                RETURN d.id as discovery_id,
                       s.length as sequence_length,
                       d.validation_score as validation_score,
                       d.energy_kcal_mol as energy,
                       d.quantum_coherence as quantum_coherence,
                       d.timestamp as discovery_date
                ORDER BY d.timestamp DESC
                LIMIT 10
                """
                
                details_result = session.run(details_query, {'last_export': last_export_time})
                
                print("\n📋 Recent Discoveries:")
                print("-" * 80)
                for record in details_result:
                    timestamp = record['discovery_date'].strftime('%Y-%m-%d %H:%M:%S') if record['discovery_date'] else 'Unknown'
                    print(f"  🧬 {record['discovery_id'][:12]}... | "
                          f"Length: {record['sequence_length']:3d} | "
                          f"Score: {record['validation_score']:.3f} | "
                          f"Energy: {record['energy']:8.2f} | "
                          f"Coherence: {record['quantum_coherence']:.3f} | "
                          f"Time: {timestamp}")
                
                print("-" * 80)
                
                # Get total count for comparison
                total_query = """
                MATCH (d:Discovery)
                RETURN count(d) as total_count
                """
                total_result = session.run(total_query)
                total_count = total_result.single()['total_count']
                
                print(f"\n📈 SUMMARY:")
                print(f"   Total discoveries in database: {total_count:,}")
                print(f"   New discoveries since export: {new_count:,}")
                print(f"   Percentage increase: {(new_count/total_count)*100:.2f}%")
                
                # Check if significant enough to warrant update
                if new_count > 1000 or (new_count/total_count) > 0.01:  # >1K new or >1% increase
                    print(f"\n✅ SIGNIFICANT UPDATE DETECTED!")
                    print(f"   Recommendation: Update Streamlit data and GitHub wiki")
                    return True, new_count, total_count
                else:
                    print(f"\n📝 Minor update ({new_count} discoveries)")
                    print(f"   Recommendation: Update not critical, but can proceed if desired")
                    return True, new_count, total_count
                    
            else:
                print("✅ No new discoveries since last export")
                print("📊 Data is up to date!")
                
                # Still get total count for info
                total_query = """
                MATCH (d:Discovery)
                RETURN count(d) as total_count
                """
                total_result = session.run(total_query)
                total_count = total_result.single()['total_count']
                
                print(f"📈 Current total discoveries: {total_count:,}")
                return False, 0, total_count
                
    except Exception as e:
        print(f"❌ Error checking for new discoveries: {str(e)}")
        import traceback
        traceback.print_exc()
        return False, 0, 0
    
    finally:
        try:
            engine.close()
        except:
            pass

def get_current_stats():
    """Get current discovery statistics for comparison"""
    
    print("\n📊 CURRENT DISCOVERY STATISTICS")
    print("=" * 50)
    
    try:
        engine = Neo4jDiscoveryEngine()
        
        with engine.driver.session() as session:
            # Get comprehensive stats
            stats_query = """
            MATCH (d:Discovery)-[:HAS_SEQUENCE]->(s:Sequence)
            RETURN count(d) as total_discoveries,
                   avg(d.validation_score) as avg_validation_score,
                   avg(d.quantum_coherence) as avg_quantum_coherence,
                   avg(s.length) as avg_sequence_length,
                   min(d.timestamp) as first_discovery,
                   max(d.timestamp) as latest_discovery
            """
            
            result = session.run(stats_query)
            record = result.single()
            
            if record:
                print(f"🔢 Total discoveries: {record['total_discoveries']:,}")
                print(f"⭐ Avg validation score: {record['avg_validation_score']:.4f}")
                print(f"⚛️  Avg quantum coherence: {record['avg_quantum_coherence']:.4f}")
                print(f"📏 Avg sequence length: {record['avg_sequence_length']:.1f}")
                
                if record['first_discovery']:
                    first = record['first_discovery'].strftime('%Y-%m-%d %H:%M:%S')
                    print(f"🎯 First discovery: {first}")
                
                if record['latest_discovery']:
                    latest = record['latest_discovery'].strftime('%Y-%m-%d %H:%M:%S')
                    print(f"🕐 Latest discovery: {latest}")
                
                # Get quality distribution
                quality_query = """
                MATCH (d:Discovery)
                RETURN 
                    count(CASE WHEN d.validation_score >= 0.9 THEN 1 END) as excellent,
                    count(CASE WHEN d.validation_score >= 0.8 AND d.validation_score < 0.9 THEN 1 END) as very_good,
                    count(CASE WHEN d.validation_score >= 0.7 AND d.validation_score < 0.8 THEN 1 END) as good,
                    count(CASE WHEN d.validation_score < 0.7 THEN 1 END) as lower_quality
                """
                
                quality_result = session.run(quality_query)
                quality_record = quality_result.single()
                
                if quality_record:
                    total = record['total_discoveries']
                    print(f"\n📈 QUALITY DISTRIBUTION:")
                    print(f"   🌟 Excellent (≥0.9): {quality_record['excellent']:,} ({quality_record['excellent']/total*100:.1f}%)")
                    print(f"   ⭐ Very Good (0.8-0.9): {quality_record['very_good']:,} ({quality_record['very_good']/total*100:.1f}%)")
                    print(f"   ✅ Good (0.7-0.8): {quality_record['good']:,} ({quality_record['good']/total*100:.1f}%)")
                    print(f"   📊 Lower (<0.7): {quality_record['lower_quality']:,} ({quality_record['lower_quality']/total*100:.1f}%)")
            else:
                print("❌ Could not retrieve discovery statistics")
                
    except Exception as e:
        print(f"❌ Error getting current stats: {str(e)}")
    
    finally:
        try:
            engine.close()
        except:
            pass

if __name__ == "__main__":
    print("🔍 NEW DISCOVERY CHECK & STATISTICS")
    print("=" * 60)
    
    # Check for new discoveries
    has_new, new_count, total_count = check_for_new_discoveries()
    
    # Get current statistics
    get_current_stats()
    
    print("\n🎯 RECOMMENDATION:")
    if has_new and new_count > 0:
        print("✅ UPDATE RECOMMENDED")
        print("   Actions:")
        print("   1. Run export_complete_dataset.py to update data files")
        print("   2. Run export_chunked_for_git.py for GitHub-compatible chunks")
        print("   3. Commit and push updated data to GitHub")
        print("   4. Update GitHub wiki with new discovery count")
        print("   5. Streamlit Cloud will automatically pick up changes")
    else:
        print("✅ NO UPDATE NEEDED")
        print("   Current data is up to date!")
    
    print(f"\n📊 Status: {total_count:,} total discoveries in database")
    print("🎉 Discovery check complete!")
