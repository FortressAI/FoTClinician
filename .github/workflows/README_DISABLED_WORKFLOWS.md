# Disabled GitHub Actions Workflows

## discovery.yml.disabled

**Status**: DISABLED (renamed from discovery.yml)  
**Reason**: Architecture mismatch with actual FoTChemistry implementation

### The Problem:
- **Workflow Expected**: Docker/AKG setup with orchestrator.py
- **Actual System**: Neo4j + Python discovery engines (continuous_chemistry_discovery.py)
- **Result**: Hourly failures causing recurring notifications

### What It Was Trying To Do:
```yaml
schedule:
  - cron: "0 * * * *"    # Every hour
  - cron: "0 6 * * *"    # Daily  
  - cron: "0 12 * * 1"   # Weekly
```

- Run autonomous discovery campaigns
- Use orchestrator.py with campaign YAML configs
- Expected Docker containerized AKG setup
- Truth-mining workflows with different architecture

### Our Actual Architecture:
- **Database**: Local Neo4j (not Docker AKG)
- **Discovery Engine**: `continuous_chemistry_discovery.py`
- **Real Molecular Generation**: `agents/alchemist/real_molecular_generator.py`
- **Quantum Substrate**: `core/chemistry_vqbit_engine.py`
- **Results**: 3,043+ real molecular discoveries

### To Re-enable (Future):
1. Adapt workflow to use our actual Python discovery system
2. Point to `continuous_chemistry_discovery.py` instead of `orchestrator.py`
3. Update to work with Neo4j connection (not Docker AKG)
4. Test manually before enabling schedule

### Alternative:
The current system runs locally and produces real results. GitHub Actions automation is optional.
