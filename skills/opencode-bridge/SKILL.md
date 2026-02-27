---
name: opencode-bridge
description: Bridge skill that forwards messages to local Opencode server (port 4096) for execution. Use when user wants to run scientific analysis, bioinformatics workflows, data processing, or any task requiring specialized scientific computing skills. Triggers: "use opencode", "opencode analysis", "analyze data", "run analysis", "bioinformatics", "data processing".
triggers:
  - "use opencode"
  - "opencode analysis"
  - "medical analysis"
  - "bioinformatics"
  - "clinical analysis"
  - "research analysis"
  - "run analysis"
  - "analyze data"
  - "perform research"
  - "using opencode"
---

# Opencode Bridge Skill

Forwards OpenClaw messages to local Opencode service for execution.

## When to Use This Skill

Trigger when users need to perform the following operations:

- **Medical Data Analysis**: "Help me analyze this clinical data"
- **Bioinformatics Workflows**: "Use opencode for single-cell analysis"
- **Literature Review**: "Search and summarize recent CRISPR papers"
- **Research Design**: "Help me design a randomized controlled trial"
- **Statistical Computing**: "Calculate sample size", "Survival analysis"
- **Drug Discovery**: "Virtual screening", "ADMET prediction"

## How to Use

### 1. Ensure Opencode Service is Running

```bash
# Check service status
curl http://localhost:4096/status

# If not running, start it
~/.opencode/bin/opencode serve --port 4096
```

### 2. Use the Skill

After the user sends a message, the skill will automatically forward it to Opencode.

### 3. View Results

- Short tasks: Returned directly in the conversation
- Long tasks: Returns session ID, user can view progress at http://localhost:4096

## Available Skills (K-Dense Scientific)

After installation, Opencode will have 140+ scientific computing skills, including:

| Category | Example Skills |
|----------|---------------|
| 🔬 Bioinformatics | Genomics, transcriptomics, proteomics, sequence analysis |
| 📊 Data Analysis | Statistical analysis, visualization, machine learning |
| 📚 Literature & Search | Literature search, paper analysis, citation management |
| 🧪 Lab Tools | Protocol design, reagent calculation, experiment tracking |
| 💻 Programming | Python, R, data processing, automation scripts |
| 🗄️ Databases | PubMed, ChEMBL, ClinicalTrials.gov integration |

## Environment Variables

```bash
OPENCODE_URL=http://localhost:4096  # Opencode service address
OPENCODE_TIMEOUT=300                 # Request timeout (seconds)
```

## Scripts

- `scripts/bridge.py` - Core bridge script

## Architecture

```
Feishu/OpenClaw
    ↓
opencode-bridge skill (this)
    ↓ HTTP API
Opencode Server (:4096)
    ↓
K-Dense Scientific Skills (140+)
    ↓
Docker Analysis Environment
```
