# AIpoch-claw Project Context

This file provides context for AI agents working on this project.

## Project Overview

AIpoch-claw is a biomedical AI research assistant that combines:
- **OpenClaw**: Conversational AI gateway for messaging platforms
- **Opencode**: Code execution environment (alternative to Claude Code)
- **AIPOCH Medical Research Skills**: 200+ curated medical research skills

## Architecture

```
User (Feishu/WhatsApp/Slack/Discord)
    ↓
OpenClaw Gateway
    ↓
opencode-bridge skill
    ↓ HTTP API
Opencode Server (localhost:4096)
    ↓
AIPOCH Medical Research Skills (200+)
    ↓
Docker Analysis Environment (R + Python)
    ├─ RStudio Server (:8787)
    └─ JupyterLab (:8888)
```

## Key Components

### 1. Bridge Skill (`skills/opencode-bridge/`)
- Located in `~/.openclaw/workspace/skills/opencode-bridge/`
- Forwards OpenClaw messages to Opencode via HTTP API
- Handles session management and response streaming

### 2. Analysis Environment (`docker/`)
- R 4.3.3 with Bioconductor packages (DESeq2, Seurat, etc.)
- Python 3 with scientific packages (scanpy, biopython, etc.)
- RStudio Server and JupyterLab for interactive analysis

### 3. Skill Library (`scientific-skills/`)
- Downloaded from https://github.com/aipoch/medical-research-skills
- 200+ skills across clinical research, bioinformatics, data analysis
- Downloaded during `setup.sh` execution

## Development Guidelines

### Adding New Skills
1. Create skill directory under appropriate category
2. Write SKILL.md following AIPOCH format
3. Add to `scientific-skills/` directory
4. Test with Opencode directly

### Modifying Bridge
- Core logic in `skills/opencode-bridge/scripts/bridge.py`
- Supports session persistence and error handling
- Environment variables control behavior

### Docker Updates
- Modify `docker/Dockerfile` for new packages
- Rebuild with `docker-compose build --no-cache`
- Data persists in `data/` and `outputs/` volumes

## Common Tasks

### Start Development Environment
```bash
docker-compose up -d
~/.opencode/bin/opencode serve --port 4096
```

### Test Bridge Connection
```bash
python skills/opencode-bridge/scripts/bridge.py "test message"
```

### Update AIPOCH Skills
```bash
cd scientific-skills
git pull
cd ..
```

## References

- OpenClaw: https://github.com/openclaw/openclaw
- Opencode: https://opencode.dev
- AIPOCH: https://AIPOCH.com
- AIPOCH Skills: https://github.com/aipoch/medical-research-skills
