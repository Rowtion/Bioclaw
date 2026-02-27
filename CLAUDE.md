# Bioclaw Project Context

This file provides context for AI agents working on this project.

## Project Overview

Bioclaw is a scientific AI research assistant that combines:
- **OpenClaw**: Conversational AI gateway for messaging platforms
- **Opencode**: Code execution environment (alternative to Claude Code)
- **K-Dense Scientific Skills**: 140+ curated scientific computing skills

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
K-Dense Scientific Skills (140+)
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
- Downloaded from https://github.com/K-Dense-AI/claude-scientific-skills
- 140+ skills across bioinformatics, data analysis, literature search, lab tools
- Downloaded during `setup.sh` execution

## Development Guidelines

### Adding New Skills
1. Create skill directory under appropriate category
2. Write SKILL.md following K-Dense format
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

### Update K-Dense Skills
```bash
cd scientific-skills
git pull
cd ..
```

## References

- OpenClaw: https://github.com/openclaw/openclaw
- Opencode: https://opencode.dev
- K-Dense: https://k-dense.ai
- K-Dense Skills: https://github.com/K-Dense-AI/claude-scientific-skills
