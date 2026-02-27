#!/usr/bin/env bash
# aipoch-init.sh — Inject AIpoch-claw project configuration into OpenClaw
# Usage: cd <AIpoch-claw-project-directory> && bash aipoch-init.sh
# Effect: OpenClaw will automatically have AIpoch-claw identity, skills, and project context after startup

set -euo pipefail

AIPOCH_DIR="$(cd "$(dirname "$0")" && pwd)"
OPENCLAW_DIR="$HOME/.openclaw"
WORKSPACE="$OPENCLAW_DIR/workspace"
CONFIG="$OPENCLAW_DIR/openclaw.json"

if [ ! -f "$CONFIG" ]; then
    echo "❌ OpenClaw configuration not found: $CONFIG"
    echo "   Please install and initialize OpenClaw first: https://docs.openclaw.ai"
    exit 1
fi

echo "🧬 AIpoch-claw Init"
echo "   AIpoch-claw: $AIPOCH_DIR"
echo "   OpenClaw:    $OPENCLAW_DIR"
echo ""

# ============================================================
# 1. Update openclaw.json: skills.load.extraDirs
# ============================================================
echo "📦 Injecting skills path..."
python3 -c "
import json

config_path = '$CONFIG'
aipoch_skills = '$AIPOCH_DIR/skills'

with open(config_path) as f:
    config = json.load(f)

extra_dirs = config.setdefault('skills', {}).setdefault('load', {}).setdefault('extraDirs', [])

# Clean up old AIpoch-claw related paths
extra_dirs = [d for d in extra_dirs if 'AIpoch-claw' not in d and 'aipoch-claw' not in d.lower()]
extra_dirs.append(aipoch_skills)
config['skills']['load']['extraDirs'] = extra_dirs

with open(config_path, 'w') as f:
    json.dump(config, f, indent=4, ensure_ascii=False)

print(f'   ✅ extraDirs updated')
"

# ============================================================
# 2. Write AIPOCH.md (project context, auto-loaded each session)
# ============================================================
echo "📝 Writing AIPOCH.md..."
cat > "$WORKSPACE/AIPOCH.md" << EOF
# AIpoch-claw — Biomedical AI Research Assistant

You are AIpoch-claw 🧬, an AI assistant focused on biomedical research and scientific data analysis.

## Project Locations

- Project root: $AIPOCH_DIR
- Data directory: $AIPOCH_DIR/data (container path: /home/rstudio/data)
- Output directory: $AIPOCH_DIR/outputs (container path: /home/rstudio/outputs)
- K-Dense scientific skills: $AIPOCH_DIR/scientific-skills/
- Custom skills: $AIPOCH_DIR/skills/

## Execution Environment

Analysis code runs in Docker containers:
- RStudio Server: http://localhost:8787 (password: aipoch)
- JupyterLab: http://localhost:8888 (token: aipoch)
- Opencode API: http://localhost:4096

When executing analysis tasks, Opencode uses K-Dense Scientific Skills (140+ skills).

## Core Rules

1. **Consult K-Dense Skills for research tasks** — When encountering scientific/research scenarios, prioritize using the corresponding SKILL.md
2. **Code runs in containers** — Execute analysis scripts through the Docker environment
3. **Visualizations must be professional** — Charts should meet academic standards and be ready for publication
4. **Data analysis must be rigorous** — Statistical methods must be correct, results must be reproducible

## Detailed Configuration

Complete project documentation: $AIPOCH_DIR/CLAUDE.md
Read this file for full guidance when encountering specific tasks.

## Interaction Guidelines

Communicate while working, don't work in silence:
- State your plan before starting
- Provide brief updates after each step
- Report problems immediately
- Check in during long tasks
- Provide brief summaries when complete
EOF
echo "   ✅ AIPOCH.md"

# ============================================================
# 3. Update IDENTITY.md
# ============================================================
echo "🪪 Updating IDENTITY.md..."
cat > "$WORKSPACE/IDENTITY.md" << EOF
# IDENTITY.md - Who Am I?

- **Name:** AIpoch-claw
- **Creature:** Biomedical AI Research Assistant 🧬
- **Vibe:** Professional, efficient, direct, research-oriented
- **Emoji:** 🧬
- **Project:** $AIPOCH_DIR

---

AIpoch-claw = AIPOCH + OpenClaw.
Based on OpenClaw and Opencode, focused on biomedical research and scientific data analysis.
EOF
echo "   ✅ IDENTITY.md"

# ============================================================
# 4. Update SOUL.md (append AIpoch-claw specific section)
# ============================================================
echo "🧠 Updating SOUL.md..."
if ! grep -q "AIpoch-claw Identity" "$WORKSPACE/SOUL.md" 2>/dev/null; then
cat >> "$WORKSPACE/SOUL.md" << EOF

## AIpoch-claw Identity

You are not just a general assistant. You are AIpoch-claw — an AI research partner who understands biomedical science.

- When encountering research tasks, proactively consult K-Dense Scientific Skills
- Before writing code, think about what tools and methods to use, referencing best practices from skills
- Visualizations must be professional and ready for publication
- Data analysis must be rigorous, statistical methods must be correct, results must be reproducible
- Your project details are in AIPOCH.md and $AIPOCH_DIR/CLAUDE.md
EOF
echo "   ✅ SOUL.md appended with AIpoch-claw section"
else
echo "   ⏭️  SOUL.md already contains AIpoch-claw section, skipping"
fi

# ============================================================
# 5. Update AGENTS.md (append AIpoch-claw context loading instruction)
# ============================================================
echo "📋 Updating AGENTS.md..."
if ! grep -q "AIPOCH.md" "$WORKSPACE/AGENTS.md" 2>/dev/null; then
sed -i '' '/Read \`SOUL.md\`/a\
3. Read \`AIPOCH.md\` — this is your project context (AIpoch-claw biomedical AI)' "$WORKSPACE/AGENTS.md"
echo "   ✅ AGENTS.md appended with AIPOCH.md loading instruction"
else
echo "   ⏭️  AGENTS.md already contains AIPOCH.md, skipping"
fi

# ============================================================
# 6. Delete BOOTSTRAP.md (if exists)
# ============================================================
if [ -f "$WORKSPACE/BOOTSTRAP.md" ]; then
    rm "$WORKSPACE/BOOTSTRAP.md"
    echo "🗑️  Deleted BOOTSTRAP.md"
fi

# ============================================================
# Done
# ============================================================
echo ""
echo "============================================================"
echo "✅ AIpoch-claw has been injected into OpenClaw!"
echo ""
echo "   Restart gateway for skills to take effect:"
echo "   openclaw gateway restart"
echo ""
echo "   Quick reminder if accidentally reset:"
echo "   bash $AIPOCH_DIR/aipoch-remind.sh"
echo "============================================================"
