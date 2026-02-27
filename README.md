<p align="center">
  <img src="https://via.placeholder.com/200x200?text=AIpoch-claw" alt="AIpoch-claw Logo" width="200">
</p>

<h1 align="center">AIpoch-claw</h1>

<p align="center">
  <strong>Open-Source Biomedical AI Research Assistant</strong><br>
  OpenClaw × Opencode × AIPOCH Medical Research Skills
</p>

<p align="center">
  <a href="https://github.com/Rowtion/AIpoch-claw/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT">
  </a>
  <a href="https://github.com/Rowtion/AIpoch-claw">
    <img src="https://img.shields.io/badge/Skills-200+-green.svg" alt="Skills: 200+">
  </a>
  <a href="https://opencode.dev">
    <img src="https://img.shields.io/badge/Powered%20by-Opencode-orange.svg" alt="Powered by Opencode">
  </a>
</p>

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Available Skills](#available-skills)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Overview

**AIpoch-claw** is an open-source biomedical AI research assistant that enables researchers to perform complex medical data analysis through natural language conversations. It integrates three powerful components:

- **OpenClaw**: Conversational AI gateway that connects to your favorite messaging platforms (Feishu, WhatsApp, Slack, Discord)
- **Opencode**: Lightweight, self-hosted code execution environment
- **AIPOCH Medical Research Skills**: A curated collection of 200+ medical research skills covering clinical research, bioinformatics, data analysis, and more

### Why AIpoch-claw?

| Feature | Benefit |
|---------|---------|
| 🔬 **Research-Ready** | 200+ pre-built skills for common biomedical tasks |
| 💬 **Conversational** | Interact via messaging apps you already use |
| 🔒 **Self-Hosted** | Your data stays on your machine |
| 🚀 **Extensible** | Easy to add custom skills |
| 🐳 **Containerized** | Reproducible analysis environment via Docker |

## Features

### Core Capabilities

- **Clinical Research**: Trial design, patient data analysis, diagnostic support
- **Bioinformatics**: Genomics, transcriptomics, proteomics workflows
- **Data Analysis**: Statistical analysis, visualization, machine learning
- **Literature Synthesis**: Literature search, systematic reviews, evidence grading
- **Pharmaceutical Research**: Drug discovery, ADMET prediction, target analysis

### Supported Platforms

- Feishu (飞书)
- WhatsApp
- Slack
- Discord
- And more via OpenClaw plugins

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Layer                               │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────────────┐    │
│  │  Feishu │  │ WhatsApp│  │  Slack  │  │     Discord     │    │
│  └────┬────┘  └────┬────┘  └────┬────┘  └────────┬────────┘    │
└───────┼────────────┼────────────┼────────────────┼─────────────┘
        │            │            │                │
        └────────────┴────────────┴────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                      OpenClaw Gateway                            │
│                   (Message Routing Layer)                        │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                   opencode-bridge Skill                          │
│              (OpenClaw → Opencode Integration)                   │
└──────────────────────────┬──────────────────────────────────────┘
                           │ HTTP API
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Opencode Server                              │
│              (localhost:4096, Code Execution)                    │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│              AIPOCH Medical Research Skills                      │
│                      (200+ Skills)                               │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                 Docker Analysis Environment                      │
│                                                                  │
│  ┌─────────────────────┐    ┌─────────────────────┐             │
│  │   RStudio Server    │    │     JupyterLab      │             │
│  │      (:8787)        │    │      (:8888)        │             │
│  │                     │    │                     │             │
│  │  • R 4.3.3          │    │  • Python 3         │             │
│  │  • Bioconductor     │    │  • scanpy           │             │
│  │  • DESeq2, Seurat   │    │  • biopython        │             │
│  │  • tidyverse        │    │  • scikit-learn     │             │
│  └─────────────────────┘    └─────────────────────┘             │
└─────────────────────────────────────────────────────────────────┘
```

## Quick Start

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) and docker-compose
- [Git](https://git-scm.com/)
- [Node.js](https://nodejs.org/) 22+ (for Opencode)
- [OpenClaw](https://github.com/openclaw/openclaw) (installed and configured)

### 1. Install Opencode

```bash
# Install Opencode
curl -fsSL https://opencode.dev/install.sh | sh

# Start the server
~/.opencode/bin/opencode serve --port 4096
```

### 2. Clone and Setup AIpoch-claw

```bash
# Clone the repository
git clone https://github.com/Rowtion/AIpoch-claw.git
cd AIpoch-claw

# Copy environment template
cp .env.template .env

# (Optional) Edit configuration
nano .env

# Run setup script
bash setup.sh
```

The setup script will:
- Download AIPOCH Medical Research Skills (200+ skills)
- Build Docker images with R and Python environments
- Start RStudio Server and JupyterLab
- Create necessary directories

### 3. Install the Bridge Skill

```bash
# Copy bridge skill to OpenClaw skills directory
cp -r skills/opencode-bridge ~/.openclaw/workspace/skills/

# Verify installation
openclaw skills check
```

### 4. Start Using

Send messages through your connected messaging platform:

**Clinical Analysis:**
```
Analyze clinical trial data comparing treatment vs control groups
```

**Bioinformatics:**
```
Perform single-cell RNA-seq analysis on my 10X data
```

**Statistics:**
```
Calculate sample size for a two-arm RCT with 80% power
```

**Literature Research:**
```
Search PubMed for recent CRISPR base editing papers and summarize top 10
```

Results are saved to `./outputs/` and can be viewed in:
- **RStudio Server**: http://localhost:8787 (password: `aipoch`)
- **JupyterLab**: http://localhost:8888 (token: `aipoch`)

## Configuration

Create a `.env` file in the project root:

```bash
# Opencode connection
OPENCODE_URL=http://localhost:4096
OPENCODE_TIMEOUT=300

# Model configuration (if using third-party APIs)
# ANTHROPIC_API_KEY=your_key_here
# ANTHROPIC_SMALL_FAST_MODEL=claude-sonnet-4-20250514

# Docker services
JUPYTER_TOKEN=aipoch
RSTUDIO_PASSWORD=aipoch
```

### Using Third-Party API Providers

If using MiniMax, GLM, DeepSeek, or other non-Anthropic endpoints:

```bash
# Required: Set a compatible model for pre-flight checks
ANTHROPIC_SMALL_FAST_MODEL=claude-sonnet-4-20250514
```

## Available Skills

AIPOCH Medical Research Skills includes 200+ skills organized by domain:

| Category | Skills | Examples |
|----------|--------|----------|
| 🏥 **Clinical Research** | 20+ | Trial design, patient data analysis, diagnostic support |
| 🔬 **Experimental Research** | 15+ | Protocol design, quality control, experiment optimization |
| 📊 **Medical Data Analysis** | 30+ | Statistical analysis, visualization, machine learning |
| 🧬 **Bioinformatics** | 25+ | Genomics, transcriptomics, proteomics, metabolomics |
| 📚 **Literature & Evidence** | 20+ | Literature search, systematic review, meta-analysis |
| 💊 **Pharmaceutical** | 15+ | Drug discovery, ADMET prediction, target analysis |
| 🎓 **Education** | 11+ | Teaching materials, training workflows |
| 🧾 **Grant & Strategy** | 10+ | Grant writing, research strategy, career development |

Browse all skills at: https://github.com/aipoch/medical-research-skills

## Troubleshooting

### Opencode Connection Issues

```bash
# Check if Opencode is running
curl http://localhost:4096/status

# View Opencode logs
~/.opencode/bin/opencode logs

# Restart Opencode
~/.opencode/bin/opencode serve --port 4096
```

### Docker Services Won't Start

```bash
# View container logs
docker-compose logs

# Rebuild images (no cache)
docker-compose down
docker-compose build --no-cache
docker-compose up -d

# Check container status
docker-compose ps
```

### Skills Not Triggering

1. Verify bridge skill is installed:
   ```bash
   ls ~/.openclaw/workspace/skills/opencode-bridge/
   ```

2. Check OpenClaw skill recognition:
   ```bash
   openclaw skills check
   ```

3. Ensure trigger words are in your message:
   - `用opencode`, `opencode分析`, `医学分析`, `生物信息学`
   - `clinical analysis`, `research analysis`, `run analysis`

### API Errors with Third-Party Providers

If you see `Pre-flight check is taking longer than expected`:

1. Add `ANTHROPIC_SMALL_FAST_MODEL` to `.env`
2. Re-run `bash setup.sh`
3. Test: `claude --dangerously-skip-permissions -p 'run: echo hello'`

## Project Structure

```
AIpoch-claw/
├── docker/
│   ├── Dockerfile              # R + Python analysis environment
│   └── entrypoint.sh           # Service startup script
│
├── skills/
│   └── opencode-bridge/        # OpenClaw → Opencode bridge
│       ├── SKILL.md
│       └── scripts/
│           └── bridge.py       # Bridge implementation
│
├── docs/
│   ├── API.md                  # API documentation
│   ├── ARCHITECTURE.md         # System architecture details
│   └── QUICKSTART_ZH.md        # 中文快速开始
│
├── data/                       # Data directory (mounted to containers)
├── outputs/                    # Analysis outputs (mounted to containers)
│
├── docker-compose.yml          # Docker services configuration
├── setup.sh                    # Installation script
├── .env.template               # Environment variables template
├── .gitignore
├── CLAUDE.md                   # AI agent context
├── CONTRIBUTING.md             # Contribution guidelines
├── LICENSE                     # MIT License
└── README.md                   # This file
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Ways to contribute:
- 🐛 Report bugs
- 💡 Suggest new features
- 📝 Improve documentation
- 🔧 Submit bug fixes
- ✨ Add new medical research skills

## Comparison with Similar Projects

| Aspect | MedgeClaw | AIpoch-claw |
|--------|-----------|-------------|
| Execution Engine | Claude Code | **Opencode** (self-hosted) |
| Skill Library | K-Dense (140 skills) | **AIPOCH (200+ skills)** |
| Gateway | OpenClaw | **OpenClaw** |
| Analysis Environment | Docker | **Docker** |
| Messaging Platforms | WhatsApp/Slack/Discord | **Feishu/WhatsApp/Slack/Discord** |
| API Flexibility | Anthropic only | **Multiple providers supported** |

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgments

- [OpenClaw](https://github.com/openclaw/openclaw) - Conversational AI gateway
- [Opencode](https://opencode.dev) - Code execution environment
- [AIPOCH](https://AIPOCH.com) - Medical research skills ecosystem
- [AIPOCH Medical Research Skills](https://github.com/aipoch/medical-research-skills) - Skill library

## Links

- 🏠 **Repository**: https://github.com/Rowtion/AIpoch-claw
- 📖 **Documentation**: See `docs/` directory
- 🐛 **Issues**: https://github.com/Rowtion/AIpoch-claw/issues
- 💬 **Discussions**: https://github.com/Rowtion/AIpoch-claw/discussions

---

<p align="center">
  Made with ❤️ for biomedical researchers worldwide
</p>
