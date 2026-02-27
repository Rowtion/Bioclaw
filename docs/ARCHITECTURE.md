# Architecture

## System Overview

AIpoch-claw integrates multiple components to create a seamless biomedical research assistant.

## Component Interaction

```
┌─────────────────────────────────────────────────────────────────┐
│                         User Layer                               │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────────────┐    │
│  │  Feishu │  │ WhatsApp│  │  Slack  │  │     Discord     │    │
│  └────┬────┘  └────┬────┘  └────┬────┘  └────────┬────────┘    │
└───────┼────────────┼────────────┼────────────────┼─────────────┘
        │            │            │                │
        └────────────┴────────────┘────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                      OpenClaw Gateway                            │
│                   (Message Routing)                              │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                   opencode-bridge Skill                          │
│              (OpenClaw → Opencode Bridge)                        │
│                                                                  │
│  • Session Management    • Message Forwarding                    │
│  • Error Handling        • Response Streaming                    │
└──────────────────────────┬──────────────────────────────────────┘
                           │ HTTP API
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Opencode Server                              │
│                     (Port 4096)                                  │
│                                                                  │
│  • Code Execution        • Skill Resolution                      │
│  • Context Management    • Tool Orchestration                    │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│              AIPOCH Medical Research Skills                      │
│                      (200+ Skills)                               │
│                                                                  │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────┐     │
│  │   Clinical   │ │Bioinformatics│ │   Data Analysis      │     │
│  │   Research   │ │              │ │                      │     │
│  └──────────────┘ └──────────────┘ └──────────────────────┘     │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────┐     │
│  │  Literature  │ │Pharmaceutical│ │   Experimental       │     │
│  │  & Evidence  │ │              │ │   Research           │     │
│  └──────────────┘ └──────────────┘ └──────────────────────┘     │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                 Docker Analysis Environment                      │
│                                                                  │
│  ┌─────────────────────┐    ┌─────────────────────┐             │
│  │   RStudio Server    │    │     JupyterLab      │             │
│  │      (Port 8787)    │    │      (Port 8888)    │             │
│  │                     │    │                     │             │
│  │  • R 4.3.3          │    │  • Python 3         │             │
│  │  • Bioconductor     │    │  • scanpy           │             │
│  │  • DESeq2, Seurat   │    │  • biopython        │             │
│  │  • tidyverse        │    │  • scikit-learn     │             │
│  └─────────────────────┘    └─────────────────────┘             │
└─────────────────────────────────────────────────────────────────┘
```

## Data Flow

1. **User Input**: User sends message via messaging platform
2. **OpenClaw Routing**: Message routed through OpenClaw gateway
3. **Skill Matching**: `opencode-bridge` skill matches trigger words
4. **API Forwarding**: Bridge forwards message to Opencode HTTP API
5. **Skill Execution**: Opencode loads and executes appropriate AIPOCH skill
6. **Code Execution**: Skill runs analysis code in Docker environment
7. **Result Return**: Results streamed back through the chain to user

## Security Considerations

- Opencode runs locally on user's machine (localhost:4096)
- Docker containers are isolated from host system
- API keys stored in `.env` file (not committed to git)
- No data leaves local network unless explicitly configured

## Scalability

Current architecture supports:
- Single-user local deployment
- Multiple concurrent sessions
- Extensible skill library
- Pluggable messaging platforms
