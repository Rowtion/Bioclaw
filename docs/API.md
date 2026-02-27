# API Reference

## Opencode Bridge API

The bridge skill communicates with Opencode via HTTP API.

### Endpoints

#### List Sessions
```http
GET /session
```

Returns all active Opencode sessions.

#### Create Session
```http
POST /session
Content-Type: application/json

{
  "title": "New Session"
}
```

#### Send Message
```http
POST /session/{session_id}/message
Content-Type: application/json

{
  "parts": [
    {"type": "text", "text": "your message here"}
  ]
}
```

#### Get Messages
```http
GET /session/{session_id}/message?limit=10
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `OPENCODE_URL` | `http://localhost:4096` | Opencode server URL |
| `OPENCODE_TIMEOUT` | `300` | Request timeout in seconds |

## Docker Services

### RStudio Server
- **Port**: 8787
- **Default Password**: `aipoch`
- **Data Directory**: `/home/rstudio/data`

### JupyterLab
- **Port**: 8888
- **Default Token**: `aipoch`
- **Notebook Directory**: `/home/rstudio/data`

## Skill Trigger Words

The bridge skill activates on these keywords:

- `用opencode`
- `opencode分析`
- `医学分析`
- `生物信息学`
- `clinical analysis`
- `research analysis`
- `run analysis`
- `analyze data`
- `perform research`
