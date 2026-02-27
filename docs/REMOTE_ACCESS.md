# Remote Access Guide

Enable remote access to AIpoch-claw services (RStudio, JupyterLab, Opencode) from anywhere.

## Quick Start

### 1. Configure Environment

```bash
cp .env.remote.template .env

# Edit .env and set secure passwords
nano .env
```

**Important:** Change default passwords for security!

### 2. Start Services with Remote Access

```bash
# Use the remote docker-compose file
docker-compose -f docker-compose-remote.yml up -d

# Or rename it as the main file:
mv docker-compose.yml docker-compose-local.yml
mv docker-compose-remote.yml docker-compose.yml
docker-compose up -d
```

### 3. Access Services

Once running, access from any device on your network:

| Service | Local URL | Remote URL (replace YOUR_IP) |
|---------|-----------|------------------------------|
| **RStudio Server** | http://localhost:8787 | http://YOUR_IP:8787 |
| **JupyterLab** | http://localhost:8888 | http://YOUR_IP:8888 |
| **Opencode** | http://localhost:4096 | http://YOUR_IP:4096 |

**Find your IP address:**
```bash
# macOS
ifconfig | grep "inet " | grep -v 127.0.0.1

# Linux
ip addr show | grep "inet " | grep -v 127.0.0.1
```

## Security Considerations

⚠️ **WARNING:** Opening ports to the internet has security risks!

### Basic Security (Minimum)

1. **Change default passwords** in `.env`:
```bash
RSTUDIO_PASSWORD=your_secure_password_123
JUPYTER_TOKEN=your_secure_token_456
```

2. **Use firewall** to restrict access:
```bash
# macOS - Only allow specific IP
sudo ipfw add allow tcp from 192.168.1.0/24 to any 8787,8888,4096

# Linux (ufw)
sudo ufw allow from 192.168.1.0/24 to any port 8788,8888,4096
```

### Advanced Security (Recommended for Internet)

#### Option A: SSH Tunnel (Safest)

Don't expose ports directly. Use SSH tunnel:

```bash
# From your local machine, create tunnels to remote server
ssh -L 8787:localhost:8787 -L 8888:localhost:8888 -L 4096:localhost:4096 user@your-server-ip

# Then access via localhost on your machine
# http://localhost:8787
# http://localhost:8888
# http://localhost:4096
```

#### Option B: Nginx Reverse Proxy with SSL

Use Nginx as a reverse proxy with HTTPS:

```nginx
server {
    listen 443 ssl;
    server_name your-domain.com;
    
    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;
    
    location /rstudio/ {
        proxy_pass http://localhost:8787/;
        proxy_set_header Host $host;
    }
    
    location /jupyter/ {
        proxy_pass http://localhost:8888/;
        proxy_set_header Host $host;
    }
}
```

#### Option C: VPN (Best for Teams)

Use Tailscale, ZeroTier, or WireGuard for secure private network:

```bash
# Install Tailscale
curl -fsSL https://tailscale.com/install.sh | sh
sudo tailscale up

# Now all devices in your Tailscale network can access
# http://your-machine-tailscale-ip:8787
```

## Cloud Deployment

### AWS/GCP/Azure

1. Open security group ports: 8787, 8888, 4096
2. Use strong passwords
3. Consider using AWS Systems Manager Session Manager instead of opening ports

### Docker with Cloudflared (Free Tunnel)

```yaml
version: '3.8'
services:
  tunnel:
    image: cloudflare/cloudflared:latest
    command: tunnel --no-autoupdate run --token YOUR_TOKEN
    restart: unless-stopped
    
  analysis-env:
    # ... same as before
```

## Troubleshooting

### Cannot access from remote

1. **Check firewall:**
```bash
# macOS
sudo /usr/libexec/ApplicationFirewall/socketfilterfw --getglobalstate

# Linux
sudo iptables -L | grep 8787
```

2. **Check service is listening on 0.0.0.0:**
```bash
netstat -an | grep 8787
# Should show 0.0.0.0:8787, not 127.0.0.1:8787
```

3. **Test from another machine:**
```bash
curl http://your-server-ip:8787
```

### Connection refused

- Ensure Docker containers are running: `docker-compose ps`
- Check logs: `docker-compose logs analysis-env`
- Verify port mapping: `docker port aipoch-analysis`

## Production Checklist

Before exposing to internet:

- [ ] Changed all default passwords
- [ ] Enabled HTTPS (via reverse proxy)
- [ ] Set up firewall rules
- [ ] Disabled root login (for SSH)
- [ ] Enabled fail2ban or similar
- [ ] Regular backups configured
- [ ] Monitoring enabled

## Support

For help with remote access configuration:
- GitHub Issues: https://github.com/Rowtion/AIpoch-claw/issues
- Documentation: See main README.md
