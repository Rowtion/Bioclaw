# 远程访问指南

启用 Bioclaw 服务（RStudio、JupyterLab、Opencode）的远程访问，从任何地方都能使用。

## 快速开始

### 1. 配置环境

```bash
cp .env.remote.template .env

# 编辑 .env 设置安全密码
nano .env
```

**重要：** 为了安全必须修改默认密码！

### 2. 启动远程访问服务

```bash
# 使用远程 docker-compose 文件
docker-compose -f docker-compose-remote.yml up -d

# 或者重命名为主文件：
mv docker-compose.yml docker-compose-local.yml
mv docker-compose-remote.yml docker-compose.yml
docker-compose up -d
```

### 3. 访问服务

启动后，从网络中的任何设备访问：

| 服务 | 本地地址 | 远程地址（替换 YOUR_IP） |
|---------|-----------|------------------------------|
| **RStudio Server** | http://localhost:8787 | http://YOUR_IP:8787 |
| **JupyterLab** | http://localhost:8888 | http://YOUR_IP:8888 |
| **Opencode** | http://localhost:4096 | http://YOUR_IP:4096 |

**查找你的 IP 地址：**
```bash
# macOS
ifconfig | grep "inet " | grep -v 127.0.0.1

# Linux
ip addr show | grep "inet " | grep -v 127.0.0.1
```

## 安全注意事项

⚠️ **警告：** 开放端口到互联网存在安全风险！

### 基础安全（最低要求）

1. **修改 .env 中的默认密码**：
```bash
RSTUDIO_PASSWORD=your_secure_password_123
JUPYTER_TOKEN=your_secure_token_456
```

2. **使用防火墙**限制访问：
```bash
# macOS - 只允许特定 IP
sudo ipfw add allow tcp from 192.168.1.0/24 to any 8787,8888,4096

# Linux (ufw)
sudo ufw allow from 192.168.1.0/24 to any port 8788,8888,4096
```

### 高级安全（互联网访问推荐）

#### 方案 A：SSH 隧道（最安全）

不要直接暴露端口，使用 SSH 隧道：

```bash
# 在你的本地机器上创建到远程服务器的隧道
ssh -L 8787:localhost:8787 -L 8888:localhost:8888 -L 4096:localhost:4096 user@your-server-ip

# 然后通过本机的 localhost 访问
# http://localhost:8787
# http://localhost:8888
# http://localhost:4096
```

#### 方案 B：Nginx 反向代理 + SSL

使用 Nginx 作为反向代理并提供 HTTPS：

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

#### 方案 C：VPN（团队最佳选择）

使用 Tailscale、ZeroTier 或 WireGuard 建立安全私有网络：

```bash
# 安装 Tailscale
curl -fsSL https://tailscale.com/install.sh | sh
sudo tailscale up

# 现在 Tailscale 网络中的所有设备都可以访问
# http://your-machine-tailscale-ip:8787
```

## 云端部署

### AWS/GCP/Azure

1. 开放安全组端口：8787、8888、4096
2. 使用强密码
3. 考虑使用 AWS Systems Manager Session Manager 代替开放端口

### Docker + Cloudflared（免费隧道）

```yaml
version: '3.8'
services:
  tunnel:
    image: cloudflare/cloudflared:latest
    command: tunnel --no-autoupdate run --token YOUR_TOKEN
    restart: unless-stopped
    
  analysis-env:
    # ... 同上
```

## 故障排除

### 无法远程访问

1. **检查防火墙：**
```bash
# macOS
sudo /usr/libexec/ApplicationFirewall/socketfilterfw --getglobalstate

# Linux
sudo iptables -L | grep 8787
```

2. **检查服务是否监听 0.0.0.0：**
```bash
netstat -an | grep 8787
# 应该显示 0.0.0.0:8787，而不是 127.0.0.1:8787
```

3. **从另一台机器测试：**
```bash
curl http://your-server-ip:8787
```

### 连接被拒绝

- 确保 Docker 容器在运行：`docker-compose ps`
- 检查日志：`docker-compose logs analysis-env`
- 验证端口映射：`docker port aipoch-analysis`

## 生产环境检查清单

在暴露到互联网之前：

- [ ] 已修改所有默认密码
- [ ] 已启用 HTTPS（通过反向代理）
- [ ] 已设置防火墙规则
- [ ] 已禁用 root 登录（对于 SSH）
- [ ] 已启用 fail2ban 或类似工具
- [ ] 已配置定期备份
- [ ] 已启用监控

## 支持

如需远程访问配置帮助：
- GitHub Issues：https://github.com/Rowtion/Bioclaw/issues
- 文档：参见主 README.md
