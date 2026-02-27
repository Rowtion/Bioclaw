# 常见问题 (FAQ)

## 安装问题

### Q: 安装到一半卡住不动了？
**A:** Docker 镜像构建需要 5-10 分钟，这是正常的。如果超过 15 分钟没动静：
1. 按 `Ctrl+C` 停止
2. 检查网络连接
3. 重新运行安装命令

### Q: 提示 "Docker 未运行" 但我已经打开了？
**A:** macOS 上 Docker Desktop 启动需要 30-60 秒，请等待状态栏显示 "Docker Desktop is running" 后再运行安装脚本。

### Q: 安装失败，如何重新开始？
**A:** 
```bash
# 清理残留
rm -rf ~/.bioclaw
rm -f ~/.local/bin/bioclaw

# 重新安装
curl -fsSL https://.../install.sh | bash
```

## 使用问题

### Q: 怎么打开 RStudio/JupyterLab？
**A:** 在浏览器中输入：
- RStudio: http://localhost:8787
- JupyterLab: http://localhost:8888

### Q: 忘记密码了？
**A:** 默认密码是 `bioclaw`。如需修改：
```bash
cd ~/.bioclaw
nano .env  # 修改 RSTUDIO_PASSWORD 和 JUPYTER_TOKEN
docker-compose restart
```

### Q: 怎么上传自己的数据？
**A:** 将文件放到 `~/.bioclaw/data/` 目录，在 RStudio/JupyterLab 中就能看到。

### Q: 分析结果保存在哪里？
**A:** 结果自动保存在 `~/.bioclaw/outputs/` 目录。

## 故障排除

### Q: 端口被占用了怎么办？
**A:** 找到占用端口的程序并关闭：
```bash
# macOS
lsof -i :8787
kill -9 <PID>

# 或修改端口
cd ~/.bioclaw
nano docker-compose.yml  # 修改 ports 部分
```

### Q: Docker 提示磁盘空间不足？
**A:** 清理 Docker 镜像：
```bash
docker system prune -a
```

### Q: 如何完全卸载？
**A:** 运行卸载脚本：
```bash
cd ~/.bioclaw
bash uninstall.sh
```

## 其他

### Q: 支持 Windows 吗？
**A:** 目前只支持 macOS 和 Linux。Windows 用户可以使用 WSL2。

### Q: 数据安全吗？
**A:** 所有数据都保存在你的本地机器上，不会上传到云端。

### Q: 怎么更新到最新版？
**A:** 
```bash
bioclaw update
```

### Q: 还有其他问题？
**A:** 提交 GitHub Issue: https://github.com/Rowtion/Bioclaw/issues
