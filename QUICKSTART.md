# Bioclaw 小白安装指南

## 是什么？

Bioclaw 是一个**开箱即用的安装包**，帮你一键搭建生物科研 AI 环境。

不用懂 Docker、不用配环境，一条命令搞定所有！

## 一键安装

打开终端，粘贴运行：

```bash
curl -fsSL https://raw.githubusercontent.com/Rowtion/Bioclaw/main/install.sh | bash
```

**安装过程全自动：**
- ✅ 自动检测并安装 Docker
- ✅ 自动下载 Bioclaw
- ✅ 自动构建环境（约 5-10 分钟）
- ✅ 自动启动所有服务

## 装完后怎么用？

### 启动 Bioclaw
```bash
bioclaw start
```

### 停止 Bioclaw
```bash
bioclaw stop
```

### 查看状态
```bash
bioclaw status
```

### 查看日志
```bash
bioclaw logs
```

## 访问你的科研环境

安装完成后，在浏览器中打开：

- **RStudio**（数据分析）：http://localhost:8787
- **JupyterLab**（Python 编程）：http://localhost:8888

**默认密码：** `bioclaw`

## 下一步

1. **打开 RStudio** → 开始分析你的数据
2. **配置 AI 模型** → 运行 `~/.opencode/bin/opencode auth login`
3. **在飞书中对话** → 发送 "帮我分析这个数据"

## 常见问题

### Q: 安装失败了怎么办？
A: 确保网络畅通，然后重新运行安装命令

### Q: Docker 是什么？
A: 不用管，Bioclaw 会自动帮你安装和配置

### Q: 如何更新 Bioclaw？
A: 运行 `bioclaw update`

### Q: 数据会丢失吗？
A: 不会，你的数据保存在 `~/.bioclaw/data` 和 `~/.bioclaw/outputs` 中

## 需要帮助？

- GitHub: https://github.com/Rowtion/Bioclaw
- 提交问题: https://github.com/Rowtion/Bioclaw/issues
