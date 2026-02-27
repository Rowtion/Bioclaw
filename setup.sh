#!/bin/bash
# AIpoch-claw Setup Script
# 设置 OpenClaw + Opencode + AIPOCH Medical Research Skills

set -e

echo "🚀 AIpoch-claw 设置脚本"
echo "======================="

# 颜色
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# 检查命令
command_exists() {
    command -v "$1" &> /dev/null
}

# 检查依赖
echo ""
echo "📋 检查依赖..."

if ! command_exists docker; then
    echo -e "${RED}❌ Docker 未安装${NC}"
    echo "   请先安装 Docker: https://docs.docker.com/get-docker/"
    exit 1
fi

if ! command_exists docker-compose; then
    echo -e "${RED}❌ docker-compose 未安装${NC}"
    exit 1
fi

if ! command_exists git; then
    echo -e "${RED}❌ Git 未安装${NC}"
    exit 1
fi

echo -e "${GREEN}✅ 所有依赖已安装${NC}"

# 创建工作目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 克隆 AIPOCH Medical Research Skills
echo ""
echo "📚 下载 AIPOCH Medical Research Skills..."

if [ -d "scientific-skills" ]; then
    echo -e "${YELLOW}⚠️  scientific-skills 目录已存在${NC}"
    read -p "是否更新? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        cd scientific-skills
        git pull
        cd ..
    fi
else
    git clone --depth 1 https://github.com/aipoch/medical-research-skills.git scientific-skills
    echo -e "${GREEN}✅ AIPOCH 技能库已下载${NC}"
fi

# 统计技能数量
SKILL_COUNT=$(find scientific-skills -name "SKILL.md" 2>/dev/null | wc -l)
echo "   发现 $SKILL_COUNT 个技能"

# 创建数据目录
echo ""
echo "📁 创建数据目录..."
mkdir -p data outputs

# 创建 .env 文件
echo ""
echo "⚙️  配置文件..."

if [ -f ".env" ]; then
    echo -e "${YELLOW}⚠️  .env 文件已存在${NC}"
else
    cat > .env << EOF
# Opencode 配置
OPENCODE_URL=http://localhost:4096
OPENCODE_TIMEOUT=300

# 模型配置 (如果使用第三方 API)
# ANTHROPIC_API_KEY=your_key_here
# ANTHROPIC_SMALL_FAST_MODEL=claude-sonnet-4-20250514

# Docker 配置
JUPYTER_TOKEN=aipoch
RSTUDIO_PASSWORD=aipoch
EOF
    echo -e "${GREEN}✅ .env 文件已创建${NC}"
    echo "   请编辑 .env 文件添加你的 API 密钥"
fi

# 构建 Docker 镜像
echo ""
echo "🐳 构建 Docker 镜像..."
docker-compose build

# 启动服务
echo ""
echo "🚀 启动服务..."
docker-compose up -d

echo ""
echo -e "${GREEN}✅ 设置完成!${NC}"
echo ""
echo "📊 服务状态:"
echo "   • RStudio Server: http://localhost:8787 (密码: aipoch)"
echo "   • JupyterLab:     http://localhost:8888 (token: aipoch)"
echo ""
echo "🔧 下一步:"
echo "   1. 确保 Opencode 已安装并运行:"
echo "      ~/.opencode/bin/opencode serve --port 4096"
echo ""
echo "   2. 在 OpenClaw 中安装 bridge skill:"
echo "      技能已位于: ~/.openclaw/workspace/skills/opencode-bridge"
echo ""
echo "   3. 开始使用:"
echo "      在 Feishu 中发送: '用opencode分析数据'"
echo ""
echo "📖 查看日志:"
echo "   docker-compose logs -f"
