#!/bin/bash
# Bioclaw 安装脚本
# 前置条件：Docker 已安装并运行

set -e

# 颜色
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_status() { echo -e "${BLUE}⏳ $1${NC}"; }
print_success() { echo -e "${GREEN}✅ $1${NC}"; }
print_error() { echo -e "${RED}❌ $1${NC}"; }
print_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; }

# 检查 Docker
print_status "检查 Docker..."
if ! command -v docker >/dev/null 2>&1; then
    print_error "Docker 未安装"
    echo ""
    echo "📥 请先安装 Docker："
    echo "   macOS: https://docs.docker.com/desktop/install/mac-install/"
    echo "   Ubuntu: sudo apt-get install docker.io"
    echo ""
    echo "安装完成后重新运行此脚本"
    exit 1
fi

# 检查 Docker 是否运行
if ! docker info >/dev/null 2>&1; then
    print_error "Docker 未运行"
    echo ""
    echo "🚀 请启动 Docker："
    echo "   macOS: 打开 Docker Desktop 应用"
    echo "   Ubuntu: sudo systemctl start docker"
    exit 1
fi
print_success "Docker 运行正常"

# 检查 Git
print_status "检查 Git..."
if ! command -v git >/dev/null 2>&1; then
    print_error "Git 未安装"
    echo "请安装 Git 后重试"
    exit 1
fi
print_success "Git 已安装"

# 设置安装目录
PROJECT_DIR="$HOME/.bioclaw"

# 下载 Bioclaw
if [ ! -d "$PROJECT_DIR" ]; then
    print_status "下载 Bioclaw..."
    git clone --depth 1 https://github.com/Rowtion/Bioclaw.git "$PROJECT_DIR"
    print_success "下载完成"
else
    print_status "更新 Bioclaw..."
    cd "$PROJECT_DIR" && git pull
    print_success "更新完成"
fi

cd "$PROJECT_DIR"

# 创建快捷命令
print_status "创建 bioclaw 命令..."
mkdir -p "$HOME/.local/bin"

cat > "$HOME/.local/bin/bioclaw" <> 'EOF'
#!/bin/bash
cd "$HOME/.bioclaw"
case "$1" in
    start)
        echo "🚀 启动 Bioclaw..."
        docker-compose up -d
        echo ""
        echo "✅ Bioclaw 已启动!"
        echo ""
        echo "📊 访问地址:"
        echo "   RStudio:    http://localhost:8787"
        echo "   JupyterLab: http://localhost:8888"
        echo ""
        echo "🔑 默认密码: bioclaw"
        ;;
    stop)
        echo "🛑 停止 Bioclaw..."
        docker-compose down
        echo "✅ Bioclaw 已停止"
        ;;
    restart)
        echo "🔄 重启 Bioclaw..."
        docker-compose restart
        echo "✅ Bioclaw 已重启"
        ;;
    status)
        docker-compose ps
        ;;
    logs)
        docker-compose logs -f
        ;;
    update)
        echo "🔄 更新 Bioclaw..."
        git pull && docker-compose pull
        echo "✅ 更新完成"
        ;;
    *)
        echo "Bioclaw 管理工具"
        echo ""
        echo "用法: bioclaw [命令]"
        echo ""
        echo "命令:"
        echo "  start    启动 Bioclaw"
        echo "  stop     停止 Bioclaw"
        echo "  restart  重启 Bioclaw"
        echo "  status   查看运行状态"
        echo "  logs     查看日志"
        echo "  update   更新 Bioclaw"
        echo ""
        echo "访问地址:"
        echo "  RStudio:    http://localhost:8787"
        echo "  JupyterLab: http://localhost:8888"
        echo ""
        echo "默认密码: bioclaw"
        ;;
esac
EOF

chmod +x "$HOME/.local/bin/bioclaw"

# 添加到 PATH
if ! echo "$PATH" | grep -q "$HOME/.local/bin"; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
    export PATH="$HOME/.local/bin:$PATH"
fi

# 首次构建
print_status "首次构建 (约需 5-10 分钟)..."
docker-compose build

# 启动
print_status "启动 Bioclaw..."
docker-compose up -d

# 显示完成信息
clear
echo ""
echo -e "${GREEN}╔════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║  🎉 Bioclaw 安装完成!                                ║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "📊 访问地址:"
echo "   RStudio Server: http://localhost:8787"
echo "   JupyterLab:     http://localhost:8888"
echo ""
echo "🔑 默认密码: bioclaw"
echo ""
echo "💡 常用命令:"
echo "   bioclaw start   - 启动"
echo "   bioclaw stop    - 停止"
echo "   bioclaw status  - 查看状态"
echo ""
echo "📁 项目位置: $PROJECT_DIR"
echo ""
