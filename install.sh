#!/bin/bash
# Bioclaw - 小白一键安装脚本
# 自动检测、自动安装、自动配置

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

# 检查命令是否存在
command_exists() {
    command -v "$1" &> /dev/null
}

# 检查系统
print_status "检查系统环境..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
else
    print_error "暂不支持此系统: $OSTYPE"
    exit 1
fi
print_success "系统检查通过: $OS"

# 检查并安装 Docker
print_status "检查 Docker..."
if ! command_exists docker; then
    print_warning "Docker 未安装"
    echo ""
    echo "📥 正在自动安装 Docker..."
    
    if [ "$OS" = "macos" ]; then
        echo "请从 https://docs.docker.com/desktop/install/mac-install/ 下载 Docker Desktop"
        echo "安装完成后重新运行此脚本"
        exit 1
    else
        # Linux
        curl -fsSL https://get.docker.com | sh
        sudo usermod -aG docker $USER
        print_warning "请注销并重新登录，或运行: newgrp docker"
        exit 1
    fi
else
    # 检查 Docker 是否运行
    if ! docker info >/dev/null 2>&1; then
        print_warning "Docker 已安装但未运行"
        echo ""
        if [ "$OS" = "macos" ]; then
            echo "🚀 正在启动 Docker Desktop..."
            open -a Docker
        else
            echo "🚀 正在启动 Docker..."
            sudo systemctl start docker
        fi
        echo "⏳ 等待 Docker 启动 (约 30 秒)..."
        sleep 30
        
        # 再次检查
        if ! docker info >/dev/null 2>&1; then
            print_error "Docker 启动失败，请手动启动后重试"
            exit 1
        fi
    fi
    print_success "Docker 运行正常"
fi

# 检查 Git
print_status "检查 Git..."
if ! command_exists git; then
    print_warning "Git 未安装，正在安装..."
    if [ "$OS" = "macos" ]; then
        if command_exists brew; then
            brew install git
        else
            echo "请先安装 Homebrew: https://brew.sh"
            exit 1
        fi
    else
        sudo apt-get update && sudo apt-get install -y git
    fi
fi
print_success "Git 已安装"

# 检查 Node.js (用于 Opencode)
print_status "检查 Node.js..."
if ! command_exists node; then
    print_warning "Node.js 未安装，正在安装..."
    if [ "$OS" = "macos" ]; then
        if command_exists brew; then
            brew install node
        else
            curl -fsSL https://nodejs.org/dist/v20.11.0/node-v20.11.0.pkg -o /tmp/node.pkg
            sudo installer -pkg /tmp/node.pkg -target /
        fi
    else
        curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
        sudo apt-get install -y nodejs
    fi
fi
print_success "Node.js 已安装 ($(node --version))"

# 克隆项目
PROJECT_DIR="$HOME/.bioclaw"
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

# 创建启动脚本
print_status "创建快捷命令..."

mkdir -p "$HOME/.local/bin"

# bioclaw 命令
cat > "$HOME/.local/bin/bioclaw" <> 'EOF'
#!/bin/bash
case "$1" in
    start)
        echo "🚀 启动 Bioclaw..."
        cd "$HOME/.bioclaw" && docker-compose up -d
        echo ""
        echo "✅ Bioclaw 已启动!"
        echo ""
        echo "📊 访问地址:"
        echo "   RStudio:    http://localhost:8787"
        echo "   JupyterLab: http://localhost:8888"
        echo ""
        ;;
    stop)
        echo "🛑 停止 Bioclaw..."
        cd "$HOME/.bioclaw" && docker-compose down
        echo "✅ Bioclaw 已停止"
        ;;
    status)
        cd "$HOME/.bioclaw" && docker-compose ps
        ;;
    logs)
        cd "$HOME/.bioclaw" && docker-compose logs -f
        ;;
    update)
        echo "🔄 更新 Bioclaw..."
        cd "$HOME/.bioclaw" && git pull && docker-compose pull
        echo "✅ 更新完成，请运行: bioclaw restart"
        ;;
    restart)
        bioclaw stop && bioclaw start
        ;;
    *)
        echo "Bioclaw - 生物科研环境管理"
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
        echo "  RStudio:    http://localhost:8787 (密码: bioclaw)"
        echo "  JupyterLab: http://localhost:8888 (token: bioclaw)"
        ;;
esac
EOF

chmod +x "$HOME/.local/bin/bioclaw"

# 添加到 PATH
if ! echo "$PATH" | grep -q "$HOME/.local/bin"; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.zshrc" 2>/dev/null || true
    export PATH="$HOME/.local/bin:$PATH"
fi

# 首次构建和启动
print_status "首次构建 (约需 5-10 分钟)..."
cd "$PROJECT_DIR"
docker-compose build --no-cache

print_status "启动服务..."
docker-compose up -d

# 等待服务就绪
print_status "等待服务就绪..."
sleep 10

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
echo "   bioclaw logs    - 查看日志"
echo ""
echo "📁 项目位置: $PROJECT_DIR"
echo ""
echo -e "${YELLOW}⚠️  重要: 请修改默认密码以确保安全!${NC}"
echo "   编辑: $PROJECT_DIR/.env"
echo ""
