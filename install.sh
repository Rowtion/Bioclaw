#!/bin/bash
# Bioclaw 带进度显示的安装脚本

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
print_progress() { echo -e "${YELLOW}📊 $1${NC}"; }

# 进度条函数
show_progress() {
    local duration=$1
    local message=$2
    local width=50
    local progress=0
    
    echo ""
    echo "$message"
    echo ""
    
    while [ $progress -le $width ]; do
        local percent=$((progress * 2))
        local filled=$(printf '%*s' "$progress" '' | tr ' ' '█')
        local empty=$(printf '%*s' $((width - progress)) '' | tr ' ' '░')
        printf "\r  [%s%s] %d%%" "$filled" "$empty" "$percent"
        sleep $((duration / width))
        ((progress++))
    done
    echo ""
}

# 第 1 步：检查 Docker
print_status "第 1/5 步：检查 Docker..."
if ! command -v docker >/dev/null 2>&1; then
    print_error "Docker 未安装"
    echo ""
    echo "📥 请先安装 Docker："
    echo "   macOS: https://docs.docker.com/desktop/install/mac-install/"
    echo "   Ubuntu: sudo apt-get install docker.io"
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    print_error "Docker 未运行"
    echo "🚀 请启动 Docker Desktop"
    exit 1
fi
print_success "Docker 检查通过"

# 第 2 步：下载代码
PROJECT_DIR="$HOME/.bioclaw"
print_status "第 2/5 步：下载 Bioclaw 代码..."

if [ ! -d "$PROJECT_DIR" ]; then
    git clone --depth 1 https://github.com/Rowtion/Bioclaw.git "$PROJECT_DIR" >/dev/null 2>&1
else
    cd "$PROJECT_DIR" && git pull >/dev/null 2>&1
fi
print_success "代码下载完成"

# 第 3 步：创建命令
print_status "第 3/5 步：创建快捷命令..."
mkdir -p "$HOME/.local/bin"

cat > "$HOME/.local/bin/bioclaw" <> 'EOF'
#!/bin/bash
cd "$HOME/.bioclaw"
case "$1" in
    start)
        echo "🚀 启动 Bioclaw..."
        docker-compose up -d
        echo ""
        echo "✅ 已启动!"
        echo ""
        echo "📊 访问地址:"
        echo "   RStudio:    http://localhost:8787"
        echo "   JupyterLab: http://localhost:8888"
        echo "   密码: bioclaw"
        ;;
    stop)
        echo "🛑 停止 Bioclaw..."
        docker-compose down
        echo "✅ 已停止"
        ;;
    status)
        docker-compose ps
        ;;
    *)
        echo "用法: bioclaw [start|stop|status]"
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

print_success "快捷命令创建完成"

# 第 4 步：构建镜像（带进度）
print_status "第 4/5 步：构建 Docker 镜像..."
print_progress "预计需要 5-10 分钟，请耐心等待..."

cd "$PROJECT_DIR"
docker-compose build --no-cache > /tmp/bioclaw-build.log 2>&1 &
BUILD_PID=$!

# 显示进度动画
spin='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'
i=0
while kill -0 $BUILD_PID 2>/dev/null; do
    i=$(( (i+1) % 10 ))
    printf "\r  ${spin:$i:1} 正在构建... (%d 秒)" $((SECONDS))
    sleep 1
done

# 检查结果
if wait $BUILD_PID; then
    echo ""
    print_success "镜像构建完成"
else
    echo ""
    print_error "构建失败"
    echo "查看日志: cat /tmp/bioclaw-build.log"
    exit 1
fi

# 第 5 步：启动服务
print_status "第 5/5 步：启动服务..."
docker-compose up -d

# 等待服务就绪
print_status "等待服务启动..."
sleep 5

# 完成
clear
echo ""
echo -e "${GREEN}╔════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║  🎉 安装完成！Bioclaw 已就绪                          ║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "📊 访问地址:"
echo "   RStudio:    http://localhost:8787"
echo "   JupyterLab: http://localhost:8888"
echo ""
echo "🔑 密码: bioclaw"
echo ""
echo "💡 常用命令:"
echo "   bioclaw start   启动"
echo "   bioclaw stop    停止"
echo "   bioclaw status  查看状态"
echo ""
