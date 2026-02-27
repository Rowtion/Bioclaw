#!/bin/bash
# Bioclaw 安装前检查

echo "🔍 Bioclaw 安装前检查"
echo "======================"
echo ""

ERRORS=0
WARNINGS=0

# 检查系统
echo "1. 检查操作系统..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "   ✅ macOS 支持"
    OS="macos"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "   ✅ Linux 支持"
    OS="linux"
else
    echo "   ❌ 暂不支持: $OSTYPE"
    ((ERRORS++))
fi

# 检查 Docker
echo ""
echo "2. 检查 Docker..."
if command -v docker >/dev/null 2>&1; then
    echo "   ✅ Docker 已安装"
    
    if docker info >/dev/null 2>&1; then
        echo "   ✅ Docker 正在运行"
    else
        echo "   ❌ Docker 未运行"
        echo "      请启动 Docker Desktop（macOS）或运行: sudo systemctl start docker（Linux）"
        ((ERRORS++))
    fi
else
    echo "   ❌ Docker 未安装"
    echo "      macOS: https://docs.docker.com/desktop/install/mac-install/"
    echo "      Ubuntu: sudo apt-get install docker.io"
    ((ERRORS++))
fi

# 检查 Git
echo ""
echo "3. 检查 Git..."
if command -v git >/dev/null 2>&1; then
    echo "   ✅ Git 已安装"
else
    echo "   ❌ Git 未安装"
    echo "      macOS: 安装 Xcode Command Line Tools: xcode-select --install"
    echo "      Ubuntu: sudo apt-get install git"
    ((ERRORS++))
fi

# 检查磁盘空间
echo ""
echo "4. 检查磁盘空间..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    AVAILABLE=$(df -g . | tail -1 | awk '{print $4}')
else
    AVAILABLE=$(df -BG . | tail -1 | awk '{print $4}' | sed 's/G//')
fi

if [ "$AVAILABLE" -ge 10 ]; then
    echo "   ✅ 磁盘空间充足 (${AVAILABLE}GB 可用，需要 10GB+)"
else
    echo "   ⚠️  磁盘空间不足 (${AVAILABLE}GB 可用，建议 10GB+)"
    ((WARNINGS++))
fi

# 检查内存
echo ""
echo "5. 检查内存..."
if [[ "$OSTYPE" == "darwin"* ]]; then
    MEM=$(sysctl -n hw.memsize | awk '{print $1/1024/1024/1024}')
else
    MEM=$(free -g | awk '/^Mem:/{print $2}')
fi

if [ "${MEM%.*}" -ge 4 ]; then
    echo "   ✅ 内存充足 (${MEM%.*}GB)"
else
    echo "   ⚠️  内存较低 (${MEM%.*}GB，建议 4GB+)"
    ((WARNINGS++))
fi

# 检查网络
echo ""
echo "6. 检查网络连接..."
if curl -s --max-time 10 https://github.com >/dev/null; then
    echo "   ✅ 可以访问 GitHub"
else
    echo "   ⚠️  无法访问 GitHub（可能需要代理）"
    ((WARNINGS++))
fi

# 检查端口
echo ""
echo "7. 检查端口占用..."
PORTS="8787 8888 4096"
PORT_OK=true
for port in $PORTS; do
    if lsof -Pi :"$port" -sTCP:LISTEN -t >/dev/null 2>&1; then
        echo "   ⚠️  端口 $port 已被占用"
        PORT_OK=false
        ((WARNINGS++))
    fi
done

if [ "$PORT_OK" = true ]; then
    echo "   ✅ 端口可用"
fi

# 总结
echo ""
echo "======================"
if [ $ERRORS -gt 0 ]; then
    echo "❌ 检查失败，有 $ERRORS 个错误需要修复"
    echo ""
    echo "请先修复上述错误，然后重新运行安装脚本。"
    exit 1
elif [ $WARNINGS -gt 0 ]; then
    echo "⚠️  检查通过，但有 $WARNINGS 个警告"
    echo ""
    read -p "是否继续安装？(y/N): " confirm
    if [ "$confirm" != "y" ] && [ "$confirm" != "Y" ]; then
        echo "已取消"
        exit 0
    fi
else
    echo "✅ 所有检查通过，可以安装！"
    echo ""
    read -p "是否立即安装？(Y/n): " confirm
    if [ "$confirm" = "n" ] || [ "$confirm" = "N" ]; then
        echo "已取消"
        exit 0
    fi
fi

# 运行安装
echo ""
echo "🚀 开始安装..."
bash install.sh
