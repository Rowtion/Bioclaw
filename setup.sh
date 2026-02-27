#!/bin/bash
# Bioclaw Setup Script v2.1
# One-click installation: OpenClaw, Opencode, RStudio, Jupyter, and K-Dense skills

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================
readonly SCRIPT_VERSION="2.1"
readonly REQUIRED_PORTS=(8787 8888 4096)
readonly MIN_DISK_GB=10
readonly LOG_FILE="setup-$(date +%Y%m%d-%H%M%S).log"

# Installation paths
readonly OPENCODE_INSTALL_DIR="$HOME/.opencode"
readonly OPENCLAW_CONFIG_DIR="$HOME/.openclaw"

# Colors
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly RED='\033[0;31m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m'

# State tracking for rollback
declare -a CREATED_CONTAINERS=()
declare -a CREATED_VOLUMES=()
declare -a INSTALLED_TOOLS=()

# =============================================================================
# Utility Functions
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

print_status() {
    echo -e "${BLUE}ℹ️  $1${NC}"
    log "INFO: $1"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
    log "SUCCESS: $1"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
    log "WARNING: $1"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
    log "ERROR: $1"
}

print_header() {
    echo -e "${CYAN}"
    echo "═══════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "═══════════════════════════════════════════════════════════════"
    echo -e "${NC}"
}

command_exists() {
    command -v "$1" &> /dev/null
}

# =============================================================================
# Cleanup and Rollback
# =============================================================================

cleanup_on_error() {
    local exit_code=$?
    print_error "Error occurred during installation (exit code: $exit_code), cleaning up..."
    log "Starting rollback procedure"
    
    # Stop and remove containers created by this script
    if [ ${#CREATED_CONTAINERS[@]} -gt 0 ]; then
        print_status "Cleaning up containers..."
        for container in "${CREATED_CONTAINERS[@]}"; do
            docker stop "$container" 2>/dev/null || true
            docker rm "$container" 2>/dev/null || true
        done
    fi
    
    log "Rollback completed"
    print_error "Installation failed. Check logs: $LOG_FILE"
    print_status "Installed components may need manual cleanup: ${INSTALLED_TOOLS[*]}"
    exit 1
}

# Set trap for cleanup
trap cleanup_on_error ERR

# =============================================================================
# Pre-flight Checks
# =============================================================================

check_prerequisites() {
    print_header "System Check"
    
    # Check OS
    if [[ "$OSTYPE" != "linux-gnu"* ]] && [[ "$OSTYPE" != "darwin"* ]]; then
        print_warning "Not tested on Linux/macOS, may not be compatible: $OSTYPE"
    fi
    
    # Check if running as root
    if [ "$EUID" -eq 0 ]; then
        print_warning "Running as root may have security risks, recommend using regular user"
    fi
    
    print_success "System check passed"
}

check_dependencies() {
    print_header "Checking Dependencies"
    local missing_deps=()
    
    # Check Docker
    if ! command_exists docker; then
        missing_deps+=("docker")
        print_error "Docker not installed"
        echo ""
        echo "   🔧 Installation methods:"
        echo "   • macOS:    brew install --cask docker"
        echo "   • Ubuntu:   sudo apt-get install docker.io"
        echo "   • Official guide: https://docs.docker.com/get-docker/"
        echo ""
    else
        # Check Docker daemon is running
        if ! docker info &>/dev/null; then
            print_error "Docker daemon not running"
            echo ""
            echo "   🔧 How to start:"
            echo "   • macOS:    open -a Docker"
            echo "   • Linux:    sudo systemctl start docker"
            echo ""
            exit 1
        fi
        local docker_version
        docker_version=$(docker --version | awk '{print $3}' | sed 's/,//')
        print_success "Docker $docker_version installed and running normally"
    fi
    
    # Check docker-compose or docker compose
    if command_exists docker-compose; then
        DOCKER_COMPOSE="docker-compose"
        print_success "docker-compose installed"
    elif docker compose version &>/dev/null; then
        DOCKER_COMPOSE="docker compose"
        print_success "docker compose (plugin) installed"
    else
        missing_deps+=("docker-compose")
        print_error "docker-compose not installed"
        echo "   Installation guide: https://docs.docker.com/compose/install/"
    fi
    
    # Check Git
    if ! command_exists git; then
        missing_deps+=("git")
        print_error "Git not installed"
    else
        local git_version
        git_version=$(git --version | awk '{print $3}')
        print_success "Git $git_version installed"
    fi
    
    # Check curl
    if ! command_exists curl; then
        missing_deps+=("curl")
        print_error "curl not installed"
    else
        print_success "curl installed"
    fi
    
    if [ ${#missing_deps[@]} -gt 0 ]; then
        print_error "Missing dependencies: ${missing_deps[*]}"
        echo ""
        echo "Please install the above dependencies before running this script."
        exit 1
    fi
    
    print_success "All base dependencies installed"
}

check_ports() {
    print_header "Checking Ports"
    local port_in_use=()
    
    for port in "${REQUIRED_PORTS[@]}"; do
        if lsof -Pi :"$port" -sTCP:LISTEN -t >/dev/null 2>&1 || \
           netstat -tuln 2>/dev/null | grep -q ":$port " || \
           ss -tuln 2>/dev/null | grep -q ":$port "; then
            port_in_use+=("$port")
            print_error "Port $port is already in use"
            
            # Try to identify what's using the port
            local process_info
            process_info=$(lsof -ti :"$port" 2>/dev/null | head -1 || echo "unknown")
            if [ "$process_info" != "unknown" ]; then
                echo "   Process using PID: $process_info"
                ps -p "$process_info" -o comm= 2>/dev/null | sed 's/^/   Process name: /' || true
            fi
        fi
    done
    
    if [ ${#port_in_use[@]} -gt 0 ]; then
        echo ""
        print_error "Following ports are in use: ${port_in_use[*]}"
        echo ""
        echo "   🔧 Solutions:"
        echo "   1. Stop programs using these ports"
        echo "   2. Or modify docker-compose.yml to use different ports"
        echo ""
        echo "   Port usage:"
        echo "   - 4096: Opencode Server"
        echo "   - 8787: RStudio Server"
        echo "   - 8888: JupyterLab"
        exit 1
    fi
    
    print_success "All ports available: ${REQUIRED_PORTS[*]}"
}

check_disk_space() {
    print_header "Checking Disk Space"
    
    # Get available space in GB
    local available_gb
    if [[ "$OSTYPE" == "darwin"* ]]; then
        available_gb=$(df -g . | tail -1 | awk '{print $4}')
    else
        available_gb=$(df -BG . | tail -1 | awk '{print $4}' | sed 's/G//')
    fi
    
    print_status "Available disk space: ${available_gb}GB (requires at least ${MIN_DISK_GB}GB)"
    
    if [ "$available_gb" -lt "$MIN_DISK_GB" ]; then
        print_error "Insufficient disk space: ${available_gb}GB available, need at least ${MIN_DISK_GB}GB"
        echo ""
        echo "   Docker images require significant space, please clean up disk and retry."
        exit 1
    fi
    
    print_success "Disk space sufficient"
}

check_in_project_root() {
    print_status "Validating project location..."
    
    if [ ! -f "docker-compose.yml" ]; then
        print_error "Please run this script from the Bioclaw project root directory"
        echo "   docker-compose.yml not found"
        exit 1
    fi
    
    if [ ! -d "docker" ]; then
        print_error "Project structure incomplete, missing docker/ directory"
        exit 1
    fi
    
    if [ ! -f ".env.template" ]; then
        print_warning ".env.template not found"
    fi
    
    print_success "Project structure validation passed"
}

# =============================================================================
# Tool Installation
# =============================================================================

check_nodejs() {
    if command_exists node; then
        local node_version
        node_version=$(node --version)
        print_success "Node.js $node_version installed"
        return 0
    else
        print_error "Node.js not installed (required for OpenClaw installation)"
        return 1
    fi
}

install_nodejs() {
    print_status "Installing Node.js..."
    
    if [[ "$OSTYPE" == "darwin"* ]]; then
        if command_exists brew; then
            brew install node
        else
            print_error "Homebrew not found, please install first: https://brew.sh"
            exit 1
        fi
    else
        # Linux - use NodeSource
        curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
        sudo apt-get install -y nodejs
    fi
    
    if command_exists node; then
        print_success "Node.js $(node --version) installation completed"
    else
        print_error "Node.js installation failed"
        exit 1
    fi
}

install_openclaw() {
    print_header "Installing OpenClaw"
    
    # Check if already installed
    if command_exists openclaw; then
        local version
        version=$(openclaw --version 2>&1 | head -1 || echo "unknown")
        print_success "OpenClaw already installed: $version"
        return 0
    fi
    
    # Check if installed but not in PATH
    if [ -f "$OPENCLAW_CONFIG_DIR/bin/openclaw" ] || [ -f "/usr/local/bin/openclaw" ]; then
        print_warning "OpenClaw installed but not in PATH"
        echo ""
        echo "   Please add the following to your ~/.bashrc or ~/.zshrc:"
        echo "   export PATH=\"$OPENCLAW_CONFIG_DIR/bin:\$PATH\""
        echo ""
        echo "   Then run: source ~/.bashrc (or ~/.zshrc)"
        echo ""
        
        read -p "Continue installing other components? (Y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Nn]$ ]]; then
            exit 1
        fi
        return 0
    fi
    
    print_status "Preparing to install OpenClaw..."
    
    # Check Node.js
    if ! check_nodejs; then
        echo ""
        read -p "Install Node.js automatically? (Y/n): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            install_nodejs
        else
            print_error "Node.js is required to install OpenClaw"
            echo "   Please install Node.js manually and retry: https://nodejs.org/"
            exit 1
        fi
    fi
    
    echo ""
    read -p "Install OpenClaw? (Y/n): " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        print_warning "Skipping OpenClaw installation"
        return 0
    fi
    
    print_status "Installing OpenClaw..."
    
    if npm install -g openclaw 2>&1 | tee -a "$LOG_FILE"; then
        INSTALLED_TOOLS+=("openclaw")
        print_success "OpenClaw installation completed"
        
        # Initialize OpenClaw
        print_status "Initializing OpenClaw..."
        if openclaw init 2>&1 | tee -a "$LOG_FILE"; then
            print_success "OpenClaw initialization completed"
        else
            print_warning "OpenClaw initialization may need to be run manually: openclaw init"
        fi
        
        # Show PATH info
        echo ""
        echo "   💡 Tip: If openclaw command not found, run:"
        echo "      export PATH=\"\$(npm bin -g):\$PATH\""
        echo ""
    else
        print_error "OpenClaw installation failed"
        echo ""
        echo "   Possible causes:"
        echo "   - npm permission issues (try: npm config set prefix ~/.npm-global)"
        echo "   - Network connectivity issues"
        echo "   - Node.js version incompatibility"
        echo ""
        echo "   Manual install: npm install -g openclaw"
        exit 1
    fi
}

install_opencode() {
    print_header "Installing Opencode"
    
    # Check if already installed
    if command_exists opencode; then
        local version
        version=$(opencode --version 2>&1 | head -1 || echo "unknown")
        print_success "Opencode already installed: $version"
        
        # Check if opencode server is already running
        if curl -s http://localhost:4096/status >/dev/null 2>&1; then
            print_success "Opencode Server already running (http://localhost:4096)"
        fi
        return 0
    fi
    
    # Check if installed but not in PATH
    if [ -f "$OPENCODE_INSTALL_DIR/bin/opencode" ]; then
        print_warning "Opencode installed but not in PATH"
        echo ""
        echo "   Please add the following to your ~/.bashrc or ~/.zshrc:"
        echo "   export PATH=\"$OPENCODE_INSTALL_DIR/bin:\$PATH\""
        echo ""
        echo "   Then run: source ~/.bashrc (or ~/.zshrc)"
        echo ""
        
        read -p "Continue installing other components? (Y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Nn]$ ]]; then
            exit 1
        fi
        return 0
    fi
    
    print_status "Preparing to install Opencode..."
    echo ""
    read -p "Install Opencode? (Y/n): " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Nn]$ ]]; then
        print_warning "Skipping Opencode installation"
        return 0
    fi
    
    print_status "Installing Opencode..."
    print_warning "This may take a few minutes..."
    
    if curl -fsSL https://opencode.dev/install.sh | sh 2>&1 | tee -a "$LOG_FILE"; then
        INSTALLED_TOOLS+=("opencode")
        print_success "Opencode installation completed"
        
        # Add to PATH for current session
        export PATH="$OPENCODE_INSTALL_DIR/bin:$PATH"
        
        # Verify installation
        if [ -f "$OPENCODE_INSTALL_DIR/bin/opencode" ]; then
            print_success "Opencode installation verification passed"
            
            echo ""
            echo "   💡 Tip: If opencode command not found, add the following to ~/.bashrc or ~/.zshrc:"
            echo "      export PATH=\"$OPENCODE_INSTALL_DIR/bin:\$PATH\""
            echo ""
        else
            print_warning "Opencode installation path not found, may need to manually add to PATH"
        fi
    else
        print_error "Opencode installation failed"
        echo ""
        echo "   Possible causes:"
        echo "   - Network connectivity issues"
        echo "   - System incompatibility"
        echo ""
        echo "   Manual install: curl -fsSL https://opencode.dev/install.sh | sh"
        echo "   Or visit: https://opencode.dev"
        exit 1
    fi
}

start_opencode_server() {
    print_header "Starting Opencode Server"
    
    # Check if already running
    if curl -s http://localhost:4096/status >/dev/null 2>&1; then
        print_success "Opencode Server already running (http://localhost:4096)"
        return 0
    fi
    
    # Check if opencode command exists
    if ! command_exists opencode; then
        print_warning "opencode command not found, skipping server startup"
        return 0
    fi
    
    print_status "Starting Opencode Server..."
    
    # Start in background
    nohup opencode serve --port 4096 > logs/opencode.log 2>&1 &
    local opencode_pid=$!
    
    # Wait for server to be ready
    print_status "Waiting for Opencode Server to be ready..."
    local retries=30
    local wait_time=1
    
    for i in $(seq 1 $retries); do
        echo -n "."
        if curl -s http://localhost:4096/status >/dev/null 2>&1; then
            echo
            print_success "Opencode Server ready (PID: $opencode_pid)"
            
            echo ""
            echo "   💡 Opencode Server is running in the background"
            echo "   View logs: tail -f logs/opencode.log"
            echo "   Stop service: kill $opencode_pid"
            echo ""
            return 0
        fi
        
        # Check if process is still running
        if ! kill -0 $opencode_pid 2>/dev/null; then
            echo
            print_error "Opencode Server failed to start (process exited)"
            echo "   View logs: logs/opencode.log"
            return 1
        fi
        
        sleep $wait_time
    done
    
    echo
    print_warning "Opencode Server startup timed out, but process is running"
    echo "   Check manually later: curl http://localhost:4096/status"
    return 0
}

# =============================================================================
# Setup Steps
# =============================================================================

init_git_submodule() {
    print_header "Initializing K-Dense Scientific Skills"
    
    # Check if .gitmodules exists
    if [ ! -f ".gitmodules" ]; then
        print_warning ".gitmodules file not found, skipping submodule initialization"
        return 0
    fi
    
    # Check network connectivity
    if ! curl -s --max-time 10 https://github.com >/dev/null 2>&1; then
        print_warning "Cannot connect to GitHub, submodule initialization may fail"
        echo "   Please check network connection, or use offline mode (--offline)"
    fi
    
    if [ -d "scientific-skills/.git" ]; then
        print_status "Submodule exists, updating..."
        if git submodule update --remote --merge 2>&1 | tee -a "$LOG_FILE"; then
            print_success "Submodule update completed"
        else
            print_warning "Submodule update failed, using local version"
        fi
    else
        print_status "Initializing submodule..."
        if git submodule init 2>&1 | tee -a "$LOG_FILE" && \
           git submodule update 2>&1 | tee -a "$LOG_FILE"; then
            print_success "Submodule initialization completed"
        else
            print_error "Submodule initialization failed"
            echo ""
            echo "   Possible causes:"
            echo "   - Network connectivity issues"
            echo "   - Incorrect submodule URL"
            echo "   - Repository access permissions"
            echo ""
            echo "   You can:"
            echo "   1. Check network and retry"
            echo "   2. Manual clone: git clone https://github.com/aipoch/medical-research-skills scientific-skills"
            exit 1
        fi
    fi
    
    # Count skills
    local skill_count
    skill_count=$(find scientific-skills -name "SKILL.md" 2>/dev/null | wc -l | tr -d ' ')
    if [ "$skill_count" -gt 0 ]; then
        print_success "Found $skill_count medical research skills"
    else
        print_warning "No skill files found, functionality may be limited"
    fi
}

create_directories() {
    print_status "Creating data directories..."
    
    mkdir -p data outputs logs
    
    # Set permissions (if not Windows)
    if [[ "$OSTYPE" != "msys" ]] && [[ "$OSTYPE" != "cygwin" ]]; then
        chmod 755 data outputs logs 2>/dev/null || true
    fi
    
    print_success "Directories created"
}

create_env_file() {
    print_header "Configuring Environment Variables"
    
    if [ -f ".env" ]; then
        print_warning ".env file already exists"
        read -p "Overwrite? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            print_status "Keeping existing .env file"
            return 0
        fi
        # Backup old .env
        cp .env ".env.backup.$(date +%Y%m%d-%H%M%S)"
        print_status "Original .env file backed up"
    fi
    
    # Generate random tokens for security
    local jupyter_token
    jupyter_token=$(openssl rand -hex 16 2>/dev/null || echo "bioclaw-$(date +%s)")
    
    cat > .env << EOF
# Bioclaw Environment Configuration
# Generated on $(date)

# =============================================================================
# Opencode Configuration
# =============================================================================
OPENCODE_URL=http://localhost:4096
OPENCODE_TIMEOUT=300

# NOTE: API keys are configured via 'opencode auth login', not in this file.
# Run: ~/.opencode/bin/opencode auth login
# Or: opencode auth login (if opencode is in your PATH)
# =============================================================================
# Docker Services Configuration
# =============================================================================
JUPYTER_TOKEN=$jupyter_token
RSTUDIO_PASSWORD=bioclaw

# =============================================================================
# Security Settings (Change in production!)
# =============================================================================
# RSTUDIO_PASSWORD=$(openssl rand -base64 32 2>/dev/null || echo "change-me-$(date +%s)")
EOF
    
    print_success ".env file created"
    print_warning "⚠️  Important: Configure API access with 'opencode auth login'"
}

build_docker_images() {
    print_header "Building Docker Images"
    print_warning "This may take 10-30 minutes depending on network speed and hardware"
    
    # Check if images already exist
    if docker images | grep -q "bioclaw"; then
        read -p "Existing images detected, rebuild? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            print_status "Skipping image build, using existing images"
            return 0
        fi
    fi
    
    # Build with progress output
    if $DOCKER_COMPOSE build --progress=plain 2>&1 | tee -a "$LOG_FILE"; then
        print_success "Docker images build completed"
    else
        print_error "Docker images build failed"
        echo ""
        echo "   Common causes:"
        echo "   - Network issues (cannot download base images)"
        echo "   - Insufficient disk space"
        echo "   - Docker daemon issues"
        echo ""
        echo "   View detailed logs: $LOG_FILE"
        exit 1
    fi
}

start_docker_services() {
    print_header "Starting Docker Services"
    
    # Stop existing containers if any
    $DOCKER_COMPOSE down 2>/dev/null || true
    
    # Start services
    if $DOCKER_COMPOSE up -d 2>&1 | tee -a "$LOG_FILE"; then
        print_success "Service startup command sent"
        CREATED_CONTAINERS+=("bioclaw-analysis")
    else
        print_error "Service startup failed"
        exit 1
    fi
    
    # Wait for services to be healthy
    print_status "Waiting for RStudio Server to be ready..."
    local retries=30
    local wait_time=2
    
    for i in $(seq 1 $retries); do
        echo -n "."
        
        # Check RStudio
        if curl -s http://localhost:8787 >/dev/null 2>&1; then
            echo
            print_success "RStudio Server ready"
            break
        fi
        
        # Check if container is still running
        if ! docker ps | grep -q "bioclaw-analysis"; then
            echo
            print_error "Container stopped unexpectedly, view logs:"
            echo "   $DOCKER_COMPOSE logs"
            exit 1
        fi
        
        sleep $wait_time
    done
    
    if [ $i -eq $retries ]; then
        echo
        print_warning "Service startup timed out, but container is running"
        echo "   Check manually later: http://localhost:8787"
    fi
}

run_aipoch_init() {
    print_header "Running Bioclaw Initialization"
    
    if [ ! -f "aipoch-init.sh" ]; then
        print_error "aipoch-init.sh script not found"
        return 1
    fi
    
    print_status "Injecting Bioclaw configuration into OpenClaw..."
    
    if bash aipoch-init.sh 2>&1 | tee -a "$LOG_FILE"; then
        print_success "Bioclaw configuration injection completed"
        
        echo ""
        echo "   💡 OpenClaw configured for Bioclaw mode"
        echo "   Takes effect after restarting gateway: openclaw gateway restart"
        echo ""
    else
        print_warning "Bioclaw initialization may not be complete"
        echo "   Run manually later: bash aipoch-init.sh"
    fi
}

health_check() {
    print_header "Health Check"
    local all_healthy=true
    
    # Check RStudio
    if curl -s http://localhost:8787 >/dev/null 2>&1; then
        print_success "RStudio Server: http://localhost:8787 ✓"
    else
        print_error "RStudio Server not responding"
        all_healthy=false
    fi
    
    # Check Jupyter
    if curl -s http://localhost:8888 >/dev/null 2>&1; then
        print_success "JupyterLab: http://localhost:8888 ✓"
    else
        print_warning "JupyterLab not responding (may need more time to start)"
    fi
    
    # Check Opencode
    if curl -s http://localhost:4096/status >/dev/null 2>&1; then
        print_success "Opencode Server: http://localhost:4096 ✓"
    else
        print_warning "Opencode Server not responding"
        all_healthy=false
    fi
    
    # Check Docker containers
    if docker ps | grep -q "bioclaw-analysis"; then
        print_success "Docker containers running normally"
    else
        print_error "Docker containers not running"
        all_healthy=false
    fi
    
    if [ "$all_healthy" = true ]; then
        print_success "All core services running normally"
        return 0
    else
        return 1
    fi
}

# =============================================================================
# Main
# =============================================================================

show_banner() {
    echo -e "${CYAN}"
    cat << 'EOF'
█████╗ ██╗██████╗  ██████╗  ██████╗██╗  ██╗
██╔══██╗██║██╔══██╗██╔═══██╗██╔════╝██║  ██║
███████║██║██████╔╝██║   ██║██║     ███████║
██╔══██║██║██╔═══╝ ██║   ██║██║     ██╔══██║
██║  ██║██║██║     ╚██████╔╝╚██████╗██║  ██║
╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝  ╚═════╝╚═╝  ╚═╝
EOF
    echo -e "${NC}"
    echo -e "${CYAN}Bioclaw One-Click Installation Script v${SCRIPT_VERSION}${NC}"
    echo -e "${BLUE}Log file: $LOG_FILE${NC}"
    echo ""
    echo "This script will install:"
    echo "  • OpenClaw (AI conversation gateway)"
    echo "  • Opencode (code execution environment)"
    echo "  • RStudio Server + JupyterLab (Docker)"
    echo "  • K-Dense scientific computing skill library (140+)"
    echo ""
}

show_summary() {
    echo ""
    echo -e "${GREEN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║                    ✅ Installation Complete!                 ║${NC}"
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    
    # Read token from .env
    local jupyter_token
    jupyter_token=$(grep "JUPYTER_TOKEN" .env | cut -d= -f2)
    
    echo -e "${CYAN}📊 Service Access URLs:${NC}"
    echo "   • RStudio Server:    http://localhost:8787"
    echo "      - Username: rstudio"
    echo "      - Password: aipoch"
    echo ""
    echo "   • JupyterLab:        http://localhost:8888"
    echo "      - Token: ${jupyter_token:-aipoch}"
    echo ""
    echo "   • Opencode Server:   http://localhost:4096"
    echo ""
    
    if [ ${#INSTALLED_TOOLS[@]} -gt 0 ]; then
        echo -e "${CYAN}🆕 Tools installed this time:${NC}"
        for tool in "${INSTALLED_TOOLS[@]}"; do
            echo "   • $tool"
        done
        echo ""
    fi
    
    echo -e "${CYAN}🚀 Quick Start:${NC}"
    echo "   1. Ensure OpenClaw gateway is running:"
    echo "      openclaw gateway start"
    echo ""
    echo "   2. Send in Feishu/WhatsApp/Slack:"
    echo "      'analyze data with opencode'"
    echo ""
    echo "   3. Or use directly in browser:"
    echo "      http://localhost:8787 (RStudio)"
    echo "      http://localhost:8888 (JupyterLab)"
    echo ""
    
    echo -e "${CYAN}📖 Common Commands:${NC}"
    echo "   View logs:    $DOCKER_COMPOSE logs -f"
    echo "   Stop services:    $DOCKER_COMPOSE down"
    echo "   Restart services:    $DOCKER_COMPOSE restart"
    echo "   View status:    docker ps"
    echo ""
    
    echo -e "${CYAN}📁 Data Directories:${NC}"
    echo "   • Input data: ./data"
    echo "   • Analysis output: ./outputs"
    echo "   • Log files: ./logs"
    echo ""
    
    echo -e "${YELLOW}⚠️  Important Reminders:${NC}"
    echo "   • Configure API access: ~/.opencode/bin/opencode auth login"
    echo "   • Default passwords are for development only, change for production"
    echo "   • View full logs: $LOG_FILE"
    echo "   • If commands not found, reload shell config or restart terminal"
    echo ""
    
    echo -e "${GREEN}🎉 Installation complete! Start your medical research journey!${NC}"
    echo ""
}

main() {
    show_banner
    
    # Parse arguments
    local offline_mode=false
    local skip_build=false
    local auto_install=true
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --offline)
                offline_mode=true
                shift
                ;;
            --skip-build)
                skip_build=true
                shift
                ;;
            --no-auto-install)
                auto_install=false
                shift
                ;;
            --help|-h)
                echo "Usage: $0 [options]"
                echo ""
                echo "Options:"
                echo "  --offline           Offline mode (skip network dependency checks)"
                echo "  --skip-build        Skip Docker image build"
                echo "  --no-auto-install   Do not auto-install OpenClaw and Opencode"
                echo "  --help, -h          Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0                              # Full installation"
                echo "  $0 --skip-build                 # Skip Docker build"
                echo "  $0 --no-auto-install            # Install Docker services only"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                exit 1
                ;;
        esac
    done
    
    log "Starting Bioclaw setup v${SCRIPT_VERSION}"
    log "Arguments: offline=${offline_mode}, skip_build=${skip_build}, auto_install=${auto_install}"
    
    # Pre-flight checks
    check_in_project_root
    check_prerequisites
    check_dependencies
    
    if [ "$offline_mode" = false ]; then
        check_ports
        check_disk_space
    fi
    
    # Install tools (OpenClaw, Opencode)
    if [ "$auto_install" = true ]; then
        install_openclaw
        install_opencode
        start_opencode_server
    fi
    
    # Setup Docker environment
    init_git_submodule
    create_directories
    create_env_file
    
    if [ "$skip_build" = false ]; then
        build_docker_images
    fi
    
    start_docker_services
    
    # Initialize Bioclaw
    if [ "$auto_install" = true ]; then
        run_aipoch_init
    fi
    
    # Health check
    if health_check; then
        show_summary
        log "Setup completed successfully"
    else
        print_warning "Some services may need more time to start"
        show_summary
        log "Setup completed with warnings"
    fi
}

# Run main
main "$@"
