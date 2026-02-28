<p align="center">
  <img src="assets/logo.svg" alt="Bioclaw Logo" width="200">
</p>

<h1 align="center">Bioclaw</h1>

<p align="center">
  <strong>Open-Source Bio-Research Integration Package</strong><br>
  One-click setup for OpenClaw + Opencode + Docker Toolkit
</p>

<p align="center">
  <a href="README_CN.md">中文</a> | <strong>English</strong>
</p>

---

## Two-Step Installation

### Step 1: Install Docker

**macOS:**
1. Visit https://docs.docker.com/desktop/install/mac-install/
2. Download and install Docker Desktop
3. Open Docker Desktop and wait for "Docker Desktop is running"

**Ubuntu:**
```bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo systemctl start docker
sudo usermod -aG docker $USER
# Log out and log back in
```

### Step 2: Install Bioclaw

```bash
curl -fsSL https://raw.githubusercontent.com/Rowtion/Bioclaw/main/install.sh | bash
```

**Installation takes about 5-10 minutes** and will automatically:
- Download Bioclaw code
- Build Docker images
- Start all services
- Create the `bioclaw` command

---

## Usage

### Start
```bash
bioclaw start
```

### Stop
```bash
bioclaw stop
```

### Check Status
```bash
bioclaw status
```

### Access Research Environment

- **RStudio** (data analysis): http://localhost:8787
- **JupyterLab** (Python coding): http://localhost:8888

**Default password:** `bioclaw`

---

## What is This?

Bioclaw is an **integration package** that helps you quickly set up:

- **OpenClaw** - AI conversation gateway (connects to Feishu/Slack/WhatsApp)
- **Opencode** - Code execution environment
- **Docker Toolkit** - Pre-installed RStudio + JupyterLab + bio-research tools

**Perfect for:**
- Users who want to quickly set up bio-research environment
- Beginners who don't want to tinker with configurations
- Teams who need standardized analysis environments

---

## Project Structure

```
~/.bioclaw/
├── docker-compose.yml     # Docker configuration
├── setup.sh               # Installation script
├── data/                  # Data directory
├── outputs/               # Analysis results
└── scientific-skills/     # Bio-research skill library
```

---

## GUI Interface (Optional)

For users who prefer graphical interface:

```bash
# Launch GUI
python3 ~/.bioclaw/bioclaw-gui.py

# Or on macOS, double-click Bioclaw.app
```

---

## Need Help?

- **Installation Check:** `bash check.sh`
- **FAQ:** See [FAQ.md](FAQ.md)
- **GitHub:** https://github.com/Rowtion/Bioclaw
- **Issues:** https://github.com/Rowtion/Bioclaw/issues

---

<p align="center">
  Made with ❤️ for researchers worldwide
</p>
