# 🔬 Antigravity — Complete Internal Architecture & Hacking Guide

## 🏗️ Architecture Overview

Antigravity is a **forked VS Code (Electron 39.2.3)** by Google. Internal codenames:
- **Jetski** = Agent/AI system
- **Exa** = Infrastructure & protobuf types
- **Cascade** = Agentic tool/plugin system

```mermaid
graph TB
    subgraph ELECTRON["Antigravity.exe (Electron 39.2.3)"]
        MAIN["out/main.js (7.6MB) — Electron Main"]
        AGENT_UI["out/jetskiAgent/main.js (10.7MB) — Chat UI (Preact+Redux)"]
        EXT["extensions/antigravity/dist/extension.js (3MB) — Core Extension"]
    end

    subgraph LS["Language Server Daemon (Go Binary)"]
        LSEXE["language_server_windows_x64.exe (160MB)"]
        HTTP["HTTP :50711"]
        HTTPS["HTTPS :50710"]
        LSP["LSP :50714"]
    end

    subgraph GOOGLE["Google Cloud"]
        API["daily-cloudcode-pa.googleapis.com<br/>/v1internal:streamGenerateContent?alt=sse"]
    end

    subgraph STORAGE["Local Data (~/.gemini/antigravity/)"]
        DAEMON["daemon/ls_*.json — port + CSRF token"]
        CONV["conversations/*.pb — protobuf history"]
        MCP["mcp_config.json — tool servers"]
        BRAIN["brain/conv-id/ — artifacts"]
    end

    EXT -->|gRPC-Web + CSRF| HTTPS
    EXT -->|gRPC-Web| HTTP
    LSEXE -->|SSE streaming| API
    LSEXE --> DAEMON
    LSEXE --> CONV
    AGENT_UI -->|webview msgs| EXT
    EXT --> AGENT_UI
```

---

## 🔑 Path A — Live gRPC API (CONFIRMED WORKING!)

### Auth Method
```
Header: x-codeium-csrf-token: <token from daemon JSON>
```
Token from: `C:\Users\HP\.gemini\antigravity\daemon\ls_*.json`

### Confirmed Live Endpoints (HTTP :50711)

| Method | Status | Notes |
|--------|--------|-------|
| `Heartbeat` | ✅ **200 OK** | Ping the server |
| `GetStatus` | ✅ **200 OK** | Server state |
| `GetUserStatus` | ✅ **200 OK** | Auth/user info |
| `GetChatHistory` | ❌ 404 | Not exposed on HTTP |

### Protocol
- **gRPC-Web** over HTTP/HTTPS
- Encode: `[0x00][length 4 bytes big-endian][protobuf bytes]`
- Auth: `x-codeium-csrf-token` header required

### Python PoC Client: [g:\plans\ag_client.py](file:///g:/plans/ag_client.py)

---

## 📡 All 200+ Discovered API Methods

### Core Agent Control
```
StartCascade / StopCascade
SendUserCascadeMessage
RevertToCascadeStep
ResolveOutstandingSteps
StreamAgentStateUpdates
StreamReactiveUpdates
UpdateCascadeMemory
UpdateCascadeTrajectorySummaries
```

### Chat & UI
```
SendActionToChatPanel
StartChatClientRequestStream
RecordChatFeedback
RecordChatPanelSession
SmartFocusConversation
UpdateDetailedViewWithCascadeInput
```

### File System
```
ReadFile / WriteFile / StatUri
SaveDocument / OpenDocument
WriteCascadeEdit
ReadNotebook / WriteNotebook
SaveMediaAsArtifact
```

### Terminal & Code
```
RunTool
RunExtensionCode
SendTerminalInput / ReadTerminal / ShowTerminal / OpenTerminal
TerminateCommand
SignalExecutableIdle
StreamTerminalShellCommand
```

### Browser
```
SmartOpenBrowser
SkipBrowserSubagent
SetBrowserOpenConversation
```

### System
```
OpenUrl
PlaySound
StartAudioRecording / StopAudioRecording
StartScreenRecording / SaveScreenRecording
SimulateSegFault (lol)
```

### MCP & Extensions
```
RefreshMcpServers
RunExtensionCode
ReconnectExtensionServer
RegisterGdmUser
```

### Config & Settings
```
GetUserSettings / SetUserSettings
RefreshCustomization
UpdateCustomizationPathsFile
SetIndexConfig
SetWorkingDirectories
GetExperiments / UpdateDevExperiments
```

---

## 🎨 UI Editing Guide — Chat Interface

### How the Chat UI Works
```
cascade-panel.html          ← Bare HTML shell (just a div#react-app)
  └── panel/chat/*.js       ← Compiled webpack chunk
       └── AntigravityPanelManager  ← Root React component from @exa/chat-client
```

### What IS Editable

**1. cascade-panel.html** — Add custom elements, styles, scripts:
```html
<!-- C:\...\extensions\antigravity\cascade-panel.html -->
<!doctype html>
<html>
  <head>
    <!-- ADD: Custom CSS here -->
    <style>
      .react-app-container { background: #0a0a0a !important; }
      /* Add TTS button, custom themes, etc */
    </style>
  </head>
  <body style="margin: 0">
    <div id="react-app" class="react-app-container"></div>
    <!-- ADD: Custom JS injection -->
    <script>
      // Intercept webview messages to add TTS, etc.
      window.addEventListener('message', e => {
        if (e.data?.type === 'response') {
          // Speak response using Web Speech API!
          speechSynthesis.speak(new SpeechSynthesisUtterance(e.data.text));
        }
      });
    </script>
  </body>
</html>
```

**2. jetskiAgent/main.css** — All the styling:
```
C:\...\resources\app\out\jetskiAgent\main.css (93KB)
```
Directly editable! CSS variables, colors, fonts — all here.

**3. jetskiAgent/main.tailwind.css** — Tailwind base styles:
```
C:\...\resources\app\out\jetskiMain.tailwind.css (99KB)
```

### Adding Voice (TTS) Button — inject via cascade-panel.html

```html
<script>
// Observe DOM for new AI messages → speak them
const observer = new MutationObserver(mutations => {
  mutations.forEach(m => {
    m.addedNodes.forEach(node => {
      if (node.textContent && node.textContent.length > 10) {
        // Could trigger TTS here
      }
    });
  });
});
observer.observe(document.body, { childList: true, subtree: true });
</script>
```

---

## 🗺️ File Map (All Editable Locations)

| File | What | Risk |
|------|------|------|
| `~/.gemini/antigravity/mcp_config.json` | Add MCP tools | ✅ Zero |
| `~/.gemini/antigravity/browserAllowlist.txt` | Allow websites | ✅ Zero |
| `g:\plans\_agent/rules/*.md` | Agent behavior | ✅ Zero |
| `g:\plans\_agent/workflows/*.md` | Automation | ✅ Zero |
| [extensions/antigravity/cascade-panel.html](file:///C:/Users/HP/AppData/Local/Programs/Antigravity/resources/app/extensions/antigravity/cascade-panel.html) | Chat UI shell | ⚠️ Low |
| [out/jetskiAgent/main.css](file:///C:/Users/HP/AppData/Local/Programs/Antigravity/resources/app/out/jetskiAgent/main.css) | All CSS styles | ⚠️ Low |
| [out/jetskiMain.tailwind.css](file:///C:/Users/HP/AppData/Local/Programs/Antigravity/resources/app/out/jetskiMain.tailwind.css) | Tailwind base | ⚠️ Low |
| [extensions/antigravity/dist/extension.js](file:///C:/Users/HP/AppData/Local/Programs/Antigravity/resources/app/extensions/antigravity/dist/extension.js) | Core logic | ❌ High (minified) |
| [out/main.js](file:///C:/Users/HP/AppData/Local/Programs/Antigravity/resources/app/out/main.js) | Electron main | ❌ Very High |
| [bin/language_server_windows_x64.exe](file:///C:/Users/HP/AppData/Local/Programs/Antigravity/resources/app/extensions/antigravity/bin/language_server_windows_x64.exe) | AI brain | ❌ Binary |

---

## 🚀 Next Steps — Build The Full External Agent

### Step 1: Decode gRPC Responses
```bash
pip install protobuf grpcio grpcio-tools
```
Run [g:\plans\ag_client.py](file:///g:/plans/ag_client.py) to see raw protobuf bytes from `GetStatus`

### Step 2: Reconstruct Proto Schema
From extension.js, we know all message type names. Use `protoc` to generate Python classes.

### Step 3: Send Real Chat Messages
`SendUserCascadeMessage` — encode message string as proto field 1 → send to HTTP:50711

### Step 4: Subscribe to Response Stream
`StreamAgentStateUpdates` — subscribe to agent state changes → your external agent gets all responses!

### Step 5: Control Antigravity Programmatically
You'll have full control: start tasks, read responses, inject files, open URLs — **all from Python!**

---

## 💡 Google's API Endpoint (Bonus Find!)

All AI calls go to:
```
POST https://daily-cloudcode-pa.googleapis.com/v1internal:streamGenerateContent?alt=sse
```
This is Google's internal CodeAssist API (not public Gemini API). Changes depending on selected model.
